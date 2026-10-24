//! Pure Rust PLINK 1.9 binary BED format decoder and streaming packer.
//!
//! PLINK 1.9 BED files store genotypes in 2-bit packed format (SNP-major mode).
//! Magic header (3 bytes): 0x6C, 0x1B, 0x01
//! Genotype bit encoding (least significant bits first):
//!   00 -> Hom 1 (Homozygote A1/A1)
//!   01 -> Missing
//!   10 -> Het (Heterozygote A1/A2)
//!   11 -> Hom 2 (Homozygote A2/A2)

pub const BED_MAGIC: [u8; 3] = [0x6C, 0x1B, 0x01];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Genotype {
    Hom1 = 0,    // 00 -> A1/A1 (dosage of A1 = 2)
    Missing = 1, // 01 -> Missing genotype
    Het = 2,     // 10 -> A1/A2 (dosage of A1 = 1)
    Hom2 = 3,    // 11 -> A2/A2 (dosage of A1 = 0)
}

impl Genotype {
    #[inline]
    pub fn from_bits(bits: u8) -> Self {
        match bits & 0x03 {
            0b00 => Genotype::Hom1,
            0b01 => Genotype::Missing,
            0b10 => Genotype::Het,
            _    => Genotype::Hom2,
        }
    }

    #[inline]
    pub fn to_bits(self) -> u8 {
        match self {
            Genotype::Hom1 => 0b00,
            Genotype::Missing => 0b01,
            Genotype::Het => 0b10,
            Genotype::Hom2 => 0b11,
        }
    }
}

pub const HAS_MISSING: [bool; 256] = {
    let mut table = [false; 256];
    let mut i = 0;
    while i < 256 {
        let b = i as u8;
        if (b & 0x03) == 0x01
            || ((b >> 2) & 0x03) == 0x01
            || ((b >> 4) & 0x03) == 0x01
            || ((b >> 6) & 0x03) == 0x01
        {
            table[i] = true;
        }
        i += 1;
    }
    table
};

pub const MISSING_COUNT: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        let b = i as u8;
        let c = ((b & 0x03) == 0x01) as u8
            + (((b >> 2) & 0x03) == 0x01) as u8
            + (((b >> 4) & 0x03) == 0x01) as u8
            + (((b >> 6) & 0x03) == 0x01) as u8;
        table[i] = c;
        i += 1;
    }
    table
};

pub const A1_DOSAGE: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        let b = i as u8;
        let mut dosage = 0u8;
        let mut r = 0;
        while r < 4 {
            let bits = (b >> (r * 2)) & 0x03;
            if bits == 0b00 {
                dosage += 2;
            } else if bits == 0b10 {
                dosage += 1;
            }
            r += 1;
        }
        table[i] = dosage;
        i += 1;
    }
    table
};

pub const NON_MISSING: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = 4 - MISSING_COUNT[i];
        i += 1;
    }
    table
};

/// Calculate the number of bytes per variant for a given sample count in SNP-major mode.
#[inline]
pub fn bytes_per_variant(n_samples: usize) -> usize {
    (n_samples + 3) / 4
}

/// Fast scanner for a single variant row of raw BED bytes.
///
/// Iterates over all `n_samples` genotypes without heap allocation.
/// Excess padding bits in the final byte are safely ignored.
#[inline]
pub fn scan_row<F>(raw: &[u8], n_samples: usize, mut f: F)
where
    F: FnMut(usize, Genotype),
{
    let full_bytes = n_samples / 4;
    let remainder = n_samples % 4;

    for byte_idx in 0..full_bytes {
        let b = raw[byte_idx];
        f(byte_idx * 4, Genotype::from_bits(b));
        f(byte_idx * 4 + 1, Genotype::from_bits(b >> 2));
        f(byte_idx * 4 + 2, Genotype::from_bits(b >> 4));
        f(byte_idx * 4 + 3, Genotype::from_bits(b >> 6));
    }

    if remainder > 0 && full_bytes < raw.len() {
        let b = raw[full_bytes];
        for r in 0..remainder {
            f(full_bytes * 4 + r, Genotype::from_bits(b >> (r * 2)));
        }
    }
}

/// Pack genotypes for retained samples into a new 2-bit packed BED byte vector.
pub fn pack_filtered_row(
    raw: &[u8],
    n_samples: usize,
    sample_keep: &[bool],
    out_buf: &mut Vec<u8>,
) {
    out_buf.clear();
    let n_kept = sample_keep.iter().take(n_samples).filter(|&&k| k).count();
    let expected_bytes = bytes_per_variant(n_kept);
    out_buf.resize(expected_bytes, 0);

    let mut out_idx = 0;
    let mut out_shift = 0;
    let mut cur_byte = 0u8;

    for (s_idx, &kept) in sample_keep.iter().take(n_samples).enumerate() {
        if kept {
            let in_byte = raw[s_idx / 4];
            let in_shift = (s_idx % 4) * 2;
            let bits = (in_byte >> in_shift) & 0x03;

            cur_byte |= bits << out_shift;
            out_shift += 2;
            if out_shift == 8 {
                out_buf[out_idx] = cur_byte;
                out_idx += 1;
                cur_byte = 0;
                out_shift = 0;
            }
        }
    }
    if out_shift > 0 && out_idx < out_buf.len() {
        out_buf[out_idx] = cur_byte;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bytes_per_variant() {
        assert_eq!(bytes_per_variant(0), 0);
        assert_eq!(bytes_per_variant(1), 1);
        assert_eq!(bytes_per_variant(4), 1);
        assert_eq!(bytes_per_variant(5), 2);
        assert_eq!(bytes_per_variant(8), 2);
    }

    #[test]
    fn test_genotype_encoding() {
        assert_eq!(Genotype::from_bits(0b00), Genotype::Hom1);
        assert_eq!(Genotype::from_bits(0b01), Genotype::Missing);
        assert_eq!(Genotype::from_bits(0b10), Genotype::Het);
        assert_eq!(Genotype::from_bits(0b11), Genotype::Hom2);

        assert_eq!(Genotype::Hom1.to_bits(), 0b00);
        assert_eq!(Genotype::Missing.to_bits(), 0b01);
        assert_eq!(Genotype::Het.to_bits(), 0b10);
        assert_eq!(Genotype::Hom2.to_bits(), 0b11);
    }

    #[test]
    fn test_scan_and_pack_roundtrip() {
        // 5 samples:
        // sample 0: Hom1 (00)
        // sample 1: Missing (01)
        // sample 2: Het (10)
        // sample 3: Hom2 (11)
        // byte 0 = 00 | (01<<2) | (10<<4) | (11<<6) = 0b11100100 = 0xE4
        // sample 4: Hom1 (00)
        // byte 1 = 00 = 0x00
        let raw = vec![0xE4, 0x00];
        let mut genotypes = Vec::new();
        scan_row(&raw, 5, |_idx, g| genotypes.push(g));

        assert_eq!(genotypes, vec![
            Genotype::Hom1,
            Genotype::Missing,
            Genotype::Het,
            Genotype::Hom2,
            Genotype::Hom1,
        ]);

        // Keep samples 0, 2, 3 (Hom1, Het, Hom2) -> count = 3
        let sample_keep = vec![true, false, true, true, false];
        let mut packed = Vec::new();
        pack_filtered_row(&raw, 5, &sample_keep, &mut packed);

        // Expected packed:
        // s0: Hom1 (00)
        // s2: Het (10)
        // s3: Hom2 (11)
        // byte 0 = 00 | (10<<2) | (11<<4) = 0b00111000 = 0x38
        assert_eq!(packed, vec![0x38]);

        let mut filtered_genotypes = Vec::new();
        scan_row(&packed, 3, |_idx, g| filtered_genotypes.push(g));
        assert_eq!(filtered_genotypes, vec![
            Genotype::Hom1,
            Genotype::Het,
            Genotype::Hom2,
        ]);
    }
}
