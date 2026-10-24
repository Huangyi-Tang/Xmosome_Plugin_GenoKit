//! Pure-Rust parser and converter for PLINK text PED and MAP formats to binary BED, BIM, and FAM.

use std::collections::HashMap;
use std::io::{BufRead, Seek, SeekFrom};
use crate::bed::{bytes_per_variant, BED_MAGIC};
use crate::bim_fam::{BimRecord, FamRecord};

pub struct ConvertedPedMap {
    pub fam_records: Vec<FamRecord>,
    pub bim_records: Vec<BimRecord>,
    pub bed_data: Vec<u8>,
}

#[derive(Debug, Clone)]
struct Marker {
    chr: String,
    id: String,
    cm: f64,
    bp: i64,
}

/// Convert PLINK .map and .ped streams into in-memory .bed, .bim, and .fam structures.
pub fn convert_ped_map<R1, R2>(
    map_reader: R1,
    mut ped_reader: R2,
) -> Result<ConvertedPedMap, String>
where
    R1: BufRead,
    R2: BufRead + Seek,
{
    // 1. Parse MAP file
    let mut markers = Vec::new();
    for (line_idx, line) in map_reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading MAP line {}: {e}", line_idx + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut tokens = trimmed.split_whitespace();
        let chr = tokens.next().ok_or_else(|| format!("Missing chr at MAP line {}", line_idx + 1))?.to_string();
        let id = tokens.next().ok_or_else(|| format!("Missing ID at MAP line {}", line_idx + 1))?.to_string();
        let cm: f64 = tokens
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let bp: i64 = tokens
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        markers.push(Marker { chr, id, cm, bp });
    }

    let n_variants = markers.len();
    if n_variants == 0 {
        return Err("No markers found in MAP file".into());
    }

    // 2. Pass 1 over PED: collect sample headers and count allele frequencies per variant
    let mut allele_counts: Vec<HashMap<u8, usize>> = vec![HashMap::new(); n_variants];
    let mut fam_records = Vec::new();

    ped_reader.seek(SeekFrom::Start(0)).map_err(|e| e.to_string())?;

    for (s_idx, line) in ped_reader.by_ref().lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading PED line {}: {e}", s_idx + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut tokens = trimmed.split_whitespace();
        let fid = tokens.next().ok_or_else(|| format!("Missing family_id at PED line {}", s_idx + 1))?.to_string();
        let iid = tokens.next().ok_or_else(|| format!("Missing individual_id at PED line {}", s_idx + 1))?.to_string();
        let pid = tokens.next().unwrap_or("0").to_string();
        let mid = tokens.next().unwrap_or("0").to_string();
        let sex = tokens.next().unwrap_or("0").to_string();
        let pheno = tokens.next().unwrap_or("-9").to_string();

        fam_records.push(FamRecord {
            fid,
            iid,
            pid,
            mid,
            sex,
            pheno,
        });

        for v_idx in 0..n_variants {
            let a1 = tokens.next();
            let a2 = tokens.next();
            if let Some(tok1) = a1 {
                if let Some(&b) = tok1.as_bytes().first() {
                    if b != b'0' && b != b'N' && b != b'?' {
                        *allele_counts[v_idx].entry(b.to_ascii_uppercase()).or_insert(0) += 1;
                    }
                }
            }
            if let Some(tok2) = a2 {
                if let Some(&b) = tok2.as_bytes().first() {
                    if b != b'0' && b != b'N' && b != b'?' {
                        *allele_counts[v_idx].entry(b.to_ascii_uppercase()).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    let n_samples = fam_records.len();
    if n_samples == 0 {
        return Err("No samples found in PED file".into());
    }

    // 3. Determine A1 (minor) and A2 (major) alleles and construct BIM records
    let mut variant_alleles: Vec<(u8, u8)> = Vec::with_capacity(n_variants);
    let mut bim_records = Vec::with_capacity(n_variants);

    for (v_idx, marker) in markers.into_iter().enumerate() {
        let counts = &allele_counts[v_idx];
        let mut sorted: Vec<(u8, usize)> = counts.iter().map(|(&b, &c)| (b, c)).collect();
        // Sort by frequency ascending (minor first), then alphabetical
        sorted.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));

        let (a1, a2) = match sorted.len() {
            0 => (b'0', b'0'),
            1 => (sorted[0].0, b'0'),
            _ => (sorted[0].0, sorted[sorted.len() - 1].0),
        };

        variant_alleles.push((a1, a2));
        bim_records.push(BimRecord {
            chr: marker.chr,
            id: marker.id,
            cm: marker.cm.to_string(),
            pos: marker.bp.max(0) as u64,
            a1: (a1 as char).to_string(),
            a2: (a2 as char).to_string(),
        });
    }

    // 4. Pass 2 over PED: pack genotypes into BED buffer
    let bytes_per_snp = bytes_per_variant(n_samples);
    let total_bed_size = 3 + n_variants * bytes_per_snp;
    let mut bed_data = vec![0x55u8; total_bed_size]; // Default to 01 (Missing) for all samples
    bed_data[0..3].copy_from_slice(&BED_MAGIC);

    ped_reader.seek(SeekFrom::Start(0)).map_err(|e| e.to_string())?;

    for (s_idx, line) in ped_reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading PED line {} in pass 2: {e}", s_idx + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut tokens = trimmed.split_whitespace();
        // Skip first 6 columns
        for _ in 0..6 {
            tokens.next();
        }

        let shift = (s_idx % 4) * 2;
        let s_byte = s_idx / 4;

        for v_idx in 0..n_variants {
            let a1_tok = tokens.next();
            let a2_tok = tokens.next();

            let a1 = a1_tok.and_then(|t| t.as_bytes().first().copied()).map(|b| b.to_ascii_uppercase());
            let a2 = a2_tok.and_then(|t| t.as_bytes().first().copied()).map(|b| b.to_ascii_uppercase());

            let (ref_a1, ref_a2) = variant_alleles[v_idx];

            let bits = match (a1, a2) {
                (Some(a), Some(b)) if a != b'0' && a != b'N' && a != b'?' && b != b'0' && b != b'N' && b != b'?' => {
                    if a == ref_a1 && b == ref_a1 {
                        0b00u8 // Hom1
                    } else if a == ref_a2 && b == ref_a2 {
                        0b11u8 // Hom2
                    } else if (a == ref_a1 && b == ref_a2) || (a == ref_a2 && b == ref_a1) {
                        0b10u8 // Het
                    } else {
                        0b01u8 // Missing (unrecognized allele)
                    }
                }
                _ => 0b01u8, // Missing
            };

            let byte_pos = 3 + v_idx * bytes_per_snp + s_byte;
            bed_data[byte_pos] &= !(0x03 << shift);
            bed_data[byte_pos] |= bits << shift;
        }
    }

    Ok(ConvertedPedMap {
        fam_records,
        bim_records,
        bed_data,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_ped_map_conversion() {
        let map_data = "\
1 rs1 0 1000
1 rs2 0 2000
";
        let ped_data = "\
FAM01 IND1 0 0 1 1 A A C C
FAM01 IND2 0 0 2 1 A G C T
FAM01 IND3 0 0 1 2 G G T T
";
        let converted = convert_ped_map(
            Cursor::new(map_data),
            Cursor::new(ped_data),
        ).expect("conversion should succeed");

        assert_eq!(converted.fam_records.len(), 3);
        assert_eq!(converted.bim_records.len(), 2);
        assert_eq!(converted.bim_records[0].id, "rs1");
        assert_eq!(converted.bim_records[1].id, "rs2");

        // BED data length: 3 header bytes + 2 variants * 1 byte per variant = 5 bytes
        assert_eq!(converted.bed_data.len(), 5);
        assert_eq!(&converted.bed_data[0..3], &BED_MAGIC);
    }
}
