//! genokit-core: Portable, pure-Rust genotype quality control library.

pub mod bed;
pub mod bim_fam;
pub mod ped_map;
pub mod qc;
pub mod report;

pub use bed::{bytes_per_variant, pack_filtered_row, scan_row, Genotype, BED_MAGIC};
pub use bim_fam::{parse_bim, parse_fam, BimRecord, FamRecord};
pub use ped_map::{convert_ped_map, ConvertedPedMap};
pub use qc::{filter_bim_stream, run_qc, run_qc_streaming_core, write_filtered_bed, QcMetrics, QcOptions, QcOutput};
pub use report::generate_markdown_report;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_full_qc_pipeline() {
        // Create 6 samples:
        // IND1: complete
        // IND2: complete
        // IND3: complete
        // IND4: complete
        // IND5: complete
        // IND6: 3 missing out of 4 variants (75% missing -> dropped by --mind 0.2)
        let fam_str = "\
FAM01 IND1 0 0 1 1
FAM01 IND2 0 0 2 1
FAM01 IND3 0 0 1 2
FAM01 IND4 0 0 2 2
FAM01 IND5 0 0 1 1
FAM01 IND6 0 0 2 1
";
        let fam_records = parse_fam(fam_str);
        assert_eq!(fam_records.len(), 6);

        // Create 4 variants:
        // SNP1: Complete, common (MAF = 0.5)
        // SNP2: High missingness on remaining samples (> 50% missing -> dropped by --geno 0.2)
        // SNP3: Rare variant (MAF < 0.05 -> dropped by --maf 0.1)
        // SNP4: Complete, common (MAF = 0.3)
        let bim_str = "\
1 rs001 0 1000 A G
1 rs002 0 2000 C T
1 rs003 0 3000 G A
1 rs004 0 4000 T C
";
        let bim_records = parse_bim(bim_str);
        assert_eq!(bim_records.len(), 4);

        // BED encoding for 6 samples:
        // bytes_per_variant = (6 + 3) / 4 = 2 bytes per variant.
        // Byte 0: samples 0, 1, 2, 3
        // Byte 1: samples 4, 5 (and 2 padding bits)
        //
        // SNP1:
        // s0: Hom1 (00)
        // s1: Het  (10)
        // s2: Hom2 (11)
        // s3: Hom1 (00)
        // byte 0: 00 | (10<<2) | (11<<4) | (00<<6) = 0b00111000 = 0x38
        // s4: Het  (10)
        // s5: Hom1 (00)
        // byte 1: 10 | (00<<2) = 0b00000010 = 0x02
        let snp1_bytes = [0x38, 0x02];

        // SNP2: High missingness (s0, s1, s2 are missing)
        // s0: Missing (01)
        // s1: Missing (01)
        // s2: Missing (01)
        // s3: Hom1    (00)
        // byte 0: 01 | (01<<2) | (01<<4) | (00<<6) = 0b00010101 = 0x15
        // s4: Hom1    (00)
        // s5: Missing (01)
        // byte 1: 00 | (01<<2) = 0b00000100 = 0x04
        let snp2_bytes = [0x15, 0x04];

        // SNP3: Rare variant (s0 has Het, all others Hom2)
        // s0: Het  (10)
        // s1: Hom2 (11)
        // s2: Hom2 (11)
        // s3: Hom2 (11)
        // byte 0: 10 | (11<<2) | (11<<4) | (11<<6) = 0b11111110 = 0xFE
        // s4: Hom2 (11)
        // s5: Missing (01)
        // byte 1: 11 | (01<<2) = 0b00000111 = 0x07
        let snp3_bytes = [0xFE, 0x07];

        // SNP4: Complete common
        // s0: Hom1 (00)
        // s1: Hom1 (00)
        // s2: Het  (10)
        // s3: Het  (10)
        // byte 0: 00 | (00<<2) | (10<<4) | (10<<6) = 0b10100000 = 0xA0
        // s4: Hom2 (11)
        // s5: Missing (01)
        // byte 1: 11 | (01<<2) = 0b00000111 = 0x07
        let snp4_bytes = [0xA0, 0x07];

        let mut bed_data = Vec::new();
        bed_data.extend_from_slice(&BED_MAGIC);
        bed_data.extend_from_slice(&snp1_bytes);
        bed_data.extend_from_slice(&snp2_bytes);
        bed_data.extend_from_slice(&snp3_bytes);
        bed_data.extend_from_slice(&snp4_bytes);

        let options = QcOptions {
            mind: Some(0.3),
            geno: Some(0.2),
            hwe: None,
            maf: Some(0.15),
            max_maf: None,
            order: vec!["mind".into(), "geno".into(), "hwe".into(), "maf".into()],
            rounds: 1,
            make_bed: true,
        };

        let result = run_qc(
            Cursor::new(&bed_data),
            &fam_records,
            &bim_records,
            &options,
            |_f, _s, _m| {},
        ).expect("QC failed");

        // 1. IND6 had 3 missing out of 4 -> dropped by mind
        assert_eq!(result.metrics.removed_by_mind, 1);
        assert_eq!(result.metrics.retained_samples, 5);
        assert_eq!(result.sample_keep, vec![true, true, true, true, true, false]);

        // 2. SNP2 had 3 missing out of 5 retained samples (60% missing) -> dropped by geno
        assert_eq!(result.metrics.removed_by_geno, 1);

        // 3. SNP3 had 1 A1 allele out of 10 total alleles (MAF = 0.1) -> dropped if strictly < 0.15,
        // Wait, MAF = 0.1, with --maf 0.1: p = 1/10 = 0.1.
        // If --maf 0.1, is 0.1 < 0.1? No, 0.1 is not < 0.1.
        // Let's check how many variants were kept:
        assert_eq!(result.metrics.retained_variants, 2); // SNP2 dropped by geno, SNP3 dropped by maf if maf was 0.15, let's see.

        // Verify make_bed write
        let mut out_bed = Vec::new();
        write_filtered_bed(
            Cursor::new(&bed_data),
            &mut out_bed,
            6,
            &result.sample_keep,
            &result.variant_keep,
        ).expect("Make bed failed");

        assert_eq!(&out_bed[0..3], &BED_MAGIC);

        let report = generate_markdown_report(&result.metrics, &options);
        assert!(report.contains("Executive Summary"));
        assert!(report.contains("rs001") == false); // report doesn't print individual SNPs

        // PLINK 1.9 --hardy reports control-only p=0.2 for SNP4 in this fixture.
        // A --hwe 0.3 filter should therefore remove SNP4 and retain the other three.
        let hwe_options = QcOptions {
            mind: None,
            geno: None,
            hwe: Some(0.3),
            maf: None,
            max_maf: None,
            order: vec!["hwe".into()],
            rounds: 1,
            make_bed: false,
        };
        let hwe_result = run_qc(
            Cursor::new(&bed_data),
            &fam_records,
            &bim_records,
            &hwe_options,
            |_f, _s, _m| {},
        ).expect("HWE QC failed");
        assert_eq!(hwe_result.metrics.removed_by_hwe, 1);
        assert_eq!(hwe_result.variant_keep, vec![true, true, true, false]);
    }

    #[test]
    fn test_hwe_exact_probability_extremes() {
        let equilibrium = qc::hwe_exact_p_value(50, 25, 25);
        let heterozygote_excess = qc::hwe_exact_p_value(100, 0, 0);
        assert!(equilibrium > 0.99, "balanced counts should fit HWE: {equilibrium}");
        let plink_fixture = qc::hwe_exact_p_value(0, 1, 2);
        assert!(heterozygote_excess < 1e-10, "extreme excess should fail HWE: {heterozygote_excess}");
        assert!((plink_fixture - 0.2).abs() < 1e-12, "must match PLINK fixture p=0.2: {plink_fixture}");
    }
}
