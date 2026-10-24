//! Core streaming quality control algorithms for genotype data (PLINK-aligned).

use std::io::{BufRead, Read, Seek, SeekFrom, Write};
use serde::{Deserialize, Serialize};
use crate::bed::{
    bytes_per_variant, pack_filtered_row, A1_DOSAGE, BED_MAGIC, HAS_MISSING, MISSING_COUNT,
    NON_MISSING,
};
use crate::bim_fam::{BimRecord, FamRecord};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcOptions {
    /// Sample missingness threshold (`--mind`). `None` disables the filter.
    pub mind: Option<f64>,
    /// Variant missingness threshold (`--geno`). `None` disables the filter.
    pub geno: Option<f64>,
    /// Founder/control HWE exact-test p-value threshold (`--hwe`). `None` disables.
    pub hwe: Option<f64>,
    /// Minimum MAF (`--maf`). `None` disables the lower bound.
    pub maf: Option<f64>,
    /// Maximum MAF (`--max-maf`). `None` disables the upper bound.
    pub max_maf: Option<f64>,
    pub order: Vec<String>,
    pub rounds: usize,
    pub make_bed: bool,
}

impl Default for QcOptions {
    fn default() -> Self {
        Self {
            mind: Some(0.1),
            geno: Some(0.1),
            hwe: None,
            maf: Some(0.01),
            max_maf: None,
            order: vec!["mind".into(), "geno".into(), "hwe".into(), "maf".into()],
            rounds: 1,
            make_bed: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcMetrics {
    pub initial_samples: usize,
    pub retained_samples: usize,
    pub removed_by_mind: usize,

    pub initial_variants: usize,
    pub retained_variants: usize,
    pub removed_by_geno: usize,
    pub removed_by_hwe: usize,
    pub removed_by_maf: usize,

    pub rounds_executed: usize,
    pub mean_sample_missingness: f64,
    pub mean_variant_missingness: f64,
    pub mean_maf: f64,
}

pub struct QcOutput {
    pub metrics: QcMetrics,
    pub sample_keep: Vec<bool>,
    pub variant_keep: Vec<bool>,
    pub filtered_fam: String,
    pub filtered_bim: String,
}

/// Run PLINK-aligned streaming genotype QC over any Seekable Read stream.
pub fn run_qc<R, F>(
    mut bed_reader: R,
    fam_records: &[FamRecord],
    bim_records: &[BimRecord],
    options: &QcOptions,
    progress: F,
) -> Result<QcOutput, String>
where
    R: Read + Seek,
    F: FnMut(f64, &str, &str),
{
    let _n_samples = fam_records.len();
    let n_variants = bim_records.len();

    let (metrics, sample_keep, variant_keep, filtered_fam) = run_qc_streaming_core(
        &mut bed_reader,
        fam_records,
        n_variants,
        options,
        progress,
    )?;

    let mut filtered_bim_lines = Vec::new();
    for (j, kept) in variant_keep.iter().enumerate() {
        if *kept && j < bim_records.len() {
            filtered_bim_lines.push(bim_records[j].to_bim_line());
        }
    }
    let filtered_bim = filtered_bim_lines.join("\n") + "\n";

    Ok(QcOutput {
        metrics,
        sample_keep,
        variant_keep,
        filtered_fam,
        filtered_bim,
    })
}

/// Core streaming QC execution without requiring full in-memory BIM structures.
pub fn run_qc_streaming_core<R, F>(
    bed_reader: &mut R,
    fam_records: &[FamRecord],
    n_variants: usize,
    options: &QcOptions,
    mut progress: F,
) -> Result<(QcMetrics, Vec<bool>, Vec<bool>, String), String>
where
    R: Read + Seek,
    F: FnMut(f64, &str, &str),
{
    let n_samples = fam_records.len();

    if n_samples == 0 {
        return Err("No samples found in FAM file".into());
    }
    if n_variants == 0 {
        return Err("No variants found in BIM file".into());
    }

    progress(0.05, "bed-validation", "Verifying PLINK BED header and file integrity");

    // 1. Verify BED Magic Header
    bed_reader.seek(SeekFrom::Start(0)).map_err(|e| e.to_string())?;
    let mut magic = [0u8; 3];
    bed_reader.read_exact(&mut magic).map_err(|e| format!("Failed to read BED header: {e}"))?;
    if magic != BED_MAGIC {
        return Err(format!(
            "Invalid BED header: expected [0x6c, 0x1b, 0x01] (SNP-major mode), got [{:#04x}, {:#04x}, {:#04x}]",
            magic[0], magic[1], magic[2]
        ));
    }

    let bytes_per_snp = bytes_per_variant(n_samples);
    let mut sample_keep = vec![true; n_samples];
    let mut variant_keep = vec![true; n_variants];

    let mut removed_by_mind = 0usize;
    let mut removed_by_geno = 0usize;
    let mut removed_by_hwe = 0usize;
    let mut removed_by_maf = 0usize;

    let mut buf = vec![0u8; bytes_per_snp];
    let max_rounds = options.rounds.max(1);
    let mut rounds_executed = 0usize;

    for round in 1..=max_rounds {
        rounds_executed = round;
        let mut changes_in_round = 0usize;

        let mut order = options.order.clone();
        if options.hwe.is_some() && !order.iter().any(|step| step.eq_ignore_ascii_case("hwe")) {
            let position = order.iter().position(|step| step.eq_ignore_ascii_case("maf"))
                .unwrap_or(order.len());
            order.insert(position, "hwe".into());
        }
        for step in &order {
            match step.trim().to_lowercase().as_str() {
                "mind" => {
                    let Some(mind_thresh) = options.mind else { continue };
                    let active_variants = variant_keep.iter().filter(|&&k| k).count();
                    if active_variants == 0 {
                        break;
                    }

                    progress(
                        0.1 + (round as f64 - 1.0) * 0.2,
                        "qc-mind",
                        &format!("Round {round}: Computing sample missingness (--mind {mind_thresh})"),
                    );

                    bed_reader.seek(SeekFrom::Start(3)).map_err(|e| e.to_string())?;
                    let mut sample_missing_counts = vec![0u32; n_samples];
                    let full_bytes = n_samples / 4;
                    let remainder = n_samples % 4;

                    for (v_idx, &is_kept) in variant_keep.iter().enumerate() {
                        bed_reader.read_exact(&mut buf).map_err(|e| format!("Error reading SNP {v_idx}: {e}"))?;
                        if !is_kept {
                            continue;
                        }
                        for byte_idx in 0..full_bytes {
                            let b = buf[byte_idx];
                            if HAS_MISSING[b as usize] {
                                let base = byte_idx * 4;
                                if (b & 0x03) == 0x01 { sample_missing_counts[base] += 1; }
                                if ((b >> 2) & 0x03) == 0x01 { sample_missing_counts[base + 1] += 1; }
                                if ((b >> 4) & 0x03) == 0x01 { sample_missing_counts[base + 2] += 1; }
                                if ((b >> 6) & 0x03) == 0x01 { sample_missing_counts[base + 3] += 1; }
                            }
                        }
                        if remainder > 0 && full_bytes < buf.len() {
                            let b = buf[full_bytes];
                            let base = full_bytes * 4;
                            for r in 0..remainder {
                                if ((b >> (r * 2)) & 0x03) == 0x01 {
                                    sample_missing_counts[base + r] += 1;
                                }
                            }
                        }
                    }

                    for s_idx in 0..n_samples {
                        if sample_keep[s_idx] {
                            let rate = sample_missing_counts[s_idx] as f64 / active_variants as f64;
                            if rate > mind_thresh {
                                sample_keep[s_idx] = false;
                                removed_by_mind += 1;
                                changes_in_round += 1;
                            }
                        }
                    }
                }
                "geno" => {
                    let Some(geno_thresh) = options.geno else { continue };
                    let active_samples = sample_keep.iter().filter(|&&k| k).count();
                    if active_samples == 0 {
                        break;
                    }

                    progress(
                        0.2 + (round as f64 - 1.0) * 0.2,
                        "qc-geno",
                        &format!("Round {round}: Computing variant missingness (--geno {geno_thresh})"),
                    );

                    bed_reader.seek(SeekFrom::Start(3)).map_err(|e| e.to_string())?;
                    let full_bytes = n_samples / 4;
                    let remainder = n_samples % 4;
                    let all_samples_kept = sample_keep.iter().all(|&k| k);

                    for v_idx in 0..n_variants {
                        bed_reader.read_exact(&mut buf).map_err(|e| format!("Error reading SNP {v_idx}: {e}"))?;
                        if !variant_keep[v_idx] {
                            continue;
                        }
                        let mut var_missing = 0usize;
                        if all_samples_kept {
                            for byte_idx in 0..full_bytes {
                                var_missing += MISSING_COUNT[buf[byte_idx] as usize] as usize;
                            }
                            if remainder > 0 && full_bytes < buf.len() {
                                let b = buf[full_bytes];
                                for r in 0..remainder {
                                    if ((b >> (r * 2)) & 0x03) == 0x01 {
                                        var_missing += 1;
                                    }
                                }
                            }
                        } else {
                            for byte_idx in 0..full_bytes {
                                let b = buf[byte_idx];
                                if HAS_MISSING[b as usize] {
                                    let base = byte_idx * 4;
                                    if sample_keep[base] && (b & 0x03) == 0x01 { var_missing += 1; }
                                    if sample_keep[base + 1] && ((b >> 2) & 0x03) == 0x01 { var_missing += 1; }
                                    if sample_keep[base + 2] && ((b >> 4) & 0x03) == 0x01 { var_missing += 1; }
                                    if sample_keep[base + 3] && ((b >> 6) & 0x03) == 0x01 { var_missing += 1; }
                                }
                            }
                            if remainder > 0 && full_bytes < buf.len() {
                                let b = buf[full_bytes];
                                let base = full_bytes * 4;
                                for r in 0..remainder {
                                    if sample_keep[base + r] && ((b >> (r * 2)) & 0x03) == 0x01 {
                                        var_missing += 1;
                                    }
                                }
                            }
                        }
                        let rate = var_missing as f64 / active_samples as f64;
                        if rate > geno_thresh {
                            variant_keep[v_idx] = false;
                            removed_by_geno += 1;
                            changes_in_round += 1;
                        }
                    }
                }
                "hwe" => {
                    let Some(threshold) = options.hwe else { continue };
                    if !(0.0..=1.0).contains(&threshold) {
                        return Err(format!("--hwe must be between 0 and 1, got {threshold}"));
                    }

                    let observed_phenotypes: Vec<f64> = fam_records.iter()
                        .filter_map(|sample| sample.pheno.parse::<f64>().ok())
                        .filter(|value| *value != 0.0 && *value != -9.0)
                        .collect();
                    let has_binary_phenotype = !observed_phenotypes.is_empty()
                        && observed_phenotypes.iter().all(|value| *value == 1.0 || *value == 2.0);
                    let hwe_samples: Vec<bool> = fam_records.iter().enumerate().map(|(index, sample)| {
                        let founder = matches!(sample.pid.as_str(), "0" | "-9")
                            && matches!(sample.mid.as_str(), "0" | "-9");
                        let is_control = sample.pheno.parse::<f64>().ok() == Some(1.0);
                        sample_keep[index]
                            && founder
                            && (!has_binary_phenotype || is_control)
                    }).collect();
                    let active_hwe_samples = hwe_samples.iter().filter(|&&keep| keep).count();
                    if active_hwe_samples == 0 {
                        continue;
                    }

                    progress(
                        0.25 + (round as f64 - 1.0) * 0.2,
                        "qc-hwe",
                        &format!("Round {round}: Computing HWE exact tests (--hwe {threshold})"),
                    );

                    bed_reader.seek(SeekFrom::Start(3)).map_err(|e| e.to_string())?;
                    for v_idx in 0..n_variants {
                        bed_reader.read_exact(&mut buf).map_err(|e| format!("Error reading SNP {v_idx}: {e}"))?;
                        if !variant_keep[v_idx] {
                            continue;
                        }

                        let mut hom_a1 = 0usize;
                        let mut het = 0usize;
                        let mut hom_a2 = 0usize;
                        for sample_idx in 0..n_samples {
                            if !hwe_samples[sample_idx] {
                                continue;
                            }
                            let bits = (buf[sample_idx / 4] >> ((sample_idx % 4) * 2)) & 0x03;
                            match bits {
                                0b00 => hom_a1 += 1,
                                0b10 => het += 1,
                                0b11 => hom_a2 += 1,
                                _ => {}
                            }
                        }
                        if hom_a1 + het + hom_a2 == 0 {
                            continue;
                        }
                        let p_value = hwe_exact_p_value(het, hom_a1, hom_a2);
                        if p_value < threshold {
                            variant_keep[v_idx] = false;
                            removed_by_hwe += 1;
                            changes_in_round += 1;
                        }
                    }
                }
                "maf" => {
                    if options.maf.is_none() && options.max_maf.is_none() {
                        continue;
                    }
                    let maf_thresh = options.maf.unwrap_or(0.0);
                    let active_samples = sample_keep.iter().filter(|&&k| k).count();
                    if active_samples == 0 {
                        break;
                    }

                    progress(
                        0.3 + (round as f64 - 1.0) * 0.2,
                        "qc-maf",
                        &format!("Round {round}: Computing minor allele frequency (--maf {maf_thresh})"),
                    );

                    bed_reader.seek(SeekFrom::Start(3)).map_err(|e| e.to_string())?;
                    let full_bytes = n_samples / 4;
                    let remainder = n_samples % 4;
                    let all_samples_kept = sample_keep.iter().all(|&k| k);

                    for v_idx in 0..n_variants {
                        bed_reader.read_exact(&mut buf).map_err(|e| format!("Error reading SNP {v_idx}: {e}"))?;
                        if !variant_keep[v_idx] {
                            continue;
                        }
                        let (a1_count, non_missing) = if all_samples_kept {
                            let mut a1 = 0usize;
                            let mut nm = 0usize;
                            for byte_idx in 0..full_bytes {
                                let b = buf[byte_idx];
                                a1 += A1_DOSAGE[b as usize] as usize;
                                nm += NON_MISSING[b as usize] as usize;
                            }
                            if remainder > 0 && full_bytes < buf.len() {
                                let b = buf[full_bytes];
                                for r in 0..remainder {
                                    let bits = (b >> (r * 2)) & 0x03;
                                    match bits {
                                        0b00 => { a1 += 2; nm += 1; }
                                        0b10 => { a1 += 1; nm += 1; }
                                        0b11 => { nm += 1; }
                                        _ => {}
                                    }
                                }
                            }
                            (a1, nm)
                        } else {
                            let mut a1 = 0usize;
                            let mut nm = 0usize;
                            for byte_idx in 0..full_bytes {
                                let b = buf[byte_idx];
                                let base = byte_idx * 4;
                                for r in 0..4 {
                                    if sample_keep[base + r] {
                                        let bits = (b >> (r * 2)) & 0x03;
                                        match bits {
                                            0b00 => { a1 += 2; nm += 1; }
                                            0b10 => { a1 += 1; nm += 1; }
                                            0b11 => { nm += 1; }
                                            _ => {}
                                        }
                                    }
                                }
                            }
                            if remainder > 0 && full_bytes < buf.len() {
                                let b = buf[full_bytes];
                                let base = full_bytes * 4;
                                for r in 0..remainder {
                                    if sample_keep[base + r] {
                                        let bits = (b >> (r * 2)) & 0x03;
                                        match bits {
                                            0b00 => { a1 += 2; nm += 1; }
                                            0b10 => { a1 += 1; nm += 1; }
                                            0b11 => { nm += 1; }
                                            _ => {}
                                        }
                                    }
                                }
                            }
                            (a1, nm)
                        };

                        if non_missing == 0 {
                            variant_keep[v_idx] = false;
                            removed_by_maf += 1;
                            changes_in_round += 1;
                            continue;
                        }

                        let total_alleles = non_missing * 2;
                        let p = a1_count as f64 / total_alleles as f64;
                        let maf = p.min(1.0 - p);

                        if maf < maf_thresh || options.max_maf.map_or(false, |m| maf > m) {
                            variant_keep[v_idx] = false;
                            removed_by_maf += 1;
                            changes_in_round += 1;
                        }
                    }
                }
                _ => {}
            }
        }

        if changes_in_round == 0 {
            break;
        }
    }

    progress(0.75, "calculating-distributions", "Collecting final summary distributions");

    // Final pass to calculate retained distribution metrics
    let retained_samples = sample_keep.iter().filter(|&&k| k).count();
    let retained_variants = variant_keep.iter().filter(|&&k| k).count();

    let mut sum_var_missing_rate = 0.0f64;
    let mut sum_maf = 0.0f64;
    let mut total_missing_calls = 0usize;

    if retained_variants > 0 && retained_samples > 0 {
        bed_reader.seek(SeekFrom::Start(3)).map_err(|e| e.to_string())?;
        let full_bytes = n_samples / 4;
        let remainder = n_samples % 4;
        let all_samples_kept = sample_keep.iter().all(|&k| k);

        for (_v_idx, &is_kept) in variant_keep.iter().enumerate() {
            bed_reader.read_exact(&mut buf).map_err(|e| format!("Error reading SNP {_v_idx}: {e}"))?;
            if !is_kept {
                continue;
            }
            let (var_missing, a1_count, non_missing) = if all_samples_kept {
                let mut var_m = 0usize;
                let mut a1_c = 0usize;
                let mut non_m = 0usize;
                for byte_idx in 0..full_bytes {
                    let b = buf[byte_idx];
                    var_m += MISSING_COUNT[b as usize] as usize;
                    a1_c += A1_DOSAGE[b as usize] as usize;
                    non_m += NON_MISSING[b as usize] as usize;
                }
                if remainder > 0 && full_bytes < buf.len() {
                    let b = buf[full_bytes];
                    for r in 0..remainder {
                        let bits = (b >> (r * 2)) & 0x03;
                        if bits == 0b01 {
                            var_m += 1;
                        } else {
                            non_m += 1;
                            if bits == 0b00 { a1_c += 2; }
                            else if bits == 0b10 { a1_c += 1; }
                        }
                    }
                }
                (var_m, a1_c, non_m)
            } else {
                let mut var_m = 0usize;
                let mut a1_c = 0usize;
                let mut non_m = 0usize;
                for byte_idx in 0..full_bytes {
                    let b = buf[byte_idx];
                    let base = byte_idx * 4;
                    for r in 0..4 {
                        if sample_keep[base + r] {
                            let bits = (b >> (r * 2)) & 0x03;
                            match bits {
                                0b01 => var_m += 1,
                                0b00 => { a1_c += 2; non_m += 1; }
                                0b10 => { a1_c += 1; non_m += 1; }
                                0b11 => { non_m += 1; }
                                _ => {}
                            }
                        }
                    }
                }
                if remainder > 0 && full_bytes < buf.len() {
                    let b = buf[full_bytes];
                    let base = full_bytes * 4;
                    for r in 0..remainder {
                        if sample_keep[base + r] {
                            let bits = (b >> (r * 2)) & 0x03;
                            match bits {
                                0b01 => var_m += 1,
                                0b00 => { a1_c += 2; non_m += 1; }
                                0b10 => { a1_c += 1; non_m += 1; }
                                0b11 => { non_m += 1; }
                                _ => {}
                            }
                        }
                    }
                }
                (var_m, a1_c, non_m)
            };

            total_missing_calls += var_missing;
            sum_var_missing_rate += var_missing as f64 / retained_samples as f64;
            if non_missing > 0 {
                let p = a1_count as f64 / (2 * non_missing) as f64;
                sum_maf += p.min(1.0 - p);
            }
        }
    }

    let mean_variant_missingness = if retained_variants > 0 {
        sum_var_missing_rate / retained_variants as f64
    } else {
        0.0
    };
    let mean_maf = if retained_variants > 0 {
        sum_maf / retained_variants as f64
    } else {
        0.0
    };
    let mean_sample_missingness = if retained_samples > 0 && retained_variants > 0 {
        total_missing_calls as f64 / (retained_samples * retained_variants) as f64
    } else {
        0.0
    };

    // Build filtered FAM content
    let mut filtered_fam_lines = Vec::new();
    for (i, kept) in sample_keep.iter().enumerate() {
        if *kept {
            filtered_fam_lines.push(fam_records[i].to_fam_line());
        }
    }
    let filtered_fam = filtered_fam_lines.join("\n") + "\n";

    let metrics = QcMetrics {
        initial_samples: n_samples,
        retained_samples,
        removed_by_mind,
        initial_variants: n_variants,
        retained_variants,
        removed_by_geno,
        removed_by_hwe,
        removed_by_maf,
        rounds_executed,
        mean_sample_missingness,
        mean_variant_missingness,
        mean_maf,
    };

    Ok((metrics, sample_keep, variant_keep, filtered_fam))
}

/// Two-sided Hardy-Weinberg exact test from the Wigginton et al. recurrence
/// used by PLINK. The p-value is the sum of genotype configurations whose
/// conditional probability is no greater than the observed configuration.
pub fn hwe_exact_p_value(obs_hets: usize, obs_hom1: usize, obs_hom2: usize) -> f64 {
    let obs_hom_r = obs_hom1.min(obs_hom2);
    let obs_hom_c = obs_hom1.max(obs_hom2);
    let rare_copies = 2 * obs_hom_r + obs_hets;
    let genotypes = obs_hets + obs_hom_c + obs_hom_r;
    if genotypes == 0 {
        return 1.0;
    }

    let mut probabilities = vec![0.0f64; rare_copies + 1];
    let mut midpoint = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
    if (rare_copies & 1) != (midpoint & 1) {
        midpoint += 1;
    }

    let mut current_hets = midpoint;
    let mut current_hom_r = (rare_copies - midpoint) / 2;
    let mut current_hom_c = genotypes - current_hets - current_hom_r;
    probabilities[midpoint] = 1.0;
    let mut probability_sum = 1.0;

    while current_hets >= 2 {
        let next = probabilities[current_hets]
            * current_hets as f64
            * (current_hets - 1) as f64
            / (4.0 * (current_hom_r + 1) as f64 * (current_hom_c + 1) as f64);
        probabilities[current_hets - 2] = next;
        probability_sum += next;
        current_hets -= 2;
        current_hom_r += 1;
        current_hom_c += 1;
    }

    current_hets = midpoint;
    current_hom_r = (rare_copies - midpoint) / 2;
    current_hom_c = genotypes - current_hets - current_hom_r;
    while current_hets + 2 <= rare_copies {
        let next = probabilities[current_hets]
            * 4.0
            * current_hom_r as f64
            * current_hom_c as f64
            / ((current_hets + 2) as f64 * (current_hets + 1) as f64);
        probabilities[current_hets + 2] = next;
        probability_sum += next;
        current_hets += 2;
        current_hom_r = current_hom_r.saturating_sub(1);
        current_hom_c = current_hom_c.saturating_sub(1);
    }

    let observed_probability = probabilities.get(obs_hets).copied().unwrap_or(0.0);
    let epsilon = 1e-12;
    let tail_sum: f64 = probabilities.iter()
        .filter(|&&probability| probability <= observed_probability + epsilon)
        .sum();
    (tail_sum / probability_sum).clamp(0.0, 1.0)
}

/// Stream-copy retained variants from a BIM file into a new BIM writer.
pub fn filter_bim_stream<R, W>(
    reader: R,
    mut writer: W,
    variant_keep: &[bool],
) -> Result<(), String>
where
    R: BufRead,
    W: Write,
{
    let mut index = 0;
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if index < variant_keep.len() && variant_keep[index] {
            writeln!(writer, "{trimmed}").map_err(|e| e.to_string())?;
        }
        index += 1;
    }
    writer.flush().map_err(|e| e.to_string())?;
    Ok(())
}

/// Stream-copy retained variants and retained samples into a new BED writer.
pub fn write_filtered_bed<R, W>(
    mut bed_reader: R,
    mut bed_writer: W,
    n_samples: usize,
    sample_keep: &[bool],
    variant_keep: &[bool],
) -> Result<(), String>
where
    R: Read + Seek,
    W: Write,
{
    bed_writer.write_all(&BED_MAGIC).map_err(|e| e.to_string())?;
    bed_reader.seek(SeekFrom::Start(3)).map_err(|e| e.to_string())?;

    let bytes_per_snp = bytes_per_variant(n_samples);
    let mut in_buf = vec![0u8; bytes_per_snp];
    let mut out_buf = Vec::new();
    let all_samples_kept = sample_keep.iter().take(n_samples).all(|&k| k);

    for &is_kept in variant_keep {
        bed_reader.read_exact(&mut in_buf).map_err(|e| e.to_string())?;
        if is_kept {
            if all_samples_kept {
                bed_writer.write_all(&in_buf).map_err(|e| e.to_string())?;
            } else {
                pack_filtered_row(&in_buf, n_samples, sample_keep, &mut out_buf);
                bed_writer.write_all(&out_buf).map_err(|e| e.to_string())?;
            }
        }
    }
    bed_writer.flush().map_err(|e| e.to_string())?;
    Ok(())
}
