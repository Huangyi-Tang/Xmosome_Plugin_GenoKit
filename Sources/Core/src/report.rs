//! QC Report generation in Markdown and structured JSON.

use crate::qc::{QcMetrics, QcOptions};

pub fn generate_markdown_report(metrics: &QcMetrics, options: &QcOptions) -> String {
    let sample_retention_pct = if metrics.initial_samples > 0 {
        (metrics.retained_samples as f64 / metrics.initial_samples as f64) * 100.0
    } else {
        0.0
    };

    let variant_retention_pct = if metrics.initial_variants > 0 {
        (metrics.retained_variants as f64 / metrics.initial_variants as f64) * 100.0
    } else {
        0.0
    };

    let mean_sample_miss = metrics.mean_sample_missingness;
    let mean_variant_miss = metrics.mean_variant_missingness;
    let mean_maf = metrics.mean_maf;

    let max_maf_str = options
        .max_maf
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "Disabled".to_string());
    let hwe_str = options
        .hwe
        .map(|v| format!("{v:.3e}"))
        .unwrap_or_else(|| "Disabled".to_string());
    let mind_str = options
        .mind
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "Disabled".to_string());
    let geno_str = options
        .geno
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "Disabled".to_string());
    let maf_str = options
        .maf
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "Disabled".to_string());

    format!(
r#"# GenoKit: Genotype Quality Control Report

This report summarizes the genotype quality control (QC) filtering performed by **GenoKit** (`genokit`), aligned with the statistical definitions and 2-bit binary specifications of PLINK 1.9.

---

## 1. Executive Summary

| Category | Initial Count | Retained Count | Removed Count | Retention Rate |
| :--- | :--- | :--- | :--- | :--- |
| **Individuals / Samples** | {init_samples} | **{ret_samples}** | {rem_mind} (via `--mind`) | **{sample_pct:.2}%** |
| **Variants / SNPs** | {init_variants} | **{ret_variants}** | {rem_variants_total} (geno/HWE/maf) | **{variant_pct:.2}%** |

### Breakdown of Variant Exclusions
- **High missingness (`--geno > {geno_str}`)**: {rem_geno} variants removed
- **Hardy-Weinberg exact test (`--hwe < {hwe_str}`)**: {rem_hwe} variants removed
- **Low minor allele frequency (`--maf < {maf_str}`)**: {rem_maf} variants removed

---

## 2. QC Parameters & Execution Details

| Parameter | Configured Value | Description |
| :--- | :--- | :--- |
| `--mind` | `{mind_str}` | Maximum allowable individual missing genotype rate; disabled when omitted |
| `--geno` | `{geno_str}` | Maximum allowable variant missing call rate; disabled when omitted |
| `--hwe` | `{hwe_str}` | Minimum founder/control HWE exact-test p-value; disabled when omitted |
| `--maf` | `{maf_str}` | Minimum minor allele frequency threshold; disabled when omitted |
| `--max-maf` | `{max_maf_str}` | Maximum minor allele frequency threshold |
| `--order` | `{filter_order}` | Sequential QC step execution pipeline |
| `--rounds` | `{rounds_exec} / {max_rounds}` | Iterations until convergence |
| `--make-bed` | `{make_bed}` | Generate filtered `.bed`, `.bim`, and `.fam` fileset |

---

## 3. Post-QC Retained Distribution Metrics

- **Mean Individual Missingness**: `{mean_sample_miss:.4}` ({pct_sample_miss:.2}%)
- **Mean Variant Missingness**: `{mean_variant_miss:.4}` ({pct_variant_miss:.2}%)
- **Mean Minor Allele Frequency (MAF)**: `{mean_maf:.4}`

---

## 4. Generated Artifacts

![QC summary](data/qc-summary.json)

![Filtered FAM](data/filtered.fam)

![Filtered BIM](data/filtered.bim)

![Filtered BED](data/filtered.bed)

> **Methodology Note**: Genotype filtering was computed using pure streaming passes over the raw PLINK BED binary with $O(N)$ memory overhead without materializing the matrix in RAM.
"#,
        init_samples = metrics.initial_samples,
        ret_samples = metrics.retained_samples,
        rem_mind = metrics.removed_by_mind,
        sample_pct = sample_retention_pct,
        init_variants = metrics.initial_variants,
        ret_variants = metrics.retained_variants,
        rem_variants_total = metrics.removed_by_geno + metrics.removed_by_hwe + metrics.removed_by_maf,
        variant_pct = variant_retention_pct,
        mind_str = mind_str,
        geno_str = geno_str,
        maf_str = maf_str,
        rem_geno = metrics.removed_by_geno,
        rem_hwe = metrics.removed_by_hwe,
        rem_maf = metrics.removed_by_maf,
        hwe_str = hwe_str,
        max_maf_str = max_maf_str,
        filter_order = options.order.join(", "),
        rounds_exec = metrics.rounds_executed,
        max_rounds = options.rounds,
        make_bed = options.make_bed,
        mean_sample_miss = mean_sample_miss,
        pct_sample_miss = mean_sample_miss * 100.0,
        mean_variant_miss = mean_variant_miss,
        pct_variant_miss = mean_variant_miss * 100.0,
        mean_maf = mean_maf,
    )
}
