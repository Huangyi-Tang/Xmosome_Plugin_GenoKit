use std::{
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Cursor, Write},
    path::{Path, PathBuf},
};

wit_bindgen::generate!({ path: "../../WIT", world: "analysis-plugin" });

use exports::xmosome::plugin::analysis::{Guest, PluginError};
use xmosome::plugin::host::{emit_progress, ProgressEvent};
use genokit_core::BimRecord;

struct GenoKit;

impl Guest for GenoKit {
    fn run(request_json: String) -> Result<String, PluginError> {
        run_analysis(&request_json).map_err(|message| PluginError {
            code: "qc-failed".into(),
            message,
        })
    }
}

enum GenotypeSource {
    Binary {
        bed_path: PathBuf,
        bim_path: PathBuf,
        fam_path: PathBuf,
    },
    Text {
        ped_path: PathBuf,
        map_path: PathBuf,
    },
}

enum BimSource {
    Path(PathBuf),
    Records(Vec<BimRecord>),
}

enum BedSource {
    Path(PathBuf, usize),
    Buffer(Vec<u8>, usize),
}

fn run_analysis(request_json: &str) -> Result<String, String> {
    let request: serde_json::Value =
        serde_json::from_str(request_json).map_err(|e| format!("Invalid JSON request: {e}"))?;

    let args = request
        .get("arguments")
        .and_then(|a| a.as_object())
        .ok_or("Missing 'arguments' in request")?;

    // Resolve file paths
    let inputs = resolve_inputs(args)?;

    // Explicit JSON null disables a filter; a missing key keeps the historical default.
    let mind = optional_threshold(args, "mind", 0.1);
    let geno = optional_threshold(args, "geno", 0.1);
    let hwe = args.get("hwe").and_then(|v| v.as_f64());
    let maf = optional_threshold(args, "maf", 0.01);
    let max_maf = args.get("max_maf").and_then(|v| v.as_f64());
    let rounds = args.get("rounds").and_then(|v| v.as_u64()).unwrap_or(1) as usize;
    let make_bed = args.get("make_bed").and_then(|v| v.as_bool()).unwrap_or(true);

    let order = if let Some(order_str) = args.get("order").and_then(|v| v.as_str()) {
        order_str
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        vec!["mind".into(), "geno".into(), "hwe".into(), "maf".into()]
    };

    let options = genokit_core::QcOptions {
        mind,
        geno,
        hwe,
        maf,
        max_maf,
        order,
        rounds,
        make_bed,
    };

    let output_dir_str = request
        .get("output")
        .and_then(|o| o.get("directory"))
        .and_then(|d| d.as_str())
        .ok_or("Missing 'output.directory' in request")?;
    let output_dir = Path::new(output_dir_str);
    fs::create_dir_all(output_dir).map_err(|e| format!("Failed to create output dir: {e}"))?;
    let data_dir = output_dir.join("data");
    fs::create_dir_all(&data_dir).map_err(|e| format!("Failed to create data dir: {e}"))?;

    let (metrics, sample_keep, variant_keep, filtered_fam, bim_source, bed_source) = match inputs {
        GenotypeSource::Binary { bed_path, bim_path, fam_path } => {
            emit_progress(&ProgressEvent {
                fraction: Some(0.02),
                stage: Some("init".into()),
                message: "Reading FAM and BIM annotation files".into(),
            });

            let fam_content = fs::read_to_string(&fam_path)
                .map_err(|e| format!("Failed to read FAM file {}: {e}", fam_path.display()))?;
            let fam_records = genokit_core::parse_fam(&fam_content);

            let bim_file = File::open(&bim_path)
                .map_err(|e| format!("Failed to open BIM file {}: {e}", bim_path.display()))?;
            let bim_reader = BufReader::new(bim_file);
            let mut n_variants = 0usize;
            for line in bim_reader.lines() {
                let l = line.map_err(|e| e.to_string())?;
                let t = l.trim();
                if !t.is_empty() && !t.starts_with('#') {
                    n_variants += 1;
                }
            }

            let bed_file = File::open(&bed_path)
                .map_err(|e| format!("Failed to open BED file {}: {e}", bed_path.display()))?;
            let mut bed_reader = BufReader::with_capacity(1024 * 1024, bed_file);

            let (metrics, sample_keep, variant_keep, filtered_fam) = genokit_core::run_qc_streaming_core(
                &mut bed_reader,
                &fam_records,
                n_variants,
                &options,
                |fraction, stage, message| {
                    emit_progress(&ProgressEvent {
                        fraction: Some(fraction),
                        stage: Some(stage.into()),
                        message: message.into(),
                    });
                },
            )?;
            (
                metrics,
                sample_keep,
                variant_keep,
                filtered_fam,
                BimSource::Path(bim_path),
                BedSource::Path(bed_path, fam_records.len()),
            )
        }
        GenotypeSource::Text { ped_path, map_path } => {
            emit_progress(&ProgressEvent {
                fraction: Some(0.02),
                stage: Some("init".into()),
                message: "Parsing and converting PED and MAP text files to binary BED".into(),
            });

            let map_file = File::open(&map_path)
                .map_err(|e| format!("Failed to open MAP file {}: {e}", map_path.display()))?;
            let ped_file = File::open(&ped_path)
                .map_err(|e| format!("Failed to open PED file {}: {e}", ped_path.display()))?;

            let converted = genokit_core::convert_ped_map(
                BufReader::new(map_file),
                BufReader::new(ped_file),
            )?;

            let n_samples = converted.fam_records.len();
            let n_variants = converted.bim_records.len();
            let mut bed_cursor = Cursor::new(converted.bed_data);

            let (metrics, sample_keep, variant_keep, filtered_fam) = genokit_core::run_qc_streaming_core(
                &mut bed_cursor,
                &converted.fam_records,
                n_variants,
                &options,
                |fraction, stage, message| {
                    emit_progress(&ProgressEvent {
                        fraction: Some(fraction),
                        stage: Some(stage.into()),
                        message: message.into(),
                    });
                },
            )?;
            (
                metrics,
                sample_keep,
                variant_keep,
                filtered_fam,
                BimSource::Records(converted.bim_records),
                BedSource::Buffer(bed_cursor.into_inner(), n_samples),
            )
        }
    };

    emit_progress(&ProgressEvent {
        fraction: Some(0.85),
        stage: Some("write-outputs".into()),
        message: "Writing filtered datasets and QC report".into(),
    });

    // Write report.md
    let report_content = genokit_core::generate_markdown_report(&metrics, &options);
    fs::write(output_dir.join("report.md"), report_content)
        .map_err(|e| format!("Failed to write report.md: {e}"))?;

    // Write data/qc-summary.json
    let summary_json = serde_json::to_string_pretty(&metrics)
        .map_err(|e| format!("Failed to serialize metrics: {e}"))?;
    fs::write(data_dir.join("qc-summary.json"), summary_json)
        .map_err(|e| format!("Failed to write qc-summary.json: {e}"))?;

    // Write data/filtered.fam
    fs::write(data_dir.join("filtered.fam"), &filtered_fam)
        .map_err(|e| format!("Failed to write filtered.fam: {e}"))?;

    // Write data/filtered.bim
    let bim_out = File::create(data_dir.join("filtered.bim")).map_err(|e| e.to_string())?;
    let mut bim_writer = BufWriter::new(bim_out);
    match bim_source {
        BimSource::Path(bim_path) => {
            let bim_in = File::open(&bim_path).map_err(|e| e.to_string())?;
            genokit_core::filter_bim_stream(BufReader::new(bim_in), bim_writer, &variant_keep)?;
        }
        BimSource::Records(records) => {
            for (j, is_kept) in variant_keep.iter().enumerate() {
                if *is_kept && j < records.len() {
                    writeln!(bim_writer, "{}", records[j].to_bim_line()).map_err(|e| e.to_string())?;
                }
            }
            bim_writer.flush().map_err(|e| e.to_string())?;
        }
    }

    let mut outputs = vec![
        serde_json::json!({
            "id": "summary",
            "type": "json",
            "path": "data/qc-summary.json"
        }),
        serde_json::json!({
            "id": "filtered_fam",
            "type": "text",
            "path": "data/filtered.fam"
        }),
        serde_json::json!({
            "id": "filtered_bim",
            "type": "text",
            "path": "data/filtered.bim"
        }),
    ];

    // Optionally write filtered BED
    if options.make_bed {
        let out_bed = File::create(data_dir.join("filtered.bed"))
            .map_err(|e| format!("Failed to create filtered.bed: {e}"))?;
        let out_bed_writer = BufWriter::with_capacity(1024 * 1024, out_bed);

        match bed_source {
            BedSource::Path(bed_path, n_samples) => {
                let in_bed = File::open(&bed_path)
                    .map_err(|e| format!("Failed to re-open BED file: {e}"))?;
                let in_bed_reader = BufReader::with_capacity(1024 * 1024, in_bed);
                genokit_core::write_filtered_bed(
                    in_bed_reader,
                    out_bed_writer,
                    n_samples,
                    &sample_keep,
                    &variant_keep,
                )?;
            }
            BedSource::Buffer(bed_data, n_samples) => {
                let in_bed_reader = Cursor::new(bed_data);
                genokit_core::write_filtered_bed(
                    in_bed_reader,
                    out_bed_writer,
                    n_samples,
                    &sample_keep,
                    &variant_keep,
                )?;
            }
        }

        outputs.push(serde_json::json!({
            "id": "filtered_bed",
            "type": "binary",
            "path": "data/filtered.bed"
        }));
    }

    // Machine-readable card manifest (Agent / GUI shared layout).
    {
        let mut arts = vec![
            serde_json::json!({"id":"report","title":"QC report","kind":"report","path":"output/report.md","role":"report"}),
            serde_json::json!({"id":"qc-summary","title":"QC metrics","kind":"qc","path":"output/data/qc-summary.json","role":"qc"}),
            serde_json::json!({"id":"filtered-fam","title":"Filtered FAM","kind":"table","path":"output/data/filtered.fam","role":"table"}),
            serde_json::json!({"id":"filtered-bim","title":"Filtered BIM","kind":"table","path":"output/data/filtered.bim","role":"table"}),
        ];
        if options.make_bed {
            arts.push(serde_json::json!({"id":"filtered-bed","title":"Filtered BED","kind":"table","path":"output/data/filtered.bed","role":"table"}));
        }
        let manifest = serde_json::json!({ "version": 1, "artifacts": arts });
        let parent = output_dir
            .parent()
            .unwrap_or(output_dir)
            .to_path_buf();
        let _ = fs::write(
            parent.join("artifacts.json"),
            serde_json::to_string_pretty(&manifest).unwrap_or_else(|_| "{}".into()),
        );
    }

    emit_progress(&ProgressEvent {
        fraction: Some(1.0),
        stage: Some("completed".into()),
        message: "Quality control completed successfully".into(),
    });

    let response = serde_json::json!({
        "status": "success",
        "message": format!(
            "QC completed: {}/{} samples, {}/{} variants retained",
            metrics.retained_samples,
            metrics.initial_samples,
            metrics.retained_variants,
            metrics.initial_variants
        ),
        "report": {
            "id": "report",
            "type": "markdown",
            "path": "report.md"
        },
        "outputs": outputs,
        "metrics": {
            "initial_samples": metrics.initial_samples,
            "retained_samples": metrics.retained_samples,
            "removed_by_mind": metrics.removed_by_mind,
            "initial_variants": metrics.initial_variants,
            "retained_variants": metrics.retained_variants,
            "removed_by_geno": metrics.removed_by_geno,
            "removed_by_hwe": metrics.removed_by_hwe,
            "removed_by_maf": metrics.removed_by_maf,
            "rounds_executed": metrics.rounds_executed
        }
    });

    Ok(response.to_string())
}

fn resolve_inputs(
    args: &serde_json::Map<String, serde_json::Value>,
) -> Result<GenotypeSource, String> {
    if let (Some(bed), Some(bim), Some(fam)) = (
        args.get("bed").and_then(|v| v.as_str()),
        args.get("bim").and_then(|v| v.as_str()),
        args.get("fam").and_then(|v| v.as_str()),
    ) {
        return Ok(GenotypeSource::Binary {
            bed_path: PathBuf::from(bed),
            bim_path: PathBuf::from(bim),
            fam_path: PathBuf::from(fam),
        });
    }

    if let (Some(ped), Some(map)) = (
        args.get("ped").and_then(|v| v.as_str()),
        args.get("map").and_then(|v| v.as_str()),
    ) {
        return Ok(GenotypeSource::Text {
            ped_path: PathBuf::from(ped),
            map_path: PathBuf::from(map),
        });
    }

    if let Some(bfile) = args.get("bfile").and_then(|v| v.as_str()) {
        let bed = PathBuf::from(format!("{bfile}.bed"));
        let bim = PathBuf::from(format!("{bfile}.bim"));
        let fam = PathBuf::from(format!("{bfile}.fam"));
        if bed.exists() && bim.exists() && fam.exists() {
            return Ok(GenotypeSource::Binary { bed_path: bed, bim_path: bim, fam_path: fam });
        }
        let ped = PathBuf::from(format!("{bfile}.ped"));
        let map = PathBuf::from(format!("{bfile}.map"));
        if ped.exists() && map.exists() {
            return Ok(GenotypeSource::Text { ped_path: ped, map_path: map });
        }
        return Err(format!(
            "Prefix '{bfile}' does not resolve to existing (.bed/.bim/.fam) or (.ped/.map) files"
        ));
    }

    if let Some(file) = args.get("file").and_then(|v| v.as_str()) {
        let ped = PathBuf::from(format!("{file}.ped"));
        let map = PathBuf::from(format!("{file}.map"));
        if ped.exists() && map.exists() {
            return Ok(GenotypeSource::Text { ped_path: ped, map_path: map });
        }
        return Err(format!(
            "Prefix '{file}' does not resolve to existing .ped and .map files"
        ));
    }

    Err("Must provide either binary filesets (--bfile, --bed/--bim/--fam) or text filesets (--file, --ped/--map)".into())
}

/// Missing key → historical default; JSON `null` → filter disabled.
fn optional_threshold(
    args: &serde_json::Map<String, serde_json::Value>,
    key: &str,
    default: f64,
) -> Option<f64> {
    match args.get(key) {
        None => Some(default),
        Some(value) if value.is_null() => None,
        Some(value) => value.as_f64(),
    }
}

export!(GenoKit);
