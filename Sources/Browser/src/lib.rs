//! genokit-browser: Browser WebAssembly adapter for Xmosome GUI.
//!
//! Communicates via JSON-in / JSON-out protocol in memory:
//!   - `pk_json_alloc(len)`: allocate input buffer
//!   - `pk_analyze_json(ptr, len)`: execute QC and return result length
//!   - `pk_json_result()`: return result pointer

#![allow(static_mut_refs)]

use std::io::Cursor;
use genokit_core::{
    generate_markdown_report, parse_bim, parse_fam, run_qc, write_filtered_bed, QcOptions,
};

static mut INPUT_BUFFER: Vec<u8> = Vec::new();
static mut OUTPUT_BUFFER: Vec<u8> = Vec::new();

#[unsafe(no_mangle)]
pub extern "C" fn pk_json_alloc(len: usize) -> *mut u8 {
    unsafe {
        INPUT_BUFFER.clear();
        INPUT_BUFFER.resize(len, 0);
        INPUT_BUFFER.as_mut_ptr()
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn pk_json_result() -> *const u8 {
    unsafe {
        OUTPUT_BUFFER.as_ptr()
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn pk_analyze_json(ptr: *const u8, len: usize) -> usize {
    let input_slice = unsafe { std::slice::from_raw_parts(ptr, len) };
    let result_str = match handle_request(input_slice) {
        Ok(s) => s,
        Err(e) => serde_json::json!({
            "status": "failure",
            "error": {
                "code": "analysis-error",
                "message": e
            }
        })
        .to_string(),
    };
    let bytes = result_str.into_bytes();
    let out_len = bytes.len();
    unsafe {
        OUTPUT_BUFFER = bytes;
    }
    out_len
}

fn handle_request(input_bytes: &[u8]) -> Result<String, String> {
    let request: serde_json::Value =
        serde_json::from_slice(input_bytes).map_err(|e| format!("Invalid JSON: {e}"))?;

    let args = request
        .get("arguments")
        .and_then(|a| a.as_object())
        .ok_or("Missing 'arguments' object")?;

    let fam_text = args
        .get("fam_text")
        .and_then(|v| v.as_str())
        .ok_or("Missing 'fam_text'")?;
    let bim_text = args
        .get("bim_text")
        .and_then(|v| v.as_str())
        .ok_or("Missing 'bim_text'")?;

    let bed_bytes = if let Some(base64_str) = args.get("bed_base64").and_then(|v| v.as_str()) {
        decode_base64(base64_str)?
    } else {
        return Err("Missing 'bed_base64'".into());
    };

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

    let options = QcOptions {
        mind,
        geno,
        hwe,
        maf,
        max_maf,
        order,
        rounds,
        make_bed,
    };

    let fam_records = parse_fam(fam_text);
    let bim_records = parse_bim(bim_text);

    let qc_output = run_qc(
        Cursor::new(&bed_bytes),
        &fam_records,
        &bim_records,
        &options,
        |_fraction, _stage, _message| {},
    )?;

    let report_content = generate_markdown_report(&qc_output.metrics, &options);

    let mut filtered_bed_base64 = None;
    if options.make_bed {
        let mut out_bed_bytes = Vec::new();
        write_filtered_bed(
            Cursor::new(&bed_bytes),
            &mut out_bed_bytes,
            fam_records.len(),
            &qc_output.sample_keep,
            &qc_output.variant_keep,
        )?;
        filtered_bed_base64 = Some(encode_base64(&out_bed_bytes));
    }

    // Agent / headless: omit filtered BED so sequential GenoKit → PedigreeKit runs
    // don't OOM WKWebView with a second multi‑MB base64 blob.
    let compact = args
        .get("compact_outputs")
        .and_then(|v| v.as_bool())
        .unwrap_or(false);
    let files = if compact {
        serde_json::json!({
            "filtered_fam": qc_output.filtered_fam,
            "filtered_bim": qc_output.filtered_bim,
            // Intentionally omit filtered_bed_base64 — host can re-run with make_bed if needed.
        })
    } else {
        serde_json::json!({
            "filtered_fam": qc_output.filtered_fam,
            "filtered_bim": qc_output.filtered_bim,
            "filtered_bed_base64": filtered_bed_base64
        })
    };

    let response = serde_json::json!({
        "status": "success",
        "message": format!(
            "QC completed: {}/{} samples, {}/{} variants retained",
            qc_output.metrics.retained_samples,
            qc_output.metrics.initial_samples,
            qc_output.metrics.retained_variants,
            qc_output.metrics.initial_variants
        ),
        "report": {
            "id": "report",
            "type": "markdown",
            "path": "report.md",
            "content": report_content
        },
        "metrics": qc_output.metrics,
        "files": files
    });

    Ok(response.to_string())
}

// Minimal standard Base64 codec
const B64_CHARS: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

fn encode_base64(data: &[u8]) -> String {
    let mut result = String::with_capacity((data.len() + 2) / 3 * 4);
    for chunk in data.chunks(3) {
        let b0 = chunk[0] as usize;
        let b1 = if chunk.len() > 1 { chunk[1] as usize } else { 0 };
        let b2 = if chunk.len() > 2 { chunk[2] as usize } else { 0 };

        let n = (b0 << 16) | (b1 << 8) | b2;
        result.push(B64_CHARS[(n >> 18) & 63] as char);
        result.push(B64_CHARS[(n >> 12) & 63] as char);
        if chunk.len() > 1 {
            result.push(B64_CHARS[(n >> 6) & 63] as char);
        } else {
            result.push('=');
        }
        if chunk.len() > 2 {
            result.push(B64_CHARS[n & 63] as char);
        } else {
            result.push('=');
        }
    }
    result
}

fn decode_base64(s: &str) -> Result<Vec<u8>, String> {
    let clean: Vec<u8> = s
        .bytes()
        .filter(|&b| !b.is_ascii_whitespace())
        .collect();

    if clean.len() % 4 != 0 {
        return Err("Invalid Base64 length".into());
    }

    let mut out = Vec::with_capacity(clean.len() / 4 * 3);
    for chunk in clean.chunks(4) {
        let mut val = 0u32;
        let mut pad = 0;
        for i in 0..4 {
            let c = chunk[i];
            let d = match c {
                b'A'..=b'Z' => c - b'A',
                b'a'..=b'z' => c - b'a' + 26,
                b'0'..=b'9' => c - b'0' + 52,
                b'+' => 62,
                b'/' => 63,
                b'=' => {
                    pad += 1;
                    0
                }
                _ => return Err(format!("Invalid Base64 character: {c}")),
            };
            val = (val << 6) | (d as u32);
        }
        out.push(((val >> 16) & 0xFF) as u8);
        if pad < 2 {
            out.push(((val >> 8) & 0xFF) as u8);
        }
        if pad < 1 {
            out.push((val & 0xFF) as u8);
        }
    }
    Ok(out)
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
