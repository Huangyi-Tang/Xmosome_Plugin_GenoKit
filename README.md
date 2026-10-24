# GenoKit.plugin — PLINK-Aligned Streaming Genotype QC

An independent, universal Xmosome plugin providing PLINK 1.9-aligned genotype quality control (`mind`, `geno`, `hwe`, `maf`, and `--make-bed`), implemented in 100% portable pure Rust with $O(N)$ streaming memory overhead.

## Architecture

Conforms strictly to the Xmosome Dual-Target Architecture:

```text
                         genokit-core
                 (Streaming BED I/O + QC Logic)
                               |
               +---------------+---------------+
               |                               |
        WASI Component                  Browser WASM
     plugin.component.wasm             analysis.wasm
               |                               |
      macOS / Linux CLI            macOS / iOS / iPadOS GUI
```

## Directory Structure

```text
3rdPartyPlugins/GenoKit.plugin/
├── Cargo.toml                 # Cargo workspace definition
├── manifest.json              # Plugin manifest (manifest_version 2)
├── build.sh                   # Builds wasm32-wasip2 and wasm32-unknown-unknown
├── WIT/
│   └── xmosome-plugin.wit     # WIT contract: xmosome:plugin/analysis-plugin@1.0.0
├── Sources/
│   ├── Core/                  # Pure Rust scientific core (lib.rs, bed.rs, bim_fam.rs, qc.rs, report.rs)
│   ├── Component/             # WASI 0.2 component adapter
│   ├── Browser/               # Browser WASM export adapter
│   └── GUI/                   # Web interface (index.html, app.js, app.css, execution.js)
├── Runtime/
│   ├── analysis.wasm          # Compiled browser WASM target
│   └── plugin.component.wasm  # Compiled WASI 0.2 component
└── Tests/
    └── Fixtures/              # Test datasets (toy.bed, toy.bim, toy.fam)
```

## Supported CLI Flags & Options

Command: `xmosome genokit [OPTIONS]`

| Flag | Type | Default | Description |
|---|---|---|---|
| `--bed` | file | required (or bfile) | PLINK 1.9 binary genotype file (`.bed`) |
| `--bim` | file | required (or bfile) | PLINK 1.9 variant coordinate file (`.bim`) |
| `--fam` | file | required (or bfile) | PLINK 1.9 pedigree & phenotype file (`.fam`) |
| `--bfile` | string | optional | Prefix for `.bed`, `.bim`, `.fam` fileset |
| `--mind` | float | `0.1` | Exclude individuals with missing genotype rate > threshold |
| `--geno` | float | `0.1` | Exclude variants with missing call rate > threshold |
| `--hwe` | float | optional | Exclude variants with founder/control HWE exact-test p-value below threshold |
| `--maf` | float | `0.01` | Exclude variants with minor allele frequency < threshold |
| `--max-maf` | float | optional | Exclude variants with minor allele frequency > threshold |
| `--order` | string | `mind,geno,hwe,maf`| Comma-separated QC filter sequence |
| `--rounds`| int | `1` | Maximum iteration rounds until convergence |
| `--make-bed` | bool | `true` | Generate filtered `.bed`, `.bim`, and `.fam` files |
| `--out` | dir | `./` | Output directory for report and filtered files |
| `--json` | flag | false | Output structured execution telemetry JSON |
| `--quiet` | flag | false | Suppress banner and progress logs |

## Building

```bash
cd 3rdPartyPlugins/GenoKit.plugin
./build.sh
```

## Running Tests

```bash
cargo test --manifest-path 3rdPartyPlugins/GenoKit.plugin/Sources/Core/Cargo.toml
```

## License

MIT License. Clean-room implementation referencing publicly documented statistical algorithms and file formats of PLINK 1.9.
