#!/usr/bin/env bash
set -euo pipefail

PLUGIN_ROOT="$(cd "$(dirname "$0")" && pwd)"
REPOSITORY_ROOT="$(cd "$PLUGIN_ROOT/../.." && pwd)"

if [[ -d "$REPOSITORY_ROOT/.toolchains/rustup" ]]; then
    export RUSTUP_HOME="${RUSTUP_HOME:-$REPOSITORY_ROOT/.toolchains/rustup}"
    export CARGO_HOME="${CARGO_HOME:-$REPOSITORY_ROOT/.toolchains/cargo}"
    export PATH="$CARGO_HOME/bin:$PATH"
fi

# Keep build output outside the plugin folder: Xcode copies the whole folder into the app.
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${TMPDIR:-/tmp}/xmosome-plugin-build/genokit}"
cd "$PLUGIN_ROOT"
mkdir -p "$PLUGIN_ROOT/Runtime"

echo "Compiling genokit-browser (wasm32-unknown-unknown)..."
cargo build --offline -p genokit-browser --release --target wasm32-unknown-unknown
cp "$CARGO_TARGET_DIR/wasm32-unknown-unknown/release/genokit_browser.wasm" "$PLUGIN_ROOT/Runtime/analysis.wasm"

echo "Compiling genokit-component (wasm32-wasip2)..."
cargo build --offline -p genokit-component --release --target wasm32-wasip2
cp "$CARGO_TARGET_DIR/wasm32-wasip2/release/genokit_component.wasm" "$PLUGIN_ROOT/Runtime/plugin.component.wasm"

echo "Successfully built genokit runtimes:"
echo "  $PLUGIN_ROOT/Runtime/analysis.wasm"
echo "  $PLUGIN_ROOT/Runtime/plugin.component.wasm"
