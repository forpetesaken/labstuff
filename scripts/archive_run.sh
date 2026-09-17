#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  printf 'Usage: %s CONFIG OUTPUT_DIR RUN_NAME\n' "$0" >&2
  exit 2
fi

config_path=$1
output_dir=$2
run_name=$3
repo_root=$(cd "$(dirname "$0")/.." && pwd)
timestamp=$(date '+%Y-%m-%d_%H%M')
archive_dir="$repo_root/outputs/archive/${timestamp}_${run_name}"

[[ -f "$config_path" ]] || { printf 'Config not found: %s\n' "$config_path" >&2; exit 1; }
[[ -d "$output_dir" ]] || { printf 'Output directory not found: %s\n' "$output_dir" >&2; exit 1; }

mkdir -p "$archive_dir/figures"
cp "$config_path" "$archive_dir/config.yaml"
cp -R "$output_dir"/. "$archive_dir/"

{
  printf 'Run date: %s\n' "$(date '+%Y-%m-%d %H:%M:%S %z')"
  printf 'Git commit: %s\n' "$(git -C "$repo_root" rev-parse HEAD 2>/dev/null || printf 'uncommitted')"
  printf 'Input config: %s\n' "$config_path"
  printf 'Source output: %s\n' "$output_dir"
  printf 'Command: %s %q %q %q\n' "$0" "$config_path" "$output_dir" "$run_name"
} > "$archive_dir/run_info.txt"

printf 'Archived run: %s\n' "$archive_dir"
