#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$script_dir/.." && pwd)"
read -r -a python_includes <<<"$(python3-config --includes)"
extension_suffix="$(python3-config --extension-suffix)"
gcc -O3 -shared -fPIC "${python_includes[@]}" \
  "$root/src/alscn_switch/_xyzstream.c" \
  -o "$root/src/alscn_switch/_xyzstream${extension_suffix}"
printf 'Built %s\n' "$root/src/alscn_switch/_xyzstream${extension_suffix}"
