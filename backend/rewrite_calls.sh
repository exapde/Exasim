#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./apply_model_to_drivers.sh               # operate on current directory
#   ./apply_model_to_drivers.sh /path/to/src # operate on given root
ROOT_DIR="${1:-.}"

# List of driver function names to make templated at call sites.
DRIVERS=(
  FluxDriver
  SourceDriver
  SourcewDriver
  OutputDriver
  MonitorDriver
  AvfieldDriver
  EosDriver
  EosduDriver
  EosdwDriver
  TdfuncDriver
  FhatDriver
  FbouDriver
  UhatDriver
  UbouDriver
  InitodgDriver
  InitqDriver
  InitudgDriver
  InituDriver
  InitwdgDriver
  cpuInitodgDriver
  cpuInitqDriver
  cpuInitudgDriver
  cpuInituDriver
  cpuInitwdgDriver
  VisScalarsDriver
  VisVectorsDriver
  VisTensorsDriver
  QoIvolumeDriver
  QoIboundaryDriver
  FintDriver
)

# Detect sed -i flavor (GNU vs BSD/macOS)
if sed --version >/dev/null 2>&1; then
  # GNU sed
  SED_INPLACE=(sed -i)
else
  # BSD / macOS sed
  SED_INPLACE=(sed -i '')
fi

# Helper to run sed -i only if we have files
run_sed_inplace() {
  local pattern="$1"
  shift
  local files=("$@")
  if ((${#files[@]} > 0)); then
    "${SED_INPLACE[@]}" "$pattern" "${files[@]}"
  fi
}

for name in "${DRIVERS[@]}"; do
  echo ">>> Processing ${name}"

  # 1) Protect declarations/definitions: change
  #    void FooDriver(
  # to
  #    void FooDriver_DECL(
  # We only touch lines that start with optional spaces + 'void FooDriver('.
  mapfile -t decl_files < <(
    grep -rl "void ${name}(" "$ROOT_DIR" --include='*.cpp' --include='*.cc' --include='*.cxx' --include='*.hpp' --include='*.h' || true
  )
  run_sed_inplace "s/^[[:space:]]*void[[:space:]]\\+${name}(/void ${name}_DECL(/" "${decl_files[@]}"

  # 2) Replace all *call sites*:
  #    FooDriver(
  # -> FooDriver<Model>(
  #
  # By this point, decl/def lines have been renamed to FooDriver_DECL, so they
  # won't be affected by this substitution.
  mapfile -t call_files < <(
    grep -rl "${name}(" "$ROOT_DIR" --include='*.cpp' --include='*.cc' --include='*.cxx' --include='*.hpp' --include='*.h' || true
  )
  run_sed_inplace "s/${name}(/${name}<Model>(/g" "${call_files[@]}"

  # 3) Restore declaration/definition name:
  #    void FooDriver_DECL(
  # -> void FooDriver(
  mapfile -t decl_rest_files < <(
    grep -rl "void ${name}_DECL(" "$ROOT_DIR" --include='*.cpp' --include='*.cc' --include='*.cxx' --include='*.hpp' --include='*.h' || true
  )
  run_sed_inplace "s/void ${name}_DECL(/void ${name}(/" "${decl_rest_files[@]}"

done

echo "Done."
