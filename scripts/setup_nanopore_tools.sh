#!/usr/bin/env bash
# Install the main-evaluation toolchain into ./tools (not tracked; large):
#   squigulator v0.5.0 (prebuilt x86_64), Dorado (ONT CDN tarball), Dorado R10.4.1 e8.2
#   400bps models, and blue-crab (BLOW5 -> POD5) in an isolated venv (tools/venv).
# Records versions and sha256 sums in tools/VERSIONS.json. Re-runnable; skips finished steps.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS="$ROOT/tools"
DL="$TOOLS/dl"
SQUIG_VER="v0.5.0"
DORADO_VER="${DORADO_VER:-2.1.2}"
MODELS="${MODELS:-hac fast}"   # add 'sup' if wanted (large; see docs/notes/nanopore_pipeline.md)
mkdir -p "$DL" "$TOOLS/models"

fetch() {  # url dest
  [ -s "$2" ] || curl -fL --retry 5 --retry-delay 5 -o "$2" "$1"
}

# squigulator
fetch "https://github.com/hasindu2008/squigulator/releases/download/${SQUIG_VER}/squigulator-${SQUIG_VER}-x86_64-linux-binaries.tar.gz" \
      "$DL/squigulator-${SQUIG_VER}-x86_64-linux-binaries.tar.gz"
[ -x "$TOOLS/squigulator-${SQUIG_VER}/squigulator" ] || tar xzf "$DL/squigulator-${SQUIG_VER}-x86_64-linux-binaries.tar.gz" -C "$TOOLS"

# Dorado (note: the CDN host is cdn.oxfordnanoportal.com)
fetch "https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VER}-linux-x64.tar.gz" \
      "$DL/dorado-${DORADO_VER}-linux-x64.tar.gz"
[ -x "$TOOLS/dorado-${DORADO_VER}-linux-x64/bin/dorado" ] || tar xzf "$DL/dorado-${DORADO_VER}-linux-x64.tar.gz" -C "$TOOLS"
DORADO="$TOOLS/dorado-${DORADO_VER}-linux-x64/bin/dorado"

# models: newest R10.4.1 e8.2 400bps model of each requested kind
for kind in $MODELS; do
  name=$("$DORADO" download --list-yaml 2>/dev/null | grep -oE "dna_r10.4.1_e8.2_400bps_${kind}@v[0-9.]+\"" | tr -d '"' | sort -uV | tail -1)
  [ -n "$name" ] || { echo "no ${kind} model listed by dorado" >&2; exit 1; }
  [ -d "$TOOLS/models/$name" ] || "$DORADO" download --model "$name" --models-directory "$TOOLS/models"
done

# blue-crab (BLOW5 -> POD5), isolated from the project venv
[ -x "$TOOLS/venv/bin/blue-crab" ] || { python3 -m venv "$TOOLS/venv"; "$TOOLS/venv/bin/pip" install -q blue-crab==0.6.0; }

python3 - "$TOOLS" "$DORADO" <<'EOF'
import hashlib, json, subprocess, sys, glob, os
tools, dorado = sys.argv[1], sys.argv[2]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
v = {"archives": {os.path.basename(p): sha(p) for p in glob.glob(os.path.join(tools, "dl", "*.tar.gz"))}}
v["squigulator"] = subprocess.run([glob.glob(os.path.join(tools, "squigulator-*", "squigulator"))[0], "--version"],
                                  capture_output=True, text=True).stdout.strip()
r = subprocess.run([dorado, "--version"], capture_output=True, text=True)
v["dorado"] = next((l.strip() for l in (r.stdout + r.stderr).splitlines() if l.strip() and not l.startswith("[")), "?")
v["models"] = sorted(os.path.basename(p) for p in glob.glob(os.path.join(tools, "models", "dna_*")) if os.path.isdir(p))
v["blue_crab"] = subprocess.run([os.path.join(tools, "venv", "bin", "pip"), "show", "blue-crab"], capture_output=True,
                                text=True).stdout.split("Version: ")[1].split()[0]
v["nvidia_smi"] = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"],
                                 capture_output=True, text=True).stdout.strip()
json.dump(v, open(os.path.join(tools, "VERSIONS.json"), "w"), indent=2)
print(json.dumps(v, indent=2))
EOF
