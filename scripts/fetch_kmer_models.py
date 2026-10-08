"""Download ONT k-mer level tables at a pinned commit and record provenance.

Usage: python scripts/fetch_kmer_models.py

Never substitutes or synthesises a table. If a download fails, the script exits non-zero.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

COMMIT = "4e56daed7fbb79b538f58e41262d5c54b07356ea"  # nanoporetech/kmer_models master, 2023-06-06
BASE = f"https://raw.githubusercontent.com/nanoporetech/kmer_models/{COMMIT}"
FILES = {
    "r10.4.1_9mer": "dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt",
    "r9.4.1_6mer": "legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model",
}
OUT = Path(__file__).resolve().parent.parent / "data" / "kmer_models"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    prov = {"repository": "https://github.com/nanoporetech/kmer_models", "commit": COMMIT, "files": {}}
    for name, rel in FILES.items():
        url = f"{BASE}/{rel}"
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                blob = r.read()
        except Exception as exc:  # noqa: BLE001
            print(f"FAILED to download {url}: {exc}. Not inventing a table; stopping.", file=sys.stderr)
            return 1
        dest = OUT / Path(rel).name
        dest.write_bytes(blob)
        prov["files"][name] = {
            "path_in_repo": rel,
            "url": url,
            "local_file": dest.name,
            "sha256": hashlib.sha256(blob).hexdigest(),
            "size_bytes": len(blob),
            "downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        }
        print(f"{name}: {len(blob)} bytes, sha256 {prov['files'][name]['sha256']}")
    (OUT / "PROVENANCE.json").write_text(json.dumps(prov, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
