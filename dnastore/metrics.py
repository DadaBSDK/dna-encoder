"""Density and mass metrics (DESIGN.md §2.5).

Definitions
-----------
* ``code_density`` (net, bits/nt): stored (post-compression) bits / payload-region nt of
  the data part of the bit-packed stream (data rows' share of oligos x payload nt).
  Excludes primers, index+seed, header oligos and RS parity. The denominator includes the
  in-frame CRC, row padding and steering bases, so codecs are charged for them.
* ``raw_rate``: frame bits / payload nt (the codec's own rate, CRC included).
* ``effective_density`` (bits/nt): original file bits / all synthesised nt.
* ``stored_effective_density``: stored bits / all synthesised nt. Compares codecs without
  the effect of compression.
* Mass: ``n_oligos * copies * length * MW / N_A`` with MW = 330 g/mol per nt (ssDNA) or
  650 g/mol per bp (dsDNA).
"""

from __future__ import annotations

from typing import Any

from .ecc import CRC_BITS

AVOGADRO = 6.02214076e23
MW_SS_NT = 330.0
MW_DS_BP = 650.0


def pool_metrics(res: Any, copies_per_oligo: int = 10) -> dict[str, Any]:
    """Summary metrics for an :class:`~dnastore.encoder.EncodeResult`."""
    geom, layout, info = res.geometry, res.layout, res.info
    n_oligos = len(res.oligos)
    total_nt = sum(len(s) for _, s in res.oligos)
    data_payload_nt = layout.data_oligo_equiv * geom.payload_nt
    stored_bits, orig_bits = 8 * info.stored_len, 8 * info.original_len
    return {
        "codec": res.codec.name,
        "codec_params": res.codec.params(),
        "file_type": info.file_type,
        "original_bytes": info.original_len,
        "stored_bytes": info.stored_len,
        "compressed": info.compressed,
        "zstd_trial_ratio": info.trial_ratio,
        "oligo_len": geom.oligo_len,
        "n_oligos": n_oligos,
        "n_header_oligos": res.n_header_oligos,
        "n_data_oligos": layout.n_oligos,
        "n_data_rows": layout.n_data_rows,
        "n_parity_rows": layout.n_parity_rows,
        "rs_blocks": len(layout.blocks),
        "payload_bits_per_oligo": layout.payload_bits,
        "raw_rate": (layout.payload_bits + CRC_BITS) / geom.payload_nt,
        "seed_exhaustion_rate": (sum(bool(i.get("seed_exhausted")) for i in res.payload_info if not i.get("header"))
                                 / max(1, sum(1 for i in res.payload_info if not i.get("header")))),
        "mean_tries": (sum(i.get("tries", 1) for i in res.payload_info if not i.get("header"))
                       / max(1, sum(1 for i in res.payload_info if not i.get("header")))),
        "payload_nt_per_oligo": geom.payload_nt,
        "total_nt": total_nt,
        "code_density": stored_bits / data_payload_nt if data_payload_nt else 0.0,
        "effective_density": orig_bits / total_nt if total_nt else 0.0,
        "stored_effective_density": stored_bits / total_nt if total_nt else 0.0,
        "mass_ssdna_g": total_nt * copies_per_oligo * MW_SS_NT / AVOGADRO,
        "mass_dsdna_g": total_nt * copies_per_oligo * MW_DS_BP / AVOGADRO,
        "copies_per_oligo": copies_per_oligo,
        **{f"t_{k}": v for k, v in res.timings.items()},
    }
