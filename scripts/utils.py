"""Shared helpers: fasta/stockholm I/O, domtblout parsing, publication plot style."""
from __future__ import annotations

import gzip
import os
import sys as _sys
from typing import Dict, Iterator, List, Tuple

# ---------------------------------------------------------------------------
# Sequence I/O
# ---------------------------------------------------------------------------

def _opener(path: str):
    return gzip.open if str(path).endswith((".gz", ".gzip")) else open


def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence). Header is everything after '>' minus newline."""
    # errors="replace": a single non-ASCII byte in a header on a LANG=C cluster
    # node must not crash a 350k-proteome read. Sequences are ASCII by construction.
    with _opener(path)(path, "rt", encoding="utf-8", errors="replace") as fh:
        name, buf = None, []
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n\r"), []
            elif line.strip():
                buf.append(line.strip())
        if name is not None:
            yield name, "".join(buf)


def write_fasta(fh, name: str, seq: str, width: int = 60) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def seq_id(header: str) -> str:
    """HMMER's 'target name' is the first whitespace-delimited token."""
    return header.split(None, 1)[0]


def read_stockholm(path: str) -> Tuple[List[str], Dict[str, str], str]:
    """Parse a (possibly multi-block) Stockholm file.

    Returns (name_order, {name: aligned_seq}, RF_string).

    hmmalign writes the reference annotation line `#=GC RF` with 'x' at profile
    match states and '.' at insert states. That line is the whole reason we use
    Stockholm rather than afa here: it gives a profile-anchored coordinate
    system that is identical for every sequence and stable across reruns.
    """
    seqs: Dict[str, List[str]] = {}
    order: List[str] = []
    rf: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line or line == "//" or line.startswith("# STOCKHOLM"):
                continue
            if line.startswith("#=GC RF"):
                rf.append(line.split()[2])
                continue
            if line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            name, chunk = parts
            if name not in seqs:
                seqs[name] = []
                order.append(name)
            seqs[name].append(chunk.strip())
    return order, {k: "".join(v) for k, v in seqs.items()}, "".join(rf)


def match_columns(rf: str) -> List[int]:
    """Indices of profile match states in the full alignment."""
    return [i for i, c in enumerate(rf) if c not in ".-~_"]


def to_match_string(aligned: str, cols: List[int]) -> str:
    """Project a Stockholm-aligned sequence onto match columns, uppercased,
    with all gap characters normalised to '-'."""
    s = "".join(aligned[i] for i in cols).upper()
    return s.replace(".", "-").replace("~", "-")


# ---------------------------------------------------------------------------
# HMMER domtblout
# ---------------------------------------------------------------------------

DOMTBL_COLS = [
    "target_name", "target_accession", "tlen",
    "query_name", "query_accession", "qlen",
    "full_evalue", "full_score", "full_bias",
    "dom_idx", "dom_n", "c_evalue", "i_evalue", "dom_score", "dom_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "acc", "description",
]
_N_FIXED = 22  # fields before the free-text description


def parse_domtblout(path: str) -> Iterator[dict]:
    """Stream a --domtblout file. Whitespace-delimited, 22 fixed fields then a
    free-text description, so split with maxsplit rather than str.split()."""
    with _opener(path)(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split(None, _N_FIXED)
            if len(f) < _N_FIXED:
                continue
            if len(f) == _N_FIXED:
                f.append("")
            yield dict(zip(DOMTBL_COLS, f))


DOMTBL_DTYPES = {
    "tlen": "int32", "qlen": "int32",
    "full_evalue": "float64", "full_score": "float32", "full_bias": "float32",
    "dom_idx": "int16", "dom_n": "int16",
    "c_evalue": "float64", "i_evalue": "float64",
    "dom_score": "float32", "dom_bias": "float32",
    "hmm_from": "int32", "hmm_to": "int32",
    "ali_from": "int32", "ali_to": "int32",
    "env_from": "int32", "env_to": "int32",
    "acc": "float32",
}


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------

def set_style() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        # Keep text as text in the SVG so it stays editable in Illustrator /
        # Inkscape, and embed TrueType (not Type-3) in any PDF.
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "axes.titleweight": "bold",
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.0,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def savefig(fig, figdir: str, name: str, formats=("png", "svg"), dpi: int = 600) -> None:
    import matplotlib.pyplot as plt
    os.makedirs(figdir, exist_ok=True)
    for ext in formats:
        fig.savefig(os.path.join(figdir, f"{name}.{ext}"), dpi=dpi)
    plt.close(fig)


# Colour-blind-safe qualitative palette (Okabe-Ito, extended).
PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00",
    "#56B4E9", "#F0E442", "#000000", "#8C564B", "#7F7F7F",
    "#17BECF", "#9467BD", "#BCBD22", "#AEC7E8", "#FFBB78",
]


def palette(n: int) -> List[str]:
    if n <= len(PALETTE):
        return PALETTE[:n]
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("nipy_spectral")
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


def load_config(path: str) -> dict:
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------------------------

_GTDB_PREFIX = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
                "f": "family", "g": "genus", "s": "species"}
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def _clean_rank(series):
    """Strip the GTDB rank prefix and turn '' / 'p__' into NA."""
    import pandas as pd
    s = series.astype("string").str.strip()
    s = s.str.replace(r"^[a-z]__", "", regex=True)
    return s.replace({"": pd.NA}).astype("string")


def resolve_taxonomy(df, taxonomy_col=None, rank="phylum", min_populated=0.5):
    """Return a Series of taxon labels aligned to `df`, or None.

    Column priority, and it matters:

      1. an explicit `taxonomy_col` from config, used verbatim;
      2. `gtdb_<rank>`;
      3. a bare `<rank>` column -- but ONLY if it is actually populated;
      4. a GTDB lineage string (`classification` / `gtdb_taxonomy` / ...).

    Rule 3's guard is not paranoia. Real metadata tables carry a block of
    host-organism fields (`kingdom`, `phylum`, `class`, `order`, `family`,
    `common_name`) alongside the microbial `gtdb_*` fields. Grabbing the first
    column named `phylum` gets you the host animal's phylum, or, when the block
    is empty, a silently blank taxonomy that propagates into every plot and
    every enrichment test without raising anything.
    """
    import pandas as pd

    if taxonomy_col:
        if taxonomy_col not in df.columns:
            raise KeyError(f"plots.taxonomy_col='{taxonomy_col}' not in the table")
        return _clean_rank(df[taxonomy_col])

    lower = {c.lower(): c for c in df.columns}

    for cand in (f"gtdb_{rank}", rank):
        col = lower.get(cand.lower())
        if col is None:
            continue
        s = _clean_rank(df[col])
        if s.notna().mean() >= min_populated:
            return s
        print(f"[taxonomy] column '{col}' is only {100 * s.notna().mean():.1f}% "
              f"populated; ignoring it and looking for a lineage string",
              file=_sys.stderr)

    want = next((k for k, v in _GTDB_PREFIX.items() if v == rank), None)
    if want:
        for cand in ("classification", "gtdb_taxonomy", "taxonomy", "lineage"):
            col = lower.get(cand)
            if col is None:
                continue

            def _pick(s):
                if not isinstance(s, str):
                    return pd.NA
                for part in s.split(";"):
                    part = part.strip()
                    if part.startswith(f"{want}__"):
                        return part[3:] or pd.NA
                return pd.NA

            s = df[col].astype("string").map(_pick).astype("string")
            if s.notna().mean() >= min_populated:
                return s
    return None


def parse_gtdb_lineage(series, rank):
    """Pull one rank out of a `d__X;p__Y;...` string, keeping the prefix."""
    import pandas as pd
    want = next(k for k, v in _GTDB_PREFIX.items() if v == rank)

    def _pick(s):
        if not isinstance(s, str):
            return pd.NA
        for part in s.split(";"):
            part = part.strip()
            if part.startswith(f"{want}__") and len(part) > 3:
                return part
        return pd.NA

    return series.astype("string").map(_pick).astype("string")


def top_n_with_other(series, n: int, other: str = "Other"):
    import pandas as pd
    counts = series.dropna().value_counts()
    keep = set(counts.index[:n])

    def _m(x):
        if pd.isna(x):
            return "Unclassified"
        return x if x in keep else other

    return series.map(_m)
