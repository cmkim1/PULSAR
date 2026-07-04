"""Feature extraction from compact dbCAN/CGCFinder output."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from .model import FAMILIES


FAMILY_RE = re.compile(r"\b(GH(?:2|16|50|86|96|117|118))(?:[_\b;,\s]|$)")


def extract_families(text: str) -> set[str]:
    return {match.group(1) for match in FAMILY_RE.finditer(text)}


def _count_hmmer_families(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    if not path.exists() or path.stat().st_size == 0:
        return counts

    seen_gene_family = set()
    with path.open() as handle:
        header = next(handle, "")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            profile, gene_id = parts[0], parts[2]
            for family in extract_families(profile):
                key = (gene_id, family)
                if key not in seen_gene_family:
                    counts[family] += 1
                    seen_gene_family.add(key)
    return counts


def _count_cgc_families(path: Path) -> tuple[int, Counter[str]]:
    if not path.exists() or path.stat().st_size == 0:
        return 0, Counter()

    loci: dict[str, Counter[str]] = defaultdict(Counter)
    seen = set()
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            cgc_id = parts[4]
            gene_id = parts[8]
            families = extract_families("\t".join(parts[10:]))
            for family in families:
                key = (cgc_id, gene_id, family)
                if key not in seen:
                    loci[cgc_id][family] += 1
                    seen.add(key)

    agar_loci = {cgc_id: counts for cgc_id, counts in loci.items() if sum(counts[family] for family in FAMILIES) > 0}
    total = Counter()
    for counts in agar_loci.values():
        total.update(counts)
    return len(agar_loci), total


def features_from_dbcan_dir(dbcan_dir: Path, genome: str | None = None, taxname: str | None = None) -> dict[str, object]:
    genome_id = genome or dbcan_dir.name
    cgc_loci, strict_counts = _count_cgc_families(dbcan_dir / "cgc.out")
    genome_counts = _count_hmmer_families(dbcan_dir / "hmmer.out")

    has_genome_wide = int(bool(genome_counts))
    row: dict[str, object] = {
        "genome": genome_id,
        "taxname": taxname or genome_id,
        "strict_n_agar_loci": cgc_loci,
        "broad_n_agar_loci": cgc_loci,
        "has_genome_wide_annotation": has_genome_wide,
    }
    for family in FAMILIES:
        strict = int(strict_counts[family])
        genome_total = max(int(genome_counts[family]), strict) if has_genome_wide else strict
        outside = max(0, genome_total - strict)
        row[f"strict_n_{family}"] = strict
        row[f"genome_n_{family}"] = genome_total
        row[f"broad_locus_n_{family}"] = strict
        row[f"outside_strict_n_{family}"] = outside
        row[f"outside_strict_has_{family}"] = int(outside > 0)
    return row


def features_from_dbcan_root(root: Path, metadata: Path | None = None) -> pd.DataFrame:
    names: dict[str, str] = {}
    if metadata is not None:
        meta = pd.read_csv(metadata, sep="\t")
        if "genome" not in meta.columns or "taxname" not in meta.columns:
            raise ValueError("metadata must contain genome and taxname columns")
        names = dict(zip(meta["genome"].astype(str), meta["taxname"].astype(str)))

    rows = []
    for child in sorted(root.iterdir()):
        if child.is_dir():
            rows.append(features_from_dbcan_dir(child, genome=child.name, taxname=names.get(child.name)))
    return pd.DataFrame(rows)
