"""Feature extraction from compact dbCAN/CGCFinder output."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from .model import FAMILIES


FAMILY_RE = re.compile(r"\b(GH(?:2|16|50|86|96|117|118))(?:[_|;,\s]|$)")
ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")


def extract_families(text: str) -> set[str]:
    return {match.group(1) for match in FAMILY_RE.finditer(text)}


def _normalize_gene_id(gene_id: str) -> str:
    gene_id = gene_id.strip()
    if "|" in gene_id:
        gene_id = gene_id.split("|")[-1]
    if "_prot_" in gene_id:
        gene_id = gene_id.split("_prot_")[-1]
    match = re.search(r"(?:WP|YP|NP)_[0-9]+\.[0-9]+", gene_id)
    if match:
        return match.group(0)
    return gene_id.split()[0]


def _merge_gene_family_maps(*maps: dict[str, set[str]]) -> dict[str, set[str]]:
    merged: dict[str, set[str]] = defaultdict(set)
    for gene_map in maps:
        for gene_id, families in gene_map.items():
            merged[_normalize_gene_id(gene_id)].update(families)
    return merged


def _count_gene_family_map(gene_families: dict[str, set[str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for families in gene_families.values():
        for family in families:
            counts[family] += 1
    return counts


def _hmmer_gene_families(path: Path) -> dict[str, set[str]]:
    gene_families: dict[str, set[str]] = defaultdict(set)
    if not path.exists() or path.stat().st_size == 0:
        return gene_families

    with path.open() as handle:
        next(handle, "")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            profile, gene_id = parts[0], parts[2]
            for family in extract_families(profile):
                gene_families[_normalize_gene_id(gene_id)].add(family)
    return gene_families


def _diamond_gene_families(path: Path) -> dict[str, set[str]]:
    gene_families: dict[str, set[str]] = defaultdict(set)
    if not path.exists() or path.stat().st_size == 0:
        return gene_families

    with path.open() as handle:
        next(handle, "")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene_id, cazy_id = parts[0], parts[1]
            for family in extract_families(cazy_id):
                gene_families[_normalize_gene_id(gene_id)].add(family)
    return gene_families


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


def _gff_gene_order(path: Path) -> list[tuple[str, str]]:
    order = []
    if not path.exists() or path.stat().st_size == 0:
        return order

    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            match = ID_RE.search(parts[8])
            if not match:
                continue
            order.append((parts[0], _normalize_gene_id(match.group(1))))
    return order


def _fallback_broad_loci(
    gff_path: Path,
    gene_families: dict[str, set[str]],
    max_gap: int = 10,
) -> tuple[int, Counter[str]]:
    gene_order = _gff_gene_order(gff_path)
    marker_positions = [
        (contig, index, gene_id)
        for index, (contig, gene_id) in enumerate(gene_order)
        if gene_families.get(gene_id)
    ]
    if not marker_positions:
        return 0, Counter()

    loci: list[list[str]] = []
    current = [marker_positions[0][2]]
    prev_contig, prev_index, _ = marker_positions[0]
    for contig, index, gene_id in marker_positions[1:]:
        if contig == prev_contig and index - prev_index <= max_gap:
            current.append(gene_id)
        else:
            loci.append(current)
            current = [gene_id]
        prev_contig, prev_index = contig, index
    loci.append(current)

    agar_loci = [locus for locus in loci if len(locus) >= 2]
    counts: Counter[str] = Counter()
    for locus in agar_loci:
        seen = set()
        for gene_id in locus:
            for family in gene_families.get(gene_id, set()):
                key = (gene_id, family)
                if key not in seen:
                    counts[family] += 1
                    seen.add(key)
    return len(agar_loci), counts


def features_from_dbcan_dir(dbcan_dir: Path, genome: str | None = None, taxname: str | None = None) -> dict[str, object]:
    genome_id = genome or dbcan_dir.name
    cgc_loci, strict_counts = _count_cgc_families(dbcan_dir / "cgc.out")
    gene_families = _merge_gene_family_maps(
        _hmmer_gene_families(dbcan_dir / "hmmer.out"),
        _diamond_gene_families(dbcan_dir / "diamond.out"),
    )
    genome_counts = _count_gene_family_map(gene_families)
    fallback_loci, fallback_counts = _fallback_broad_loci(dbcan_dir / "cgc.gff", gene_families)
    broad_loci = max(cgc_loci, fallback_loci)

    has_genome_wide = int(bool(genome_counts))
    row: dict[str, object] = {
        "genome": genome_id,
        "taxname": taxname or genome_id,
        "strict_n_agar_loci": cgc_loci,
        "broad_n_agar_loci": broad_loci,
        "has_genome_wide_annotation": has_genome_wide,
    }
    for family in FAMILIES:
        strict = int(strict_counts[family])
        broad = max(strict, int(fallback_counts[family]))
        genome_total = max(int(genome_counts[family]), strict) if has_genome_wide else strict
        outside = max(0, genome_total - strict)
        row[f"strict_n_{family}"] = strict
        row[f"genome_n_{family}"] = genome_total
        row[f"broad_locus_n_{family}"] = broad
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
