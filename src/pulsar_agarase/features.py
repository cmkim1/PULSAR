"""Feature extraction from compact dbCAN/CGCFinder output."""

from __future__ import annotations

import re
import random
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from .model import FAMILIES


FAMILY_RE = re.compile(r"\b(GH(?:2|16|50|86|96|117|118))(?![0-9])")
GH16_AGAR_RE = re.compile(r"\bGH16[_-](?:14|15|16)\b")
ORDERED_ID_RE = re.compile(r"(.+)_([0-9]+)$")
LAHG_RE = re.compile(
    r"(?:\bL-?AHG\b|3,6-anhydro|anhydro-L-galactonate|dau[ABCD]\b|ahg\b)",
    re.IGNORECASE,
)
DETECTION_MARKERS = ["GH16", "GH50", "GH86", "GH96", "GH117", "GH118", "L_AHG"]
SUPPORT_MARKERS = ["GH2"]
DEFAULT_GENE_WINDOWS = [5, 10, 15, 20, 30, 50]
DEFAULT_BP_WINDOWS = [5000, 10000, 20000, 50000, 100000]


def extract_families(text: str) -> set[str]:
    families = {match.group(1) for match in FAMILY_RE.finditer(text)}
    if "GH16" in families and not GH16_AGAR_RE.search(text):
        families.remove("GH16")
    if GH16_AGAR_RE.search(text):
        families.add("GH16")
    return families


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


def _modern_tsv_gene_families(path: Path) -> dict[str, set[str]]:
    gene_families: dict[str, set[str]] = defaultdict(set)
    if not path.exists() or path.stat().st_size == 0:
        return gene_families

    with path.open() as handle:
        header = next(handle, "").rstrip("\n").split("\t")
        normalized = [item.strip().lower().replace("_", " ") for item in header]
        gene_columns = [
            index
            for index, name in enumerate(normalized)
            if name in {"gene id", "protein id", "query id", "query name", "sequence id"}
            or ("gene" in name and "cluster" not in name)
            or ("protein" in name and "family" not in name)
            or name == "query"
        ]

        for line in handle:
            parts = line.rstrip("\n").split("\t")
            families = extract_families(line)
            if not families:
                continue

            candidates = [parts[index] for index in gene_columns if index < len(parts)]
            if not candidates and parts:
                if extract_families(parts[0]) and len(parts) > 2:
                    candidates = [parts[2]]
                else:
                    candidates = [parts[0]]
            if not candidates:
                continue

            gene_families[_normalize_gene_id(candidates[0])].update(families)
    return gene_families


def _count_cgc_families(path: Path) -> tuple[int, Counter[str]]:
    if not path.exists() or path.stat().st_size == 0:
        return 0, Counter()

    loci: dict[str, Counter[str]] = defaultdict(Counter)
    seen = set()
    with path.open() as handle:
        first = next(handle, "")
        header = first.rstrip("\n").split("\t")
        normalized = [item.strip().lower().replace("_", " ") for item in header]
        has_header = any(
            name in {"cgc", "cgc#", "cgc id", "gene id", "protein id", "gene type", "protein family"}
            or name.startswith("cgc ")
            for name in normalized
        )

        cgc_column = 4
        gene_column = 8
        if has_header:
            for index, name in enumerate(normalized):
                if "cgc" in name and ("id" in name or name in {"cgc", "cgc#"}) and "gene" not in name:
                    cgc_column = index
                    break
            for index, name in enumerate(normalized):
                if name in {"gene id", "protein id", "query id", "sequence id"} or (
                    ("gene" in name or "protein" in name) and "type" not in name and "family" not in name
                ):
                    gene_column = index
                    break
            lines = handle
        else:
            lines = [first, *handle]

        for line in lines:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(cgc_column, gene_column):
                continue
            cgc_id = parts[cgc_column]
            gene_id = parts[gene_column]
            families = extract_families(line)
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


def _first_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists() and path.stat().st_size > 0:
            return path
    return paths[0]


def _gff_attributes(text: str) -> dict[str, str]:
    attrs = {}
    for item in text.split(";"):
        item = item.strip()
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def _gff_gene_records(path: Path | None) -> tuple[dict[str, dict[str, object]], dict[str, int], dict[str, set[str]]]:
    records: dict[str, dict[str, object]] = {}
    contig_lengths: dict[str, int] = defaultdict(int)
    lahg_families: dict[str, set[str]] = defaultdict(set)
    if path is None or not path.exists() or path.stat().st_size == 0:
        return records, dict(contig_lengths), lahg_families

    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig = parts[0]
            attrs = _gff_attributes(parts[8])
            ids = [
                attrs.get("ID", ""),
                attrs.get("Name", ""),
                attrs.get("protein_id", ""),
                attrs.get("locus_tag", ""),
            ]
            feature_type = parts[2].lower()
            if feature_type not in {"cds", "gene"} and not any(ids):
                continue
            contig_lengths[contig] += 1
            gene_index = contig_lengths[contig]
            product_text = " ".join([attrs.get("product", ""), attrs.get("Note", ""), parts[8]])
            if LAHG_RE.search(product_text):
                for gene_id in ids:
                    if gene_id:
                        lahg_families[_normalize_gene_id(gene_id)].add("L_AHG")
            for gene_id in ids:
                if not gene_id:
                    continue
                normalized = _normalize_gene_id(gene_id)
                records[normalized] = {
                    "gene_id": normalized,
                    "contig": contig,
                    "gene_index": gene_index,
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "annotation": product_text,
                }
    return records, dict(contig_lengths), lahg_families


def _marker_rows(
    gene_families: dict[str, set[str]],
    gff_records: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows = []
    for gene_id, families in sorted(gene_families.items()):
        location = gff_records.get(gene_id)
        if location is None:
            parsed = _gene_id_order(gene_id)
            if parsed is None:
                continue
            contig, gene_index = parsed
            location = {
                "gene_id": gene_id,
                "contig": contig,
                "gene_index": gene_index,
                "start": gene_index,
                "end": gene_index,
                "strand": ".",
                "annotation": "",
            }
        detection = sorted(family for family in families if family in DETECTION_MARKERS)
        support = sorted(family for family in families if family in SUPPORT_MARKERS)
        if not detection and not support:
            continue
        rows.append(
            {
                **location,
                "markers": ",".join(sorted(detection + support)),
                "detection_markers": ",".join(detection),
                "support_markers": ",".join(support),
                "is_detection_marker": int(bool(detection)),
                "is_support_marker": int(bool(support)),
            }
        )
    return rows


def _infer_contig_lengths(marker_rows: list[dict[str, object]], contig_lengths: dict[str, int]) -> dict[str, int]:
    inferred = dict(contig_lengths)
    for row in marker_rows:
        contig = str(row["contig"])
        gene_index = int(row["gene_index"])
        inferred[contig] = max(inferred.get(contig, 0), gene_index)
    return inferred


def _infer_contig_bp_lengths(marker_rows: list[dict[str, object]], gff_records: dict[str, dict[str, object]]) -> dict[str, int]:
    inferred: dict[str, int] = defaultdict(int)
    for record in gff_records.values():
        inferred[str(record["contig"])] = max(inferred[str(record["contig"])], int(record["end"]))
    for row in marker_rows:
        contig = str(row["contig"])
        inferred[contig] = max(inferred[contig], int(row["end"]))
    return dict(inferred)


def _scan_intervals_for_counts(
    marker_positions_by_contig: dict[str, set[int]],
    contig_lengths: dict[str, int],
    windows: list[int],
    total_markers: int,
    total_genes: int,
    min_markers: int,
) -> list[dict[str, object]]:
    if total_markers == 0 or total_genes == 0:
        return []

    intervals = []
    marker_rate = total_markers / total_genes
    for contig, length in contig_lengths.items():
        if length <= 0:
            continue
        marker_positions = marker_positions_by_contig.get(contig, set())
        counts = [0] * (length + 1)
        for position in marker_positions:
            if 1 <= position <= length:
                counts[position] = 1
        prefix = [0] * (length + 1)
        for index in range(1, length + 1):
            prefix[index] = prefix[index - 1] + counts[index]

        for width in windows:
            if width > length:
                continue
            expected = width * marker_rate
            variance = max(expected * (1 - marker_rate), 1e-9)
            for start in range(1, length - width + 2):
                end = start + width - 1
                observed = prefix[end] - prefix[start - 1]
                if observed < min_markers:
                    continue
                score = (observed - expected) / (variance**0.5)
                if score <= 0:
                    continue
                intervals.append(
                    {
                        "contig": contig,
                        "start_gene": start,
                        "end_gene": end,
                        "n_genes": width,
                        "n_markers": observed,
                        "expected_markers": round(expected, 4),
                        "scan_score": round(score, 4),
                    }
                )
    return intervals


def _scan_intervals_for_sparse_positions(
    marker_positions_by_contig: dict[str, set[int]],
    contig_lengths: dict[str, int],
    windows: list[int],
    total_markers: int,
    total_positions: int,
    min_markers: int,
) -> list[dict[str, object]]:
    if total_markers == 0 or total_positions == 0:
        return []

    intervals = []
    marker_rate = total_markers / total_positions
    for contig, length in contig_lengths.items():
        if length <= 0:
            continue
        positions = sorted(marker_positions_by_contig.get(contig, set()))
        if not positions:
            continue
        for width in windows:
            if width > length:
                continue
            starts = set()
            for position in positions:
                starts.add(max(1, position - width // 2))
                starts.add(max(1, position - width + 1))
                starts.add(position)
            expected = width * marker_rate
            variance = max(expected * (1 - marker_rate), 1e-9)
            for start in sorted(starts):
                end = min(length, start + width - 1)
                observed = sum(1 for position in positions if start <= position <= end)
                if observed < min_markers:
                    continue
                score = (observed - expected) / (variance**0.5)
                if score <= 0:
                    continue
                intervals.append(
                    {
                        "contig": contig,
                        "start_bp": start,
                        "end_bp": end,
                        "n_bp": end - start + 1,
                        "n_markers": observed,
                        "expected_markers": round(expected, 4),
                        "scan_score": round(score, 4),
                        "coordinate_unit": "bp",
                    }
                )
    return intervals


def _overlaps(a: dict[str, object], b: dict[str, object]) -> bool:
    if a["contig"] != b["contig"]:
        return False
    if a.get("coordinate_unit") == "bp" or b.get("coordinate_unit") == "bp":
        return int(a["start_bp"]) <= int(b["end_bp"]) and int(b["start_bp"]) <= int(a["end_bp"])
    return int(a["start_gene"]) <= int(b["end_gene"]) and int(b["start_gene"]) <= int(a["end_gene"])


def _select_nonoverlapping_intervals(intervals: list[dict[str, object]], max_loci: int) -> list[dict[str, object]]:
    selected = []
    ranked = sorted(
        intervals,
        key=lambda row: (-float(row["scan_score"]), -int(row["n_markers"]), int(row.get("n_genes", row.get("n_bp", 0)))),
    )
    for interval in ranked:
        if any(_overlaps(interval, kept) for kept in selected):
            continue
        selected.append(interval)
        if len(selected) >= max_loci:
            break
    return selected


def _random_marker_positions(
    contig_lengths: dict[str, int],
    total_markers: int,
    rng: random.Random,
) -> dict[str, set[int]]:
    total_positions = sum(contig_lengths.values())
    if total_markers >= total_positions:
        sample_numbers = list(range(total_positions))
    else:
        sample_numbers = rng.sample(range(total_positions), total_markers)
    offsets = []
    cursor = 0
    for contig, length in contig_lengths.items():
        offsets.append((cursor, cursor + length, contig))
        cursor += length
    positions: dict[str, set[int]] = defaultdict(set)
    for number in sample_numbers:
        for start, end, contig in offsets:
            if start <= number < end:
                positions[contig].add(number - start + 1)
                break
    return positions


def detect_agar_puls_by_scan(
    marker_rows: list[dict[str, object]],
    contig_lengths: dict[str, int],
    contig_bp_lengths: dict[str, int] | None = None,
    windows: list[int] | None = None,
    scan_unit: str = "gene",
    permutations: int = 999,
    seed: int = 1,
    min_markers: int = 2,
    alpha: float = 0.05,
    max_loci: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if scan_unit not in {"gene", "bp"}:
        raise ValueError("scan_unit must be 'gene' or 'bp'")
    windows = windows or (DEFAULT_BP_WINDOWS if scan_unit == "bp" else DEFAULT_GENE_WINDOWS)
    contig_lengths = _infer_contig_lengths(marker_rows, contig_lengths)
    scan_lengths = contig_bp_lengths if scan_unit == "bp" and contig_bp_lengths else contig_lengths
    detection_rows = [row for row in marker_rows if int(row["is_detection_marker"]) == 1]
    total_positions = sum(scan_lengths.values())
    marker_positions_by_contig: dict[str, set[int]] = defaultdict(set)
    for row in detection_rows:
        if scan_unit == "bp":
            marker_positions_by_contig[str(row["contig"])].add((int(row["start"]) + int(row["end"])) // 2)
        else:
            marker_positions_by_contig[str(row["contig"])].add(int(row["gene_index"]))
    total_markers = sum(len(values) for values in marker_positions_by_contig.values())

    if scan_unit == "bp":
        intervals = _scan_intervals_for_sparse_positions(
            marker_positions_by_contig=marker_positions_by_contig,
            contig_lengths=scan_lengths,
            windows=windows,
            total_markers=total_markers,
            total_positions=total_positions,
            min_markers=min_markers,
        )
    else:
        intervals = _scan_intervals_for_counts(
            marker_positions_by_contig=marker_positions_by_contig,
            contig_lengths=scan_lengths,
            windows=windows,
            total_markers=total_markers,
            total_genes=total_positions,
            min_markers=min_markers,
        )
        for interval in intervals:
            interval["coordinate_unit"] = "gene"
    selected = _select_nonoverlapping_intervals(intervals, max_loci=max_loci)

    if permutations > 0 and selected:
        rng = random.Random(seed)
        rank_scores = [[] for _ in selected]
        for _ in range(permutations):
            random_positions = _random_marker_positions(scan_lengths, total_markers, rng)
            if scan_unit == "bp":
                random_intervals = _scan_intervals_for_sparse_positions(
                    marker_positions_by_contig=random_positions,
                    contig_lengths=scan_lengths,
                    windows=windows,
                    total_markers=total_markers,
                    total_positions=total_positions,
                    min_markers=min_markers,
                )
            else:
                random_intervals = _scan_intervals_for_counts(
                    marker_positions_by_contig=random_positions,
                    contig_lengths=scan_lengths,
                    windows=windows,
                    total_markers=total_markers,
                    total_genes=total_positions,
                    min_markers=min_markers,
                )
                for interval in random_intervals:
                    interval["coordinate_unit"] = "gene"
            random_selected = _select_nonoverlapping_intervals(random_intervals, max_loci=len(selected))
            for rank_index in range(len(selected)):
                score = float(random_selected[rank_index]["scan_score"]) if rank_index < len(random_selected) else 0.0
                rank_scores[rank_index].append(score)
        for rank_index, interval in enumerate(selected):
            observed_score = float(interval["scan_score"])
            exceed = sum(score >= observed_score for score in rank_scores[rank_index])
            interval["empirical_p"] = round((exceed + 1) / (permutations + 1), 6)
            interval["rank"] = rank_index + 1
    else:
        for rank_index, interval in enumerate(selected):
            interval["empirical_p"] = 0.0
            interval["rank"] = rank_index + 1

    significant = [row for row in selected if float(row.get("empirical_p", 1.0)) <= alpha]
    for index, interval in enumerate(significant, start=1):
        interval["pul_id"] = f"PUL_{index}"

    candidates = pd.DataFrame(intervals).sort_values(["scan_score", "n_markers"], ascending=[False, False]) if intervals else pd.DataFrame()
    puls = pd.DataFrame(significant)
    return candidates, puls


def _counts_in_scan_puls(marker_rows: list[dict[str, object]], puls: pd.DataFrame, scan_unit: str = "gene") -> Counter[str]:
    counts: Counter[str] = Counter()
    if puls.empty:
        return counts
    for row in marker_rows:
        contig = str(row["contig"])
        gene_index = int(row["gene_index"])
        midpoint = (int(row["start"]) + int(row["end"])) // 2
        in_pul = False
        for _, pul in puls.iterrows():
            if scan_unit == "bp":
                in_pul = contig == str(pul["contig"]) and int(pul["start_bp"]) <= midpoint <= int(pul["end_bp"])
            else:
                in_pul = contig == str(pul["contig"]) and int(pul["start_gene"]) <= gene_index <= int(pul["end_gene"])
            if in_pul:
                break
        if not in_pul:
            continue
        for family in str(row.get("markers", "")).split(","):
            if family in FAMILIES:
                counts[family] += 1
    return counts


def _gene_id_order(gene_id: str) -> tuple[str, int] | None:
    match = ORDERED_ID_RE.match(gene_id)
    if not match:
        return None
    return match.group(1), int(match.group(2))

def dbcan_feature_tables(
    dbcan_dir: Path,
    genome: str | None = None,
    taxname: str | None = None,
    gff_path: Path | None = None,
    scan_windows: list[int] | None = None,
    scan_permutations: int = 999,
    scan_seed: int = 1,
    scan_alpha: float = 0.05,
    scan_unit: str = "gene",
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genome_id = genome or dbcan_dir.name
    cgc_loci, strict_counts = _count_cgc_families(
        _first_existing(
            [
                dbcan_dir / "cgc.out",
                dbcan_dir / "cgc_standard_out.tsv",
                dbcan_dir / "total_cgc_info.tsv",
            ]
        )
    )
    gff_for_order = gff_path or _first_existing([dbcan_dir / "cgc.gff", dbcan_dir / "uniInput.gff"])
    gff_records, contig_lengths, lahg_families = _gff_gene_records(gff_for_order)
    gene_families = _merge_gene_family_maps(
        _hmmer_gene_families(dbcan_dir / "hmmer.out"),
        _diamond_gene_families(dbcan_dir / "diamond.out"),
        _modern_tsv_gene_families(dbcan_dir / "dbCAN_hmm_results.tsv"),
        _modern_tsv_gene_families(dbcan_dir / "dbCANsub_hmm_results.tsv"),
        _modern_tsv_gene_families(dbcan_dir / "overview.tsv"),
        lahg_families,
    )
    marker_rows = _marker_rows(gene_families, gff_records)
    marker_table = pd.DataFrame(marker_rows)
    contig_bp_lengths = _infer_contig_bp_lengths(marker_rows, gff_records)
    candidate_windows, scan_puls = detect_agar_puls_by_scan(
        marker_rows=marker_rows,
        contig_lengths=contig_lengths,
        contig_bp_lengths=contig_bp_lengths,
        windows=scan_windows,
        scan_unit=scan_unit,
        permutations=scan_permutations,
        seed=scan_seed,
        alpha=scan_alpha,
    )
    genome_counts = _count_gene_family_map(gene_families)
    scan_counts = _counts_in_scan_puls(marker_rows, scan_puls, scan_unit=scan_unit)
    broad_loci = max(cgc_loci, len(scan_puls))

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
        broad = max(strict, int(scan_counts[family]))
        genome_total = max(int(genome_counts[family]), strict) if has_genome_wide else strict
        outside = max(0, genome_total - strict)
        row[f"strict_n_{family}"] = strict
        row[f"genome_n_{family}"] = genome_total
        row[f"broad_locus_n_{family}"] = broad
        row[f"outside_strict_n_{family}"] = outside
        row[f"outside_strict_has_{family}"] = int(outside > 0)
    row["scan_n_candidate_windows"] = int(len(candidate_windows))
    row["scan_n_significant_puls"] = int(len(scan_puls))
    row["scan_min_empirical_p"] = float(scan_puls["empirical_p"].min()) if not scan_puls.empty else 1.0
    return row, marker_table, candidate_windows, scan_puls


def features_from_dbcan_dir(
    dbcan_dir: Path,
    genome: str | None = None,
    taxname: str | None = None,
    gff_path: Path | None = None,
    scan_windows: list[int] | None = None,
    scan_permutations: int = 999,
    scan_seed: int = 1,
    scan_alpha: float = 0.05,
    scan_unit: str = "gene",
) -> dict[str, object]:
    row, _, _, _ = dbcan_feature_tables(
        dbcan_dir=dbcan_dir,
        genome=genome,
        taxname=taxname,
        gff_path=gff_path,
        scan_windows=scan_windows,
        scan_permutations=scan_permutations,
        scan_seed=scan_seed,
        scan_alpha=scan_alpha,
        scan_unit=scan_unit,
    )
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
