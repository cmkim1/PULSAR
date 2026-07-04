"""PULSAR architecture-based agarase/PUL scoring model."""

from __future__ import annotations

from typing import Any

import pandas as pd


FAMILIES = ["GH2", "GH16", "GH50", "GH86", "GH96", "GH117", "GH118"]
CORE_OPENERS = ["GH16", "GH86", "GH118"]
AUX_OPENERS = ["GH96", "GH2"]


def _val(row: pd.Series, col: str, default: float = 0.0) -> float:
    value = row.get(col, default)
    if pd.isna(value):
        return default
    return float(value)


def _counts(row: pd.Series, prefix: str) -> dict[str, int]:
    return {family: int(_val(row, f"{prefix}_{family}")) for family in FAMILIES}


def _add(scores: dict[str, float], reasons: dict[str, list[str]], family: str, amount: float, reason: str) -> None:
    scores[family] += amount
    reasons[family].append(reason)


def score_row(row: pd.Series) -> dict[str, Any]:
    """Score one genome-level agar-PUL feature row.

    The model does not rank GH16, GH86, and GH118 a priori. These core opener
    families are scored from locus architecture: absent from the GH117-centered
    strict PUL, present outside the strict PUL, or missing as an unresolved
    core-opener group.
    """

    scores = {family: 0.0 for family in FAMILIES}
    reasons = {family: [] for family in FAMILIES}

    strict_loci = int(_val(row, "strict_n_agar_loci"))
    broad_loci = int(_val(row, "broad_n_agar_loci", strict_loci))
    has_global = int(_val(row, "has_genome_wide_annotation"))

    if strict_loci == 0 and broad_loci == 0:
        return {
            **{f"{family}_score": 0.0 for family in FAMILIES},
            "recommended_GH_group": "none",
            "top_recommended_GH": "none",
            "prediction_class": "no_detected_agar_PUL",
            "model_confidence": "low",
            "model_rationale": "No strict or broad agar-PUL context was detected.",
            "central_pul_status": "no_detected_agar_PUL",
            "core_pathway_status": "not_detected",
            "auxiliary_pathway_status": "not_detected",
            "core_opener_status": "not_detected",
        }

    context_prefix = "strict_n" if strict_loci > 0 else "broad_locus_n"
    locus_count = strict_loci if strict_loci > 0 else broad_loci
    context_label = "strict PUL" if strict_loci > 0 else "broad locus"

    strict = _counts(row, context_prefix)
    genome = _counts(row, "genome_n")
    outside = _counts(row, "outside_strict_n")

    strict_core_openers = sum(strict[family] for family in CORE_OPENERS)
    genome_core_openers = sum(genome[family] for family in CORE_OPENERS)
    strict_aux_openers = sum(strict[family] for family in AUX_OPENERS)
    central = strict["GH117"] > 0

    central_status = f"GH117_positive_{context_label.replace(' ', '_')}" if central else f"no_GH117_positive_{context_label.replace(' ', '_')}"
    if not central and (strict_core_openers + strict["GH50"] + strict_aux_openers) > 0:
        _add(scores, reasons, "GH117", 4.0, f"GH117-positive central context is absent from {context_label}")

    if strict_core_openers > 0 and strict["GH50"] > 0 and strict["GH117"] > 0:
        core_status = "complete_detected_core"
    elif strict["GH50"] > 0 and strict["GH117"] > 0:
        core_status = "core_missing_polysaccharide_opener"
    elif strict_core_openers > 0 and strict["GH117"] > 0:
        core_status = "core_missing_GH50"
    elif strict_core_openers > 0 or strict["GH50"] > 0:
        core_status = "partial_core"
    else:
        core_status = "not_detected"

    if central and strict["GH50"] > 0 and strict_core_openers == 0:
        if genome_core_openers == 0:
            for family in CORE_OPENERS:
                _add(
                    scores,
                    reasons,
                    family,
                    3.0,
                    f"GH117/GH50 {context_label} lacks any detected GH16/GH86/GH118 core opener",
                )
            core_opener_status = "unresolved_core_opener_missing"
        else:
            candidates = []
            for family in CORE_OPENERS:
                if strict[family] == 0 and outside[family] > 0:
                    _add(
                        scores,
                        reasons,
                        family,
                        3.0,
                        f"{family} is present genome-wide but absent from the GH117-centered {context_label}",
                    )
                    candidates.append(family)
            core_opener_status = "outside_core_opener_candidate" if candidates else "strict_core_opener_gap"
    elif central and strict["GH50"] > 0:
        candidates = []
        for family in CORE_OPENERS:
            if strict[family] == 0 and outside[family] > 0:
                _add(
                    scores,
                    reasons,
                    family,
                    3.0,
                    f"{family} is present genome-wide but absent from the GH117-centered {context_label}",
                )
                candidates.append(family)
        core_opener_status = "outside_core_opener_candidate" if candidates else "strict_core_opener_present"
    else:
        core_opener_status = "not_GH117_GH50_core_context"

    if strict["GH50"] == 0 and strict["GH117"] > 0 and genome_core_openers > 0:
        _add(scores, reasons, "GH50", 3.5, f"Core pathway has GH16/GH86/GH118 context and GH117 but lacks {context_label} GH50")

    terminal_substrate_context = strict["GH50"] + strict["GH2"]
    if broad_loci >= 3 and strict["GH117"] <= 2 and terminal_substrate_context >= 6:
        _add(scores, reasons, "GH117", 3.0, f"GH117 is limited relative to {context_label} GH50/GH2 terminal-substrate context")
    elif strict["GH117"] == 1 and locus_count >= 2:
        _add(scores, reasons, "GH117", 1.5, f"GH117 is single-copy across multiple {context_label}s")

    if strict_aux_openers > 0 and strict["GH117"] > 0:
        aux_status = "complete_or_overlapping_auxiliary"
    elif strict_aux_openers > 0:
        aux_status = "auxiliary_missing_GH117"
        _add(scores, reasons, "GH117", 3.0, f"auxiliary GH96/GH2 context lacks GH117 in {context_label}")
    elif strict["GH117"] > 0:
        aux_status = f"auxiliary_openers_outside_{context_label.replace(' ', '_')}"
    else:
        aux_status = "not_detected"

    if strict["GH96"] == 0 and strict["GH2"] > 0 and strict["GH117"] > 0:
        _add(scores, reasons, "GH96", 1.5, f"auxiliary pathway has GH2/GH117 but no {context_label} GH96")
    if strict["GH2"] == 0 and strict["GH96"] > 0 and strict["GH117"] > 0:
        _add(scores, reasons, "GH2", 1.0, f"auxiliary pathway has GH96/GH117 but no {context_label} GH2")

    threshold = 3.0
    recommended = [family for family in FAMILIES if scores[family] >= threshold]
    recommended = sorted(recommended, key=lambda family: (-scores[family], FAMILIES.index(family)))[:3]
    recommended = sorted(recommended, key=FAMILIES.index)

    if recommended:
        recommendation = "+".join(recommended)
        prediction_class = "likely"
    else:
        top = max(FAMILIES, key=lambda family: scores[family])
        if scores[top] >= 2.0:
            recommended = [top]
            recommendation = f"{top}_possible"
            prediction_class = "possible"
        else:
            recommendation = "none_low_priority"
            prediction_class = "low_priority"

    top = recommended[0] if recommended else "none"
    max_score = max(scores.values())
    if has_global and max_score >= threshold:
        confidence = "medium_high"
    elif max_score >= threshold:
        confidence = "medium"
    else:
        confidence = "low"

    rationale_bits = []
    for family in recommended or [max(FAMILIES, key=lambda item: scores[item])]:
        if reasons[family]:
            rationale_bits.append(f"{family}: " + "; ".join(dict.fromkeys(reasons[family])))
    if not rationale_bits:
        rationale_bits.append("No strong pathway-centered GH bottleneck was detected.")

    return {
        **{f"{family}_score": round(scores[family], 2) for family in FAMILIES},
        "recommended_GH_group": recommendation,
        "top_recommended_GH": top,
        "prediction_class": prediction_class,
        "model_confidence": confidence,
        "model_rationale": " | ".join(rationale_bits),
        "central_pul_status": central_status,
        "core_pathway_status": core_status,
        "auxiliary_pathway_status": aux_status,
        "core_opener_status": core_opener_status,
    }


def score_dataframe(features: pd.DataFrame) -> pd.DataFrame:
    """Return input features plus score columns."""

    scored = pd.DataFrame([score_row(row) for _, row in features.iterrows()])
    return pd.concat([features.reset_index(drop=True), scored], axis=1)
