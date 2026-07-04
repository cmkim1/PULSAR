import pandas as pd

from pulsar_agarase.model import score_dataframe


def test_unresolved_core_openers_are_not_ranked():
    row = {
        "genome": "g1",
        "taxname": "Example",
        "strict_n_agar_loci": 1,
        "broad_n_agar_loci": 1,
        "has_genome_wide_annotation": 1,
        "strict_n_GH2": 0,
        "genome_n_GH2": 0,
        "broad_locus_n_GH2": 0,
        "outside_strict_n_GH2": 0,
        "strict_n_GH16": 0,
        "genome_n_GH16": 0,
        "broad_locus_n_GH16": 0,
        "outside_strict_n_GH16": 0,
        "strict_n_GH50": 2,
        "genome_n_GH50": 2,
        "broad_locus_n_GH50": 2,
        "outside_strict_n_GH50": 0,
        "strict_n_GH86": 0,
        "genome_n_GH86": 0,
        "broad_locus_n_GH86": 0,
        "outside_strict_n_GH86": 0,
        "strict_n_GH96": 0,
        "genome_n_GH96": 0,
        "broad_locus_n_GH96": 0,
        "outside_strict_n_GH96": 0,
        "strict_n_GH117": 1,
        "genome_n_GH117": 1,
        "broad_locus_n_GH117": 1,
        "outside_strict_n_GH117": 0,
        "strict_n_GH118": 0,
        "genome_n_GH118": 0,
        "broad_locus_n_GH118": 0,
        "outside_strict_n_GH118": 0,
    }
    scored = score_dataframe(pd.DataFrame([row]))
    assert scored.loc[0, "recommended_GH_group"] == "GH16+GH86+GH118"
    assert scored.loc[0, "core_opener_status"] == "unresolved_core_opener_missing"


def test_broad_locus_context_is_scored_when_strict_cgc_is_empty():
    row = {
        "genome": "g1",
        "taxname": "Example",
        "strict_n_agar_loci": 0,
        "broad_n_agar_loci": 1,
        "has_genome_wide_annotation": 1,
    }
    for family in ["GH2", "GH16", "GH50", "GH86", "GH96", "GH117", "GH118"]:
        row[f"strict_n_{family}"] = 0
        row[f"genome_n_{family}"] = 0
        row[f"broad_locus_n_{family}"] = 0
        row[f"outside_strict_n_{family}"] = 0
    row["genome_n_GH50"] = 1
    row["genome_n_GH117"] = 1
    row["broad_locus_n_GH50"] = 1
    row["broad_locus_n_GH117"] = 1

    scored = score_dataframe(pd.DataFrame([row]))

    assert scored.loc[0, "recommended_GH_group"] == "GH16+GH86+GH118"
    assert scored.loc[0, "central_pul_status"] == "GH117_positive_broad_locus"
