from pathlib import Path

from pulsar_agarase.features import features_from_dbcan_dir


def test_features_use_diamond_and_gff_for_broad_locus(tmp_path: Path):
    dbcan = tmp_path / "dbcan"
    dbcan.mkdir()
    (dbcan / "hmmer.out").write_text("header\n")
    (dbcan / "cgc.out").write_text("")
    (dbcan / "diamond.out").write_text(
        "Gene ID\tCAZy ID\t% Identical\n"
        "1_1\tAAA|GH50|\t90\n"
        "1_5\tBBB|GH117|\t85\n"
    )
    (dbcan / "cgc.gff").write_text(
        "ctg\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1\n"
        "ctg\tProdigal\tCDS\t91\t180\t.\t+\t0\tID=1_2\n"
        "ctg\tProdigal\tCDS\t181\t270\t.\t+\t0\tID=1_3\n"
        "ctg\tProdigal\tCDS\t271\t360\t.\t+\t0\tID=1_4\n"
        "ctg\tProdigal\tCDS\t361\t450\t.\t+\t0\tID=1_5\n"
    )

    row = features_from_dbcan_dir(dbcan, genome="g1")

    assert row["has_genome_wide_annotation"] == 1
    assert row["genome_n_GH50"] == 1
    assert row["genome_n_GH117"] == 1
    assert row["strict_n_agar_loci"] == 0
    assert row["broad_n_agar_loci"] == 1
    assert row["broad_locus_n_GH50"] == 1
    assert row["broad_locus_n_GH117"] == 1
