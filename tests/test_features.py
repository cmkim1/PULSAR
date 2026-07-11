from pathlib import Path

from pulsar_agarase.features import features_from_dbcan_dir


def test_features_use_diamond_and_gff_for_broad_locus(tmp_path: Path):
    dbcan = tmp_path / "dbcan"
    dbcan.mkdir()
    (dbcan / "hmmer.out").write_text("header\n")
    (dbcan / "cgc.out").write_text("")
    (dbcan / "diamond.out").write_text(
        "Gene ID\tCAZy ID\t% Identical\n"
        "1_50\tAAA|GH50|\t90\n"
        "1_52\tBBB|GH117|\t85\n"
        "1_53\tCCC|GH86|\t80\n"
    )
    gff_lines = []
    for index in range(1, 501):
        gff_lines.append(f"ctg\tProdigal\tCDS\t{index * 100}\t{index * 100 + 89}\t.\t+\t0\tID=1_{index}\n")
    (dbcan / "cgc.gff").write_text("".join(gff_lines))

    row = features_from_dbcan_dir(dbcan, genome="g1", scan_permutations=99)

    assert row["has_genome_wide_annotation"] == 1
    assert row["genome_n_GH50"] == 1
    assert row["genome_n_GH117"] == 1
    assert row["strict_n_agar_loci"] == 0
    assert row["broad_n_agar_loci"] == 1
    assert row["broad_locus_n_GH50"] == 1
    assert row["broad_locus_n_GH86"] == 1
    assert row["broad_locus_n_GH117"] == 1


def test_features_read_modern_dbcan_easy_cgc_outputs(tmp_path: Path):
    dbcan = tmp_path / "dbcan"
    dbcan.mkdir()
    (dbcan / "cgc_standard_out.tsv").write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tDirection\tProtein Family\n"
        "CGC1\tCAZyme\tctg\tctg_1\t1\t90\t+\tGH50.hmm\n"
        "CGC1\tCAZyme\tctg\tctg_2\t91\t180\t+\tGH117_e1\n"
        "CGC1\tCAZyme\tctg\tctg_3\t181\t270\t+\tGH2\n"
    )
    (dbcan / "dbCAN_hmm_results.tsv").write_text(
        "HMM Profile\tProfile Length\tGene ID\tGene Length\tE Value\n"
        "GH50.hmm\t300\tctg_1\t350\t1e-50\n"
        "GH117.hmm\t400\tctg_2\t420\t1e-60\n"
    )
    (dbcan / "overview.tsv").write_text(
        "Gene ID\tHMMER\tdbCAN_sub\tDIAMOND\n"
        "ctg_3\tGH2.hmm\t-\t-\n"
    )
    (dbcan / "cgc.gff").write_text(
        "ctg\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=ctg_1\n"
        "ctg\tProdigal\tCDS\t91\t180\t.\t+\t0\tID=ctg_2\n"
        "ctg\tProdigal\tCDS\t181\t270\t.\t+\t0\tID=ctg_3\n"
    )

    row = features_from_dbcan_dir(dbcan, genome="modern")

    assert row["strict_n_agar_loci"] == 1
    assert row["broad_n_agar_loci"] == 1
    assert row["strict_n_GH50"] == 1
    assert row["strict_n_GH117"] == 1
    assert row["strict_n_GH2"] == 1
    assert row["genome_n_GH2"] == 1


def test_features_infer_broad_locus_from_ordered_gene_ids_without_cgc_gff(tmp_path: Path):
    dbcan = tmp_path / "dbcan"
    dbcan.mkdir()
    (dbcan / "cgc_standard_out.tsv").write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n"
    )
    (dbcan / "overview.tsv").write_text(
        "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n"
        "NC_016613.1_100\t-\t-\t-\tGH50\t1\t-\t-\n"
        "NC_016613.1_101\t-\t-\t-\tGH117\t1\t-\t-\n"
        "NC_016613.1_102\t-\t-\t-\tGH2\t1\t-\t-\n"
        "NC_016613.1_103\t-\t-\t-\tGH86\t1\t-\t-\n"
        "NC_016613.1_1000\t-\t-\t-\tGH2\t1\t-\t-\n"
    )

    row = features_from_dbcan_dir(dbcan, genome="ordered", scan_permutations=99)

    assert row["strict_n_agar_loci"] == 0
    assert row["broad_n_agar_loci"] == 1
    assert row["broad_locus_n_GH50"] == 1
    assert row["broad_locus_n_GH86"] == 1
    assert row["broad_locus_n_GH117"] == 1
    assert row["broad_locus_n_GH2"] == 1
