from pulsar_agarase.pipeline import dbcan_database_present, dbcan_tools_args, sanitize_gff_for_dbcan


def test_modern_dbcan_tools_arg_keeps_comma_list():
    assert dbcan_tools_args("hmmer,diamond", legacy_script=False) == ["hmmer,diamond"]


def test_legacy_dbcan_tools_arg_splits_comma_list():
    assert dbcan_tools_args("hmmer,diamond", legacy_script=True) == ["hmmer", "diamond"]


def test_legacy_dbcan_tools_arg_accepts_space_list():
    assert dbcan_tools_args("hmmer diamond", legacy_script=True) == ["hmmer", "diamond"]


def test_sanitize_gff_for_dbcan_removes_empty_attributes(tmp_path):
    src = tmp_path / "input.gff"
    dst = tmp_path / "output.gff"
    src.write_text("ctg\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1;partial=00;\n")

    sanitize_gff_for_dbcan(src, dst)

    assert dst.read_text() == "ctg\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1;partial=00\n"


def test_modern_dbcan_database_markers_are_detected(tmp_path):
    for marker in ["CAZy.dmnd", "dbCAN.hmm", "dbCAN-sub.hmm", "TCDB.dmnd", "TF.hmm", "STP.hmm"]:
        (tmp_path / marker).write_text("")

    assert dbcan_database_present(tmp_path)


def test_partial_dbcan_download_is_not_detected_as_complete(tmp_path):
    for marker in ["CAZy.dmnd", "dbCAN.hmm", "TCDB.dmnd", "TF.hmm", "STP.hmm"]:
        (tmp_path / marker).write_text("")
    (tmp_path / "dbCAN-sub.hmm.part").write_text("")

    assert not dbcan_database_present(tmp_path)
