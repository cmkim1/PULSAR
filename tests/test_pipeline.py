from pulsar_agarase.pipeline import dbcan_tools_args, sanitize_gff_for_dbcan


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
