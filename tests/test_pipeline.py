from pathlib import Path

from pulsar_agarase.pipeline import (
    build_modern_run_dbcan_command,
    dbcan_database_present,
    dbcan_tools_args,
    sanitize_gff_for_dbcan,
)


def test_modern_dbcan_tools_arg_keeps_comma_list():
    assert dbcan_tools_args("hmmer,diamond", legacy_script=False) == ["hmm,diamond"]


def test_legacy_dbcan_tools_arg_splits_comma_list():
    assert dbcan_tools_args("hmmer,diamond", legacy_script=True) == ["hmmer", "diamond"]


def test_legacy_dbcan_tools_arg_accepts_space_list():
    assert dbcan_tools_args("hmmer diamond", legacy_script=True) == ["hmmer", "diamond"]


def test_modern_dbcan_with_gff_uses_easy_cgc():
    command = build_modern_run_dbcan_command(
        run_dbcan_exec="run_dbcan",
        faa=Path("proteins.faa"),
        dbcan_out=Path("out/dbcan"),
        db_dir=Path("db"),
        gff=Path("genes.gff"),
        tools="hmmer,diamond",
        cpus=8,
    )

    assert command[:2] == ["run_dbcan", "easy_CGC"]
    assert "--input_gff" in command
    assert command[command.index("--methods") + 1] == "hmm,diamond"
    assert command[command.index("--threads") + 1] == "8"


def test_modern_dbcan_without_gff_uses_cazyme_annotation():
    command = build_modern_run_dbcan_command(
        run_dbcan_exec="run_dbcan",
        faa=Path("proteins.faa"),
        dbcan_out=Path("out/dbcan"),
        db_dir=Path("db"),
        gff=None,
        tools="hmmer,diamond",
        cpus=4,
    )

    assert command[:2] == ["run_dbcan", "CAZyme_annotation"]
    assert "--input_gff" not in command


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
