from pulsar_agarase.pipeline import dbcan_tools_args


def test_modern_dbcan_tools_arg_keeps_comma_list():
    assert dbcan_tools_args("hmmer,diamond", legacy_script=False) == ["hmmer,diamond"]


def test_legacy_dbcan_tools_arg_splits_comma_list():
    assert dbcan_tools_args("hmmer,diamond", legacy_script=True) == ["hmmer", "diamond"]


def test_legacy_dbcan_tools_arg_accepts_space_list():
    assert dbcan_tools_args("hmmer diamond", legacy_script=True) == ["hmmer", "diamond"]
