"""Command-line interface for agarase-score."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .features import features_from_dbcan_root
from .model import score_dataframe
from .pipeline import score_genome, setup_dbcan_database


def score_table(args: argparse.Namespace) -> None:
    features = pd.read_csv(args.input, sep="\t")
    scored = score_dataframe(features)
    scored.to_csv(args.output, sep="\t", index=False)
    if args.summary:
        summary = (
            scored.groupby("recommended_GH_group")
            .agg(n=("genome", "count"))
            .reset_index()
            .sort_values(["n", "recommended_GH_group"], ascending=[False, True])
        )
        summary.to_csv(args.summary, sep="\t", index=False)


def features_from_dbcan(args: argparse.Namespace) -> None:
    features = features_from_dbcan_root(Path(args.dbcan_dir), Path(args.metadata) if args.metadata else None)
    features.to_csv(args.output, sep="\t", index=False)


def run_genome(args: argparse.Namespace) -> None:
    score_genome(
        genome=Path(args.genome),
        out_dir=Path(args.out_dir),
        dbcan_db=Path(args.dbcan_db),
        input_type=args.input_type,
        gff=Path(args.gff) if args.gff else None,
        genome_id=args.genome_id,
        taxname=args.taxname,
        prodigal_bin=args.prodigal_bin,
        run_dbcan_bin=args.run_dbcan_bin,
        tools=args.tools,
        dbcan_file=args.dbcan_file,
        cpus=args.cpus,
        auto_setup_dbcan=not args.skip_dbcan_setup,
        min_free_gb=args.min_free_gb,
    )


def setup_dbcan(args: argparse.Namespace) -> None:
    db_dir = setup_dbcan_database(
        db_dir=Path(args.db_dir),
        run_dbcan_bin=args.run_dbcan_bin,
        min_free_gb=args.min_free_gb,
        force=args.force,
    )
    print(f"dbCAN database is ready: {db_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="agarase-score",
        description="Score agarolytic PUL architecture and recommend candidate GH families.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_score = sub.add_parser("score-table", help="Score a genome-level feature table.")
    p_score.add_argument("-i", "--input", required=True, help="Input feature table TSV.")
    p_score.add_argument("-o", "--output", required=True, help="Output scored prediction TSV.")
    p_score.add_argument("--summary", help="Optional recommendation summary TSV.")
    p_score.set_defaults(func=score_table)

    p_features = sub.add_parser("features-from-dbcan", help="Extract features from per-genome dbCAN output directories.")
    p_features.add_argument("--dbcan-dir", required=True, help="Directory containing one dbCAN output subdirectory per genome.")
    p_features.add_argument("-o", "--output", required=True, help="Output feature table TSV.")
    p_features.add_argument("--metadata", help="Optional TSV with genome and taxname columns.")
    p_features.set_defaults(func=features_from_dbcan)

    p_setup = sub.add_parser("setup-dbcan", help="Download/prepare the dbCAN database using run_dbcan.")
    p_setup.add_argument("--db-dir", required=True, help="Directory where dbCAN database files will be stored.")
    p_setup.add_argument("--run-dbcan-bin", default="run_dbcan", help="run_dbcan executable name or path.")
    p_setup.add_argument("--min-free-gb", type=float, default=20.0, help="Minimum free disk space required before download. Default: 20.")
    p_setup.add_argument("--force", action="store_true", help="Run database setup even if a dbCAN-HMMdb file already exists.")
    p_setup.set_defaults(func=setup_dbcan)

    p_run = sub.add_parser("run-genome", help="Run Prodigal/dbCAN/CGCFinder and score one genome.")
    p_run.add_argument("--genome", required=True, help="Input genome .fna/.fa/.fasta or protein .faa file.")
    p_run.add_argument("--out-dir", required=True, help="Output directory for work files, dbCAN output, and predictions.")
    p_run.add_argument("--dbcan-db", required=True, help="dbCAN database directory used by run_dbcan.")
    p_run.add_argument("--input-type", choices=["auto", "fna", "faa"], default="auto", help="Input type. Default: auto.")
    p_run.add_argument("--gff", help="GFF file for --input-type faa if CGCFinder clustering is desired.")
    p_run.add_argument("--genome-id", help="Genome ID to write in the output table.")
    p_run.add_argument("--taxname", help="Taxon/strain name to write in the output table.")
    p_run.add_argument("--prodigal-bin", default="prodigal", help="Prodigal executable name or path.")
    p_run.add_argument("--run-dbcan-bin", default="run_dbcan", help="run_dbcan executable name or path.")
    p_run.add_argument("--tools", default="hmmer,diamond", help="dbCAN tools argument. Default: hmmer,diamond.")
    p_run.add_argument("--dbcan-file", help="Optional dbCAN HMM database filename, e.g. dbCAN-HMMdb-V9.txt.")
    p_run.add_argument("--cpus", type=int, default=4, help="CPU threads for HMMER/DIAMOND. Default: 4.")
    p_run.add_argument("--min-free-gb", type=float, default=20.0, help="Minimum free disk space before automatic dbCAN setup. Default: 20.")
    p_run.add_argument("--skip-dbcan-setup", action="store_true", help="Do not automatically run run_dbcan database if --dbcan-db is missing.")
    p_run.set_defaults(func=run_genome)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
