"""Pipeline wrappers for Prodigal and dbCAN/CGCFinder."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from .features import features_from_dbcan_dir
from .model import score_dataframe


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")
    return path


def run_command(command: list[str], log_path: Path, allow_outputs: list[Path] | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write("$ " + " ".join(command) + "\n\n")
        proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        if allow_outputs and any(path.exists() and path.stat().st_size > 0 for path in allow_outputs):
            with log_path.open("a") as log:
                log.write(
                    f"\nPULSAR warning: command exited with {proc.returncode}, "
                    "but expected dbCAN output files were created; continuing.\n"
                )
            return
        raise RuntimeError(f"Command failed with exit code {proc.returncode}. See log: {log_path}")


def free_gb(path: Path) -> float:
    path.mkdir(parents=True, exist_ok=True)
    usage = shutil.disk_usage(path)
    return usage.free / (1024**3)


def setup_dbcan_database(
    db_dir: Path,
    run_dbcan_bin: str = "run_dbcan",
    min_free_gb: float = 20.0,
    force: bool = False,
) -> Path:
    """Download/prepare the dbCAN database through run_dbcan."""

    run_dbcan_exec = require_executable(run_dbcan_bin)
    db_dir.mkdir(parents=True, exist_ok=True)
    available = free_gb(db_dir)
    if available < min_free_gb:
        raise RuntimeError(
            f"Not enough free disk space for dbCAN database setup: {available:.1f} GB available, "
            f"{min_free_gb:.1f} GB required. Use --min-free-gb to change this threshold."
        )

    marker_files = ["dbCAN-HMMdb-V12.txt", "dbCAN-HMMdb-V11.txt", "dbCAN-HMMdb-V10.txt", "dbCAN-HMMdb-V9.txt"]
    if not force and any((db_dir / marker).exists() for marker in marker_files):
        return db_dir

    log = db_dir / "setup_dbcan.log"
    command = [run_dbcan_exec, "database", "--db_dir", str(db_dir)]
    run_command(command, log)
    return db_dir


def dbcan_database_present(db_dir: Path) -> bool:
    marker_files = ["dbCAN-HMMdb-V12.txt", "dbCAN-HMMdb-V11.txt", "dbCAN-HMMdb-V10.txt", "dbCAN-HMMdb-V9.txt"]
    return any((db_dir / marker).exists() for marker in marker_files)


def detect_input_type(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith((".fna", ".fa", ".fasta", ".fasta.gz", ".fa.gz", ".fna.gz")):
        return "fna"
    if suffixes.endswith((".faa", ".faa.gz")):
        return "faa"
    raise ValueError("Could not infer input type from extension. Use --input-type fna or --input-type faa.")


def run_prodigal(fna: Path, work_dir: Path, prodigal_bin: str = "prodigal") -> tuple[Path, Path]:
    prodigal = require_executable(prodigal_bin)
    faa = work_dir / "prodigal.faa"
    gff = work_dir / "prodigal.gff"
    log = work_dir / "logs" / "prodigal.log"
    command = [
        prodigal,
        "-i",
        str(fna),
        "-a",
        str(faa),
        "-o",
        str(gff),
        "-f",
        "gff",
        "-p",
        "single",
    ]
    run_command(command, log)
    return faa, gff


def dbcan_tools_args(tools: str, legacy_script: bool = False) -> list[str]:
    if not legacy_script:
        return [tools]
    return [tool.strip() for tool in tools.replace(",", " ").split() if tool.strip()]


def run_dbcan(
    faa: Path,
    out_dir: Path,
    db_dir: Path,
    gff: Path | None,
    run_dbcan_bin: str = "run_dbcan",
    run_dbcan_script: Path | None = None,
    tools: str = "hmmer,diamond",
    dbcan_file: str | None = None,
    cpus: int = 4,
) -> Path:
    dbcan_out = out_dir / "dbcan"
    log = out_dir / "logs" / "run_dbcan.log"
    if run_dbcan_script is not None:
        command = [sys.executable, str(run_dbcan_script)]
    else:
        run_dbcan_exec = require_executable(run_dbcan_bin)
        command = [run_dbcan_exec]
    command.extend(
        [
            str(faa),
            "protein",
            "--out_dir",
            str(dbcan_out),
            "--db_dir",
            str(db_dir),
            "--tools",
        ]
        + dbcan_tools_args(tools, legacy_script=run_dbcan_script is not None)
        + [
            "--hmm_cpu",
            str(cpus),
            "--dia_cpu",
            str(cpus),
        ]
    )
    if dbcan_file:
        command.extend(["--dbCANFile", dbcan_file])
    if gff is not None:
        command.extend(["--cluster", str(gff)])
    run_command(command, log, allow_outputs=[dbcan_out / "cgc.out", dbcan_out / "hmmer.out", dbcan_out / "diamond.out"])
    return dbcan_out


def score_genome(
    genome: Path,
    out_dir: Path,
    dbcan_db: Path,
    input_type: str = "auto",
    gff: Path | None = None,
    genome_id: str | None = None,
    taxname: str | None = None,
    prodigal_bin: str = "prodigal",
    run_dbcan_bin: str = "run_dbcan",
    run_dbcan_script: Path | None = None,
    tools: str = "hmmer,diamond",
    dbcan_file: str | None = None,
    cpus: int = 4,
    auto_setup_dbcan: bool = True,
    min_free_gb: float = 20.0,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = out_dir / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    kind = detect_input_type(genome, input_type)
    if kind == "fna":
        faa, called_gff = run_prodigal(genome, work_dir, prodigal_bin=prodigal_bin)
        cluster_gff = called_gff
    elif kind == "faa":
        faa = genome
        cluster_gff = gff
    else:
        raise ValueError(f"Unsupported input type: {kind}")

    if auto_setup_dbcan and not dbcan_database_present(dbcan_db):
        setup_dbcan_database(
            db_dir=dbcan_db,
            run_dbcan_bin=run_dbcan_bin,
            min_free_gb=min_free_gb,
            force=False,
        )

    dbcan_out = run_dbcan(
        faa=faa,
        out_dir=out_dir,
        db_dir=dbcan_db,
        gff=cluster_gff,
        run_dbcan_bin=run_dbcan_bin,
        run_dbcan_script=run_dbcan_script,
        tools=tools,
        dbcan_file=dbcan_file,
        cpus=cpus,
    )
    row = features_from_dbcan_dir(dbcan_out, genome=genome_id or genome.stem, taxname=taxname or genome_id or genome.stem)
    features = pd.DataFrame([row])
    scored = score_dataframe(features)
    features.to_csv(out_dir / "features.tsv", sep="\t", index=False)
    scored.to_csv(out_dir / "predictions.tsv", sep="\t", index=False)
    return scored
