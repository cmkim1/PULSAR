"""Pipeline wrappers for Prodigal and dbCAN/CGCFinder."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from .features import dbcan_feature_tables
from .model import score_dataframe


LEGACY_DBCAN_MARKERS = ["dbCAN-HMMdb-V12.txt", "dbCAN-HMMdb-V11.txt", "dbCAN-HMMdb-V10.txt", "dbCAN-HMMdb-V9.txt"]
MODERN_DBCAN_MARKERS = ["CAZy.dmnd", "dbCAN.hmm", "dbCAN-sub.hmm", "TCDB.dmnd", "TF.hmm", "STP.hmm"]


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")
    return path


def run_command(
    command: list[str],
    log_path: Path,
    allow_outputs: list[Path] | None = None,
    cwd: Path | None = None,
    label: str | None = None,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if label:
        print(f"PULSAR: starting {label}; log={log_path}", flush=True)
    with log_path.open("w") as log:
        log.write("$ " + " ".join(command) + "\n\n")
        if cwd is not None:
            log.write(f"# cwd: {cwd}\n\n")
        proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, cwd=cwd)
    if proc.returncode != 0:
        if allow_outputs and any(path.exists() and path.stat().st_size > 0 for path in allow_outputs):
            with log_path.open("a") as log:
                log.write(
                    f"\nPULSAR warning: command exited with {proc.returncode}, "
                    "but expected dbCAN output files were created; continuing.\n"
                )
            return
        raise RuntimeError(f"Command failed with exit code {proc.returncode}. See log: {log_path}")
    if label:
        print(f"PULSAR: finished {label}", flush=True)


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

    if not force and dbcan_database_present(db_dir):
        return db_dir

    log = db_dir / "setup_dbcan.log"
    command = [run_dbcan_exec, "database", "--db_dir", str(db_dir)]
    run_command(command, log)
    return db_dir


def dbcan_database_present(db_dir: Path) -> bool:
    has_legacy = any((db_dir / marker).is_file() for marker in LEGACY_DBCAN_MARKERS)
    has_modern = all((db_dir / marker).is_file() for marker in MODERN_DBCAN_MARKERS)
    has_partial = any(db_dir.glob("*.part"))
    return (has_legacy or has_modern) and not has_partial


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
    run_command(command, log, label="Prodigal gene prediction")
    return faa, gff


def sanitize_gff_for_dbcan(gff: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gff.open() as source, destination.open("w") as sink:
        for line in source:
            if line.startswith("#") or not line.strip():
                sink.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                sink.write(line)
                continue
            attrs = [item.strip() for item in parts[8].split(";") if item.strip() and "=" in item]
            if not attrs:
                attrs = [f"ID={parts[0]}_{parts[3]}_{parts[4]}"]
            parts[8] = ";".join(attrs)
            sink.write("\t".join(parts) + "\n")
    return destination


def dbcan_tools_args(tools: str, legacy_script: bool = False) -> list[str]:
    parsed = [tool.strip() for tool in tools.replace(",", " ").split() if tool.strip()]
    if legacy_script:
        return parsed

    modern_names = {"hmmer": "hmm", "hmm": "hmm", "diamond": "diamond", "dbCANsub": "dbCANsub", "dbcansub": "dbCANsub"}
    modern = [modern_names.get(tool, tool) for tool in parsed]
    return [",".join(modern)]


def build_modern_run_dbcan_command(
    run_dbcan_exec: str,
    faa: Path,
    dbcan_out: Path,
    db_dir: Path,
    gff: Path | None,
    tools: str,
    cpus: int,
) -> list[str]:
    if gff is not None:
        return [
            run_dbcan_exec,
            "easy_CGC",
            "--mode",
            "protein",
            "--input_raw_data",
            str(faa),
            "--input_gff",
            str(gff),
            "--gff_type",
            "prodigal",
            "--output_dir",
            str(dbcan_out),
            "--db_dir",
            str(db_dir),
            "--methods",
            *dbcan_tools_args(tools, legacy_script=False),
            "--threads",
            str(cpus),
        ]

    return [
        run_dbcan_exec,
        "CAZyme_annotation",
        "--mode",
        "protein",
        "--input_raw_data",
        str(faa),
        "--output_dir",
        str(dbcan_out),
        "--db_dir",
        str(db_dir),
        "--methods",
        *dbcan_tools_args(tools, legacy_script=False),
        "--threads",
        str(cpus),
    ]


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
    faa = faa.resolve()
    out_dir = out_dir.resolve()
    db_dir = db_dir.resolve()
    gff = gff.resolve() if gff is not None else None

    dbcan_out = out_dir / "dbcan"
    log = out_dir / "logs" / "run_dbcan.log"
    if gff is not None:
        gff = sanitize_gff_for_dbcan(gff, out_dir / "work" / "dbcan_cluster.gff")

    command_cwd = None
    if run_dbcan_script is not None:
        run_dbcan_script = run_dbcan_script.resolve()
        command = [sys.executable, str(run_dbcan_script)]
        command_cwd = run_dbcan_script.parent
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
            + dbcan_tools_args(tools, legacy_script=True)
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
    else:
        run_dbcan_exec = require_executable(run_dbcan_bin)
        command = build_modern_run_dbcan_command(
            run_dbcan_exec=run_dbcan_exec,
            faa=faa,
            dbcan_out=dbcan_out,
            db_dir=db_dir,
            gff=gff,
            tools=tools,
            cpus=cpus,
        )
    allow_outputs = [dbcan_out / "cgc.out", dbcan_out / "cgc.gff"] if gff is not None else [dbcan_out / "hmmer.out", dbcan_out / "diamond.out"]
    run_command(command, log, allow_outputs=allow_outputs, cwd=command_cwd, label="dbCAN/CGCFinder")
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
    scan_windows: list[int] | None = None,
    scan_permutations: int = 999,
    scan_seed: int = 1,
    scan_alpha: float = 0.05,
    scan_unit: str = "gene",
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = out_dir / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    print(f"PULSAR: output directory={out_dir}", flush=True)

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
    row, marker_table, candidate_windows, scan_puls = dbcan_feature_tables(
        dbcan_out,
        genome=genome_id or genome.stem,
        taxname=taxname or genome_id or genome.stem,
        gff_path=cluster_gff,
        scan_windows=scan_windows,
        scan_permutations=scan_permutations,
        scan_seed=scan_seed,
        scan_alpha=scan_alpha,
        scan_unit=scan_unit,
    )
    features = pd.DataFrame([row])
    scored = score_dataframe(features)
    features.to_csv(out_dir / "features.tsv", sep="\t", index=False)
    marker_table.to_csv(out_dir / "marker_genes.tsv", sep="\t", index=False)
    candidate_windows.to_csv(out_dir / "scan_candidate_windows.tsv", sep="\t", index=False)
    scan_puls.to_csv(out_dir / "scan_agar_puls.tsv", sep="\t", index=False)
    scored.to_csv(out_dir / "predictions.tsv", sep="\t", index=False)
    print(f"PULSAR: wrote {out_dir / 'predictions.tsv'}", flush=True)
    return scored
