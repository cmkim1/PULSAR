# PULSAR

PULSAR (PUL-based Selection of AgaRase) is a small command-line tool for
architecture-based scoring of agarolytic PULs and candidate GH family additions.

The model is intentionally rule-based and interpretable. It uses the local PUL
architecture around agarolytic pathway genes rather than experimental outcomes.
In particular, GH16, GH86, and GH118 are treated as neutral core-opener
candidates unless the genome architecture itself supports one family over
another.

## Model Concept

The score is based on two partially overlapping agarolytic pathways:

- Core pathway: `polysaccharide -> GH16/GH86/GH118 -> NAOS -> GH50 -> NA2 -> GH117 -> monosaccharide`
- Auxiliary pathway: `polysaccharide -> GH96/GH2 -> AOS -> GH117`

Key interpretation rules:

- A strict agar-PUL containing GH117 is treated as a central PUL context.
- GH2 is kept as a local strict-PUL auxiliary context, not as a stand-alone
  agar-PUL marker.
- Non-strict GH2 signals are not used to drive recommendation scores.
- GH16, GH86, and GH118 are not ranked by fixed family preference.
- If GH16/GH86/GH118 are all absent from a GH50/GH117 core context, the model
  reports them as an unresolved equivalent core-opener group.
- If one of GH16/GH86/GH118 is present outside strict PULs but absent from the
  GH117-centered strict PUL, that family is reported as a split-locus candidate.
- If agarase-family genes are detected genome-wide but no strict CGC/PUL or
  broad colocalized locus is detected, PULSAR reports
  `genome_wide_agarase_without_locus_context` and does not recommend a GH
  addition from PUL architecture alone.

## How PULSAR Works

PULSAR runs the analysis in four layers:

1. Annotation layer
   - For nucleotide FASTA input, Prodigal predicts proteins and GFF coordinates.
   - dbCAN/CGCFinder annotates CAZyme, transporter, TF, STP, and CGC context.

2. Feature layer
   - `cgc.out` defines strict CGC/PUL-localized agarase-family counts.
   - `hmmer.out` and `diamond.out` define genome-wide agarase-family counts.
   - If `cgc.out` is empty but `cgc.gff` and DIAMOND/HMMER hits are available,
     PULSAR reconstructs a broad colocalized locus from gene order.

3. Context layer
   - Strict CGC/PUL context is preferred.
   - Broad colocalized locus context is used only when strict CGC/PUL context is
     absent.
   - Genome-wide hits without locus context are reported but not used to make a
     PUL-based GH recommendation.

4. Scoring layer
   - GH117-positive context is treated as the central agar-PUL context.
   - GH16/GH86/GH118 are neutral core-openers unless architecture points to a
     specific split-locus or missing-context pattern.
   - GH2 is considered only as local auxiliary context, not as a stand-alone
     agar-PUL marker.

## Installation

Install the external annotation tools first. With conda/mamba:

```bash
git clone https://github.com/YOUR_USER/PULSAR.git
cd PULSAR
mamba env create -f environment.yml
mamba activate pulsar
```

If you already have Prodigal and dbCAN installed, a lightweight editable install
is enough:

```bash
git clone https://github.com/YOUR_USER/PULSAR.git
cd PULSAR
python3 -m pip install -e .
```

Check commands:

```bash
prodigal -v
run_dbcan --help
pulsar --help
```

Or run PULSAR's environment check:

```bash
pulsar doctor --dbcan-db dbcan_db
```

Prepare the dbCAN database once, or let `run-genome` do this automatically
when the database directory is empty:

```bash
pulsar setup-dbcan \
  --db-dir dbcan_db \
  --min-free-gb 20
```

## Usage

### 1. Start from a genome FASTA

If the input is nucleotide FASTA (`.fna`, `.fa`, `.fasta`), PULSAR
runs Prodigal first, then runs dbCAN/CGCFinder, extracts features, and scores
the genome:

```bash
pulsar run-genome \
  --genome genome.fna \
  --out-dir output/GenomeA \
  --dbcan-db dbcan_db \
  --genome-id GenomeA \
  --taxname "Example species A" \
  --cpus 8
```

If `dbcan_db` is empty or missing, `run-genome` automatically runs:

```bash
run_dbcan database --db_dir dbcan_db
```

Use `--skip-dbcan-setup` to disable automatic database setup.

Main outputs:

```text
output/GenomeA/
  work/prodigal.faa
  work/prodigal.gff
  dbcan/
  features.tsv
  predictions.tsv
  logs/
```

If the input is already a protein FASTA (`.faa`), provide a matching GFF if
CGCFinder clustering is desired:

```bash
pulsar run-genome \
  --genome proteins.faa \
  --input-type faa \
  --gff genes.gff \
  --out-dir output/GenomeA \
  --dbcan-db /path/to/dbcan_db
```

Legacy `run_dbcan.py` installations can be used directly:

```bash
pulsar run-genome \
  --genome genome.fna \
  --out-dir output/GenomeA \
  --dbcan-db /path/to/dbcan_db \
  --run-dbcan-script /path/to/run_dbcan.py \
  --dbcan-file dbCAN-HMMdb-V9.txt \
  --skip-dbcan-setup
```

Check a legacy installation before running:

```bash
pulsar doctor \
  --run-dbcan-script /path/to/run_dbcan.py \
  --dbcan-db /path/to/dbcan_db
```

External programs required for this command:

- `prodigal` for nucleotide genome input
- `run_dbcan`
- a prepared dbCAN database from `pulsar setup-dbcan`

### 2. Score an existing feature table

```bash
pulsar score-table \
  --input examples/example_features.tsv \
  --output output/example_predictions.tsv \
  --summary output/example_summary.tsv
```

### 3. Extract features from dbCAN/CGCFinder output

The input directory should contain one dbCAN output directory per genome:

```text
dbcan_outputs/
  GenomeA/
    cgc.out
    hmmer.out
  GenomeB/
    cgc.out
    hmmer.out
```

Then run:

```bash
pulsar features-from-dbcan \
  --dbcan-dir dbcan_outputs \
  --output output/features.tsv

pulsar score-table \
  --input output/features.tsv \
  --output output/predictions.tsv
```

Optional metadata:

```bash
pulsar features-from-dbcan \
  --dbcan-dir dbcan_outputs \
  --metadata metadata.tsv \
  --output output/features.tsv
```

`metadata.tsv` must contain:

```text
genome  taxname
GenomeA Example_species_A
GenomeB Example_species_B
```

For one existing dbCAN output directory:

```bash
pulsar score-dbcan \
  --dbcan-dir output/GenomeA/dbcan \
  --genome-id GenomeA \
  --taxname "Example species A" \
  --features-output output/GenomeA/features.tsv \
  --output output/GenomeA/predictions.tsv \
  --print-summary
```

## New Server Quickstart

```bash
git clone https://github.com/cmkim1/PULSAR.git
cd PULSAR
mamba env create -f environment.yml
mamba activate pulsar

pulsar doctor
pulsar setup-dbcan --db-dir dbcan_db --min-free-gb 20
pulsar doctor --dbcan-db dbcan_db

pulsar run-genome \
  --genome /path/to/genome.fna \
  --out-dir output/GenomeA \
  --dbcan-db dbcan_db \
  --genome-id GenomeA \
  --taxname "Example species A" \
  --cpus 8
```

If the repository is private, clone with SSH or a GitHub token and make sure the
server account has access.

## Required Feature Table Columns

For each genome, the scoring command expects:

- `genome`
- `taxname`
- `strict_n_agar_loci`
- `broad_n_agar_loci`
- `has_genome_wide_annotation`
- For each family in `GH2, GH16, GH50, GH86, GH96, GH117, GH118`:
  - `strict_n_<FAMILY>`
  - `genome_n_<FAMILY>`
  - `broad_locus_n_<FAMILY>`
  - `outside_strict_n_<FAMILY>`

The output appends:

- `<FAMILY>_score`
- `recommended_GH_group`
- `top_recommended_GH`
- `prediction_class`
- `model_confidence`
- `model_rationale`
- `central_pul_status`
- `core_pathway_status`
- `auxiliary_pathway_status`
- `core_opener_status`

## Important Notes

This is a hypothesis-generation model, not a validated classifier of growth
phenotype. The score should be used to prioritize GH families or PUL
architectures for follow-up analysis.

For publication-scale use, report the exact version, input feature table, and
the rule descriptions in this README.
