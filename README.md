# agarase-score

`agarase-score` is a small command-line tool for architecture-based scoring of
agarolytic PULs and candidate GH family additions.

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

## Installation

Install the external annotation tools first. With conda/mamba:

```bash
git clone https://github.com/YOUR_USER/agarase-score.git
cd agarase-score
mamba env create -f environment.yml
mamba activate agarase-score
```

If you already have Prodigal and dbCAN installed, a lightweight editable install
is enough:

```bash
git clone https://github.com/YOUR_USER/agarase-score.git
cd agarase-score
python3 -m pip install -e .
```

Check commands:

```bash
prodigal -v
run_dbcan --help
agarase-score --help
```

Prepare the dbCAN database once, or let `run-genome` do this automatically
when the database directory is empty:

```bash
agarase-score setup-dbcan \
  --db-dir dbcan_db \
  --min-free-gb 20
```

## Usage

### 1. Start from a genome FASTA

If the input is nucleotide FASTA (`.fna`, `.fa`, `.fasta`), `agarase-score`
runs Prodigal first, then runs dbCAN/CGCFinder, extracts features, and scores
the genome:

```bash
agarase-score run-genome \
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
agarase-score run-genome \
  --genome proteins.faa \
  --input-type faa \
  --gff genes.gff \
  --out-dir output/GenomeA \
  --dbcan-db /path/to/dbcan_db
```

External programs required for this command:

- `prodigal` for nucleotide genome input
- `run_dbcan`
- a prepared dbCAN database from `agarase-score setup-dbcan`

### 2. Score an existing feature table

```bash
agarase-score score-table \
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
agarase-score features-from-dbcan \
  --dbcan-dir dbcan_outputs \
  --output output/features.tsv

agarase-score score-table \
  --input output/features.tsv \
  --output output/predictions.tsv
```

Optional metadata:

```bash
agarase-score features-from-dbcan \
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
