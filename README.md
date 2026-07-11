# PULSAR

PULSAR (PUL-based Selection of AgaRase) is a small command-line tool for
architecture-based scoring of agarolytic PULs and candidate GH family additions
that can enhance fitness.

The model is rule-based and interpretable. PULSAR detects agar-PUL candidates
from the genome-wide distribution of agar-related marker genes, then scores
candidate GH family additions from the detected locus architecture.

## Model Concept

The score is based on two partially overlapping agarolytic pathways:

- Core pathway: `polysaccharide -> GH16/GH86/GH118 -> NAOS -> GH50 -> NA2 -> GH117 -> monosaccharide`
- Auxiliary pathway: `polysaccharide -> GH96/GH2 -> AOS -> GH117`

Key interpretation rules:

- Agar-PUL detection is not gated by CGCFinder. CGCFinder output is kept as a
  strict comparison layer, but the default detector uses a multiscale scan over
  marker-gene positions.
- The default scan uses gene-count windows. Base-pair windows are available with
  `--scan-unit bp`.
- Scan detection counts agarase and L-AHG metabolism marker genes together.
  Pathway-combination logic is applied only after candidate PULs are detected.
- GH2 is not a stand-alone agar-PUL detection marker. If GH2 falls inside a
  detected agar-PUL, it is reported as an auxiliary/supporting feature.
- If agarase-family genes are detected genome-wide but no statistically
  supported marker-gene cluster is detected, PULSAR reports
  `genome_wide_agarase_without_locus_context` and does not recommend a GH
  addition from PUL architecture alone.

## How PULSAR Works

PULSAR runs the analysis in four layers:

1. Annotation layer
   - For nucleotide FASTA input, Prodigal predicts proteins and GFF coordinates.
   - dbCAN/CGCFinder annotates CAZyme, transporter, TF, STP, and CGC context.

2. Feature layer
   - dbCAN outputs define genome-wide marker genes.
   - PULSAR builds a marker table with gene order/coordinates where available.
   - Detection markers include agarolytic GH16 subtypes, GH50, GH86, GH96,
     GH117, GH118, and L-AHG metabolism genes detected from GFF annotations.
   - GH2 is retained only as an auxiliary/supporting marker inside detected
     PULs.

3. Context layer
   - PULSAR scans multiple gene-count window sizes across each contig.
   - Windows are scored by marker-gene enrichment relative to genome-wide marker
     density.
   - Monte Carlo simulations randomly redistribute the same number of detection
     marker genes and estimate empirical p-values for the top non-overlapping
     intervals.
   - Significant intervals are reported as scan-detected agar-PULs.
   - Strict CGC/PUL context is still reported separately for comparison.
   - Genome-wide hits without locus context are reported but not used to make a
     PUL-based GH recommendation.

4. Scoring layer
   For each family g ∈ F = {GH2,GH16,GH50,GH86,GH96,GH117,GH118}, PULSAR calculates three types of counts:
   - The number of g genes located within strict PULs detected by CGCFinder.
   - The number of g genes located within broad loci reconstructed based on gene positions when no strict PUL is detected.
   - The total number of g genes detected across the entire genome.

## Installation

Install the external annotation tools first. With mamba (or conda):

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

With the current `run_dbcan` CLI, PULSAR calls `run_dbcan easy_CGC` for
genome/GFF-aware runs. Older single-script installations can still be used with
`--run-dbcan-script`.

Main outputs:

```text
output/GenomeA/
  work/prodigal.faa
  work/prodigal.gff
  dbcan/
  features.tsv
  marker_genes.tsv
  scan_candidate_windows.tsv
  scan_agar_puls.tsv
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
  --marker-output output/GenomeA/marker_genes.tsv \
  --candidate-windows-output output/GenomeA/scan_candidate_windows.tsv \
  --pul-output output/GenomeA/scan_agar_puls.tsv \
  --output output/GenomeA/predictions.tsv \
  --print-summary
```

Useful scan options:

```bash
pulsar score-dbcan \
  --dbcan-dir output/GenomeA/dbcan \
  --gff output/GenomeA/work/prodigal.gff \
  --scan-unit gene \
  --scan-windows 5,10,15,20,30,50 \
  --scan-permutations 999 \
  --output output/GenomeA/predictions.tsv
```

For base-pair windows:

```bash
pulsar score-dbcan \
  --dbcan-dir output/GenomeA/dbcan \
  --gff output/GenomeA/work/prodigal.gff \
  --scan-unit bp \
  --scan-windows 5000,10000,20000,50000 \
  --output output/GenomeA/predictions.tsv
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
