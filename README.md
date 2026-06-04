# GLADE

GLADE: Accurate inference of Gains, Losses, Ancestral genomes, and Duplication Events for comparative genomics

GLADE is a Python tool for reconstructing the full evolutionary history of orthogroups — including gene gains, losses, duplications, and ancestral gene repertoires — using only an OrthoFinder v3 results directory as input. GLADE maps every event onto the species tree and produces rich output for comparative genomics.

<p align="center">
  <img src="https://github.com/lauriebelch/GLADE/blob/main/GLADE.png" width="50%">
</p>


## Table of contents
- [What is GLADE?](#What-is-GLADE)
- [Installation](#Installation)
- [How-to-use](#Simple-usage)
- [Output files](#Output-files)
- [Example-data](#Example-data)
- [Citation](#Citation)
- [Compatibility notes](#Compatibility-notes)

## What is GLADE?

GLADE reconstructs the history of orthogroups — defined as sets of genes descended from a single gene in the most recent common ancestor — across a species tree.

Given a complete OrthoFinder v3 run, GLADE:
- Identifies where each orthogroup first appeared (gain)
- Detects losses
- Identifies gene duplication events
- Reconstructs ancestral gene content at every internal node
- Quantifies orthogroup size changes along every branch
- Outputs complete evolutionary histories for all orthogroups

<p align="center">
<img src="glade_workflow_.png" alt="workflow" width="700"/>
</p>

## Installation

GLADE requires Python 3.9 or later and the same dependencies as OrthoFinder. We recommend running GLADE in an OrthoFinder conda environment.

See the OrthoFinder GitHub for installation details: https://github.com/OrthoFinder/OrthoFinder?tab=readme-ov-file#installation

## Simple usage

```
python GLADE.py -f path/to/orthofinder/results -t threads [default=8]
```

**All options:**

| Flag | Description | Default |
|---|---|---|
| `-f` / `--folder` | Path to OrthoFinder results directory | required |
| `-t` / `--threads` | Number of threads | 8 |
| `-m` / `--min-genes` | Minimum genes per OG for gain/loss analysis | 4 |
| `-o` / `--output` | Output folder for final results | same as `-f` |
| `-s` / `--species-tree` | Custom species tree in Newick format | OrthoFinder tree |
| `-g` / `--gene-trees-dir` | Directory of per-OG IQ-TREE gene trees | OrthoFinder trees |

### Custom species tree (`-s`)

**We strongly recommend providing a custom species tree via `-s`.** The OrthoFinder species tree is inferred from gene tree topologies using the STRIDE algorithm, which can place the root incorrectly when gene trees are noisy. A dedicated phylogenetic analysis (e.g. IQ-TREE on a supermatrix of single-copy orthologs) produces more reliable branch lengths and topology. The root position directly determines the direction of all inferred gains and losses.

The custom tree must be:
- **Rooted** — an unrooted tree (e.g. a raw IQ-TREE `.treefile` with a trifurcating root) will cause errors
- In **Newick format**
- Using **leaf names that match OrthoFinder species names** (FASTA filenames without extension, e.g. `Schizosaccharomyces_pombe`)

Internal node labels are added automatically; you do not need to provide them.

### External gene trees (`-g`)

By default GLADE uses the gene trees produced by OrthoFinder (distance-based). You can substitute higher-quality ML gene trees built with IQ-TREE or another tool via `-g`.

```bash
python GLADE.py -f Results/ -g iqtree_trees/ -s species.treefile -t 16 -o output/
```

**Preparing gene trees for `-g`:**

1. Extract per-OG sequences from OrthoFinder's `Orthogroups/Orthogroup_Sequences/`
2. Align each OG (e.g. with MAFFT or MUSCLE)
3. Run IQ-TREE on each alignment; no bootstrap needed:
   ```bash
   iqtree -s OG0000000.aln -m MFP -T AUTO --prefix OG0000000
   ```
4. Collect all `.treefile` outputs in one directory and pass it to `-g`

**Requirements:**
- File naming: `{OG}.treefile` (e.g. `OG0000000.treefile`)
- Leaf names: `{Species}_{geneID}` (e.g. `Debaryomyces_hansenii_DEHA2A00748g`)
- OG names must match `Orthogroups.tsv`
- Trees do not need to be pre-rooted — GLADE roots each tree automatically using the OrthoFinder S_IO/S_AD algorithm
- OGs without a `.treefile` are skipped; GLADE prints a coverage summary at startup

### Other options

Use `-m 1` to include all OGs (including single-copy and small families) in gain/loss analysis.

Use `-o` to write `GainsLossDuplication/` and `AncestralGenomes/` to a custom path. Intermediate files in `WorkingDirectory/GladeWD/` are unaffected.

If you are running GLADE on an OrthoFinder assign run, copy the proteomes from the core run into the assign results `WorkingDirectory/`:

```bash
cp core/WorkingDirectory/*.fa assign/WorkingDirectory/
```

## Output files

GLADE produces a structured directory containing:

**GainsLossDuplication/**
- `Gains.tsv` — where each orthogroup first appeared
- `Loss_speciation.tsv` — orthogroup losses due to speciation
- `Loss_postduplication.tsv` — losses after duplication events
- `Duplications.tsv` — duplication events with support values
- `Branch_statistics.tsv` — event counts per species-tree branch
- `*_bybranch.tsv` — expanded lists of orthogroups per event type
- `OrthogroupBranchChange.tsv` — size changes per branch

**AncestralGenomes/**
- One FASTA per internal node containing reconstructed ancestral sequences
- `AncestralGenomes.txt` — summary statistics
- `Ancestral_HOG_counts.csv` — orthogroup copy numbers for all nodes

## Example data

Unzip `ExampleData.zip`, which contains an OrthoFinder results directory on a small dataset. Then run:

```bash
python GLADE.py -f ExampleData/OrthoFinder/Results_ExampleDataGLADE/
```

## Citation

Belcher L.J. & Kelly S. (2026) GLADE: Accurate inference of Gains, Losses, Ancestral genomes, and Duplication Events for comparative genomics. [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.27.702036v1)

---

## Compatibility notes

The following issues were identified when running GLADE with OrthoFinder v3 output and have been patched.

### OrthoFinder output directory structure

Newer versions of OrthoFinder place `Resolved_Gene_Trees/` inside `WorkingDirectory/` rather than in the results root. Create the required input file manually before running GLADE:

```bash
mkdir -p path/to/Results/Resolved_Gene_Trees
bash -c 'for f in path/to/Results/WorkingDirectory/Resolved_Gene_Trees/OG*.txt; do
  og=$(basename "$f" .txt)
  printf "%s: %s\n" "$og" "$(cat "$f")"
done > path/to/Results/Resolved_Gene_Trees/Resolved_Gene_Trees.txt'
```

Use `bash -c '...'` (non-interactive) to avoid terminal escape sequences being written into the file.

### ete3 visualization (PyQt5)

`ete3.TreeStyle` requires Qt. If running in a conda environment without Qt:

```bash
pip install PyQt5
```

### Bug fixes (v3.1.3 compatibility)

Several issues were identified and fixed when running GLADE against OrthoFinder v3.1.3:

- **KeyError in `GainAndLossAndDuplication.py`**: OGs with ≥4 genes but no resolved gene tree (e.g. all genes from a single species) caused a `KeyError` in `FindDuplicationsParallel`. Fixed by skipping such OGs — they cannot have duplication events.
- **ValueError in `AncestralGenome.py`**: The same class of OGs caused a `ValueError` in `ProcessOrthogroupCurrent`. Fixed by checking for the OG's presence in `Resolved_Gene_Trees.txt` before calling `GetAncestralGenes`.
- **TreeError in `AncestralGenome.py`**: `GetAncestralGenes` could pass an empty leaf list to `gene_tree.prune()` when all remaining leaves were flagged as post-duplication duplicates. Fixed by returning `None` when `keep_leaves` is empty.
- **OG–gene mismatch in `Duplications.tsv`**: OrthoFinder v3 uses different OG naming between `Orthogroups.tsv` and `Resolved_Gene_Trees.txt` — a single gene-family tree spans genes from multiple fine-grained orthogroups. GLADE previously assumed a 1:1 correspondence, producing duplication records where the OG label and gene names referred to completely different groups. Fixed by building a gene→gene-tree reverse index at startup and pruning each tree to its OG's actual members before duplication analysis.
