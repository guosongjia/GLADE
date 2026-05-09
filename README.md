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
- [Notes on compatibility and bug fixes](#Notes-on-compatibility-and-bug-fixes)

## What is GLADE?

GLADE reconstructs the history of orthogroups — defined as sets of genes descended from a single gene in the most recent common ancestor — across a species tree.

Given a complete OrthoFinder v3 run, GLADE:
- Identifies where each orthogroup first appeared (gain)
- Detects losses
- Identifies gene duplication events
- Reconstructs ancestral gene content at every internal node,
- Quantifies orthogroup size changes along every branch,
- Outputs complete evolutionary histories for all orthogroups.

<p align="center">
<img src="glade_workflow_.png" alt="workflow" width="700"/>
</p>

## Installation

GLADE requires Python 3.9 or later

GLADE requires the same dependencies as OrthoFinder. We reccommend that you run GLADE in an orthofinder conda environment.

See the OrthoFinder github for details on how to set this up https://github.com/OrthoFinder/OrthoFinder?tab=readme-ov-file#installation

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
| `-s` / `--species-tree` | Path to a custom species tree in Newick format (e.g. IQ-TREE treefile); overrides OrthoFinder species tree; internal node labels added automatically | OrthoFinder tree |

Use `-m 1` to include all OGs (including single-copy and small families) in gain/loss analysis. Duplication analysis is only performed for OGs with a resolved gene tree regardless of this setting.

Use `-o` to write `GainsLossDuplication/` and `AncestralGenomes/` to a custom path. Intermediate files in `WorkingDirectory/GladeWD/` are unaffected.

If you are running GLADE on an OrthoFinder assign run - you need to add the proteomes from the core run to the assign results directory

e.g. cd to core/WorkingDirectory and cp *.fa to the assign/WorkingDirectory


## Output files

GLADE produces a structured directory containing:

1. Gains, Losses, and Duplications (GainsLossDuplication/)
- Gains.tsv — where each orthogroup first appeared
- Loss_speciation.tsv — orthogroup losses due to speciation
- Loss_postduplication.tsv — losses after duplication events
- Duplications.tsv — duplication events with support values
- Branch_statistics.tsv — event counts per species-tree branch
- *_bybranch.tsv — expanded lists of orthogroups per event type
- extant_OG_counts.tsv — gene counts in extant species
- OrthogroupBranchChange.tsv — size changes per branch

2. AncestralGenomes/

One FASTA file per internal node containing reconstructed ancestral sequences
- AncestralGenomes.txt — summary statistics
- Ancestral_HOG_counts.csv — orthogroup copy numbers for all nodes

## Example data

Unzip the ExampleData.zip file, which contains an OrthoFinder results directory on a small dataset.
Then run:

```
python GLADE.py -f ExampleData/OrthoFinder/Results_ExampleDataGLADE/
```

## Citation

Belcher L.J. & Kelly S. (2026) GLADE: Accurate inference of Gains, Losses, Ancestral genomes, and Duplication Events for comparative genomics. [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.27.702036v1)

---

## Notes on compatibility and bug fixes

The following issues were identified when running GLADE with OrthoFinder v3 output and have been patched in this copy of the scripts.

### OrthoFinder output directory structure

Newer versions of OrthoFinder place `Resolved_Gene_Trees/` inside `WorkingDirectory/` rather than in the results root. GLADE expects it at the results root. Fix: create the required input file manually before running GLADE:

```bash
mkdir -p path/to/Results/Resolved_Gene_Trees
bash -c 'for f in path/to/Results/WorkingDirectory/Resolved_Gene_Trees/OG*.txt; do og=$(basename "$f" .txt); printf "%s: %s\n" "$og" "$(cat "$f")"; done > path/to/Results/Resolved_Gene_Trees/Resolved_Gene_Trees.txt'
```

Note: use `bash -c '...'` (non-interactive) to avoid terminal escape sequences being written into the file.

### ete3 visualization (PyQt5)

`ete3.TreeStyle` and related classes require Qt. If running in a conda environment without Qt, install PyQt5 via pip:

```bash
pip install PyQt5
```

### Bug fix: OGs without resolved gene trees cause KeyError (`GainAndLossAndDuplication.py`)

Some OGs with ≥4 genes have no resolved gene tree (e.g. all genes from a single species). These caused a `KeyError` in `FindDuplicationsParallel`. Fixed by skipping such OGs in duplication analysis — they cannot have duplication events by definition.

### Bug fix: OGs without resolved gene trees cause ValueError (`AncestralGenome.py`)

The same class of OGs caused a `ValueError` in `ProcessOrthogroupCurrent`. Fixed by checking for the OG's presence in `Resolved_Gene_Trees.txt` before calling `GetAncestralGenes`, returning `None` if absent.

### Bug fix: empty leaf set after pruning causes TreeError (`AncestralGenome.py`)

In `GetAncestralGenes`, the second `gene_tree.prune()` call could receive an empty leaf list when all remaining leaves were flagged as post-duplication duplicates, causing `ete3.coretype.tree.TreeError: Nodes are not connected!`. Fixed by returning `None` when `keep_leaves` is empty.

### New feature: `--min-genes` parameter

The original script hard-codes `min_genes=4`, excluding small OGs from gain/loss analysis entirely. The threshold is now a command-line parameter (`-m`/`--min-genes`, default 4). Use `-m 1` to include all OGs.

### New feature: `--species-tree` parameter

A custom species tree in Newick format (e.g. from IQ-TREE) can be provided via `-s`/`--species-tree`, overriding the OrthoFinder species tree. Internal node labels are added automatically in postorder traversal (root = `N0`, others `N1`, `N2`, ...). Leaf names in the custom tree must match OrthoFinder species names (i.e. the FASTA filenames without extension, e.g. `Schizosaccharomyces_pombe`). All downstream scripts are unaffected — only `ConvertFiles.py` is modified.

### New feature: `--output` parameter

Final results (`GainsLossDuplication/` and `AncestralGenomes/`) are now written to a user-specified directory via `-o`/`--output`. Intermediate files in `WorkingDirectory/GladeWD/` are unaffected.

