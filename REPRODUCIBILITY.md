# Reproducing this analysis

The analysis runs from one entry point:

```r
source("run_full_analysis.R")
```

That creates the output tree, checks the inputs, and runs all three strands in
dependency order. It restores your working directory when it finishes, and
writes a provenance record to `Outputs/_provenance/`.

The manuscript figures are built separately, afterwards:

```r
source("Figures_paper.R")
```

`Figures_paper.R` only *reads* the caches and result tables the pipeline wrote,
so it never re-runs any analysis and cannot change a number. Keeping it out of
`run_full_analysis.R` is deliberate: figure formatting is iterated far more
often than the analysis, and a figure edit should never require a re-run.

---

## 1. Quick start

```bash
git clone <repo-url>
cd Git_hub_C4-OE-CorticalRNAseq
```

```r
# 1. packages -- pinned versions if a lockfile exists, otherwise bootstrap
renv::restore()                                    # if renv.lock is present
# or
source("scripts/00_setup/install_dependencies.R")  # if it is not

# 2. confirm the environment
source("scripts/00_setup/check_environment.R")

# 3. obtain the input data (section 3 below), then
source("run_full_analysis.R")

# 4. build the manuscript figures and supplementary tables
source("Figures_paper.R")
```

Expect roughly 10–15 minutes for the analysis (measured on an Apple-silicon
laptop: strand 01 ~5 min, strand 02 ~2 min, strand 03 ~4 min), dominated by
WGCNA, the GSEA permutations and loading the single-cell reference.
`Figures_paper.R` adds about a minute.

To run one strand at a time:

```r
RUN_STRANDS <- c("01")
source("run_full_analysis.R")
```

---

## 2. Repository layout

`Outputs/` mirrors `scripts/` one-for-one, so every result has an obvious
matching piece of code.

```
scripts/
  00_setup/                          dependency install, environment check, renv
  01_bulk_rnaseq_DE/                 QC -> DESeq2 -> GO/GSEA -> WGCNA
  02_cross_species_comparison/       mouse C4-OE vs human SCZ synapse proteome
  03_scRNAseq_reference_projection/  DEG projection onto an Allen sc reference

run_full_analysis.R                  runs the three strands
Figures_paper.R                      manuscript figures + supplementary tables,
                                     from cached results only

Outputs/
  01_bulk_rnaseq_DE/
    01_qc/{tables,figures}
    02_differential_expression/{tables,figures}
    03_enrichment/{tables,figures}
    04_wgcna/{tables,figures}
    rdata/                           per-step caches, so a step can re-run alone
  02_cross_species_comparison/
    data_raw/                        assembled by run_full_analysis.R
    data_processed/
    reference/                       frozen Ensembl ortholog table
    results/{tables,figures,logs}
  03_scRNAseq_reference_projection/{tables,figures}
  paper_figures/                     Figures_paper.R
    supplementary/                   supplementary figures and .xlsx tables
    projection_panels/               strand-03 panels, written by Projection.R
  _provenance/                       sessionInfo, timings, output manifests
```

`Outputs/` is git-ignored: it is fully regenerable, and tracking it invites the
results on disk to drift from the code that made them. If you want the final
supplementary tables in version control, un-ignore those specific files.

**Figure code is separate from analysis code.** `Figures_paper.R` registers each
figure in `PAPER_FIGURES` at the bottom of the file; set `FIGURES_TO_MAKE` before
sourcing to rebuild a subset:

```r
FIGURES_TO_MAKE <- "figure6"
source("Figures_paper.R")
```

The one exception is the strand-03 projection panels, which are written by
`Projection.R` during the run itself rather than re-drawn afterwards, so that
they are exactly the plots the analysis produced.

**Strand 02 is run by pointing the working directory at its output folder.**
Its 15 scripts use paths relative to the working directory, so this routes
their output without editing ~2,000 lines. `run_full_analysis.R` handles it and
restores the previous working directory afterwards.

**Each strand runs in its own R process.** Strand 01 attaches DESeq2,
AnnotationDbi, clusterProfiler and WGCNA, which export functions named `count`,
`select`, `filter`, `rename`, `first` and `slice`. Strand 02 is written in dplyr
and calls those verbs unqualified 107 times. Since `library(dplyr)` does nothing
when dplyr is already attached, dplyr is never promoted back to the front of the
search path and strand 02 silently calls the wrong function — the observed
failure was:

```
Error in count(., gene_symbol, name = "n_protein_entries") :
  Argument 'x' is not a vector: list
```

Running each strand in a clean process removes the entire class of problem and
makes each strand independently reproducible rather than dependent on what ran
before it in the session. `ISOLATE_STRANDS <- FALSE` disables it for debugging.

---

## 3. Input data

Four groups of files are needed. Only `raw_counts.csv` is in the repository; the
rest are not, being either too large for git or the property of another
publication.

| File | Goes in | Source |
|---|---|---|
| `raw_counts.csv` | `Input files/` | This study. Gene × sample counts, Ensembl IDs. Deposit with GEO accession before publication. |
| `m*.gmt`, `mh.all.*.gmt` (6 files) | `Input files/` | [MSigDB mouse collections](https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp), **v2023.2**, symbol versions |
| `Table S2/S3/S4 *.xlsx` | `Input files/cross_species/` | Supplementary tables of Aryal et al., human DLPFC synapse proteome |
| `sc_data_C4OE_PCA.rds` | `Input files/` | Allen Institute single-cell reference, as a Seurat object |

`run_full_analysis.R` checks all of these before doing any work, and names
exactly what is missing rather than failing partway through.

The mouse DEG table that strand 02 consumes is **not** an input — it is strand
01's output, copied into place automatically, so the cross-species comparison
always uses the DEG table produced by this pipeline.

---

## 4. What is pinned, and why

Reproducibility here means another group gets the same numbers, not merely a
similar story. Four things had to be pinned for that to hold.

### Package versions — `renv.lock`

Generated by `scripts/00_setup/snapshot_renv.R` on the machine that produced
the published results, then committed. Others run `renv::restore()`.

Versions that genuinely move numbers, rather than cosmetics:
`DESeq2` (dispersion estimation and independent filtering), `fgsea` (GSEA
p values), `clusterProfiler`, `org.Mm.eg.db`/`org.Hs.eg.db` (symbol mapping,
which sets the tested universe), `GO.db` (term membership), `msigdbr` (gene set
contents), `WGCNA` (module detection).

### Ensembl release — fixed at 116

`scripts/02_cross_species_comparison/R/04_map_human_mouse_orthologs.R` sets
`ENSEMBL_RELEASE <- 116`, so the human–mouse ortholog set does not depend on
whatever release Ensembl happens to serve on a given day. That set determines
the shared measurable universe, every Fisher test denominator, and every GSEA
gene-set size, so pinning it is essential for identical numbers across machines.
The methods text quotes release 116 accordingly.

### The ortholog table itself — frozen to CSV

Pinning the release is not quite enough, because it still needs the network.
On the first run, R/04 writes
`scripts/02_cross_species_comparison/reference/ensembl_human_mouse_orthologs_frozen.csv`.
**Commit that file.** When it is present, Ensembl is never contacted and the
gene universe is fixed for everyone. To refresh it deliberately, delete it,
re-run, and commit the new one with a note saying why.

### Random seeds

`set.seed(42)` in `00_config.R` for strand 01; `set.seed(20260723)` before each
`fgsea` call in strand 02; `set.seed(1234)` in strand 03. GSEA and WGCNA are
the stochastic steps.

---

## 5. Known limits on reproducibility

Honest caveats, all of which belong in the manuscript rather than being
discovered by a reader.

**MSigDB `.gmt` files are inputs, not downloads.** Version 2023.2 is recorded
in `check_environment.R`. Using a different release changes gene set membership
and therefore every enrichment result.

**`msigdbr` is called live in strand 02** (scripts 06, 07, 08, 10, 13) rather
than reading the `.gmt` files. Its gene sets come from the installed package
version, which `renv.lock` pins — but note that strands 01 and 02 therefore
draw gene sets from two different sources. Worth unifying.

**WGCNA on n = 8, and how the power was chosen.** `WGCNA_POWER <- 15` is fixed
rather than auto-selected, because on 8 samples no criterion-based choice
exists: signed R² never reaches the conventional 0.90 threshold at any power
from 1 to 20, peaking at 0.875 (power 20). Powers 13–20 all clear 0.80, so that
criterion alone does not discriminate between them either. Power 15 was chosen
as the most *representative* partition of that range — across powers 12–18 it
has the highest mean adjusted Rand index against the other powers (0.356, versus
0.289–0.339 for the rest), and it clears R² = 0.80 with margin (0.833). It
yields 25 modules with 203 of 15,729 genes unassigned.

`WGCNA_MIN_MODULE_SIZE` is likewise fixed, at 200 rather than WGCNA's default
of 30, so that modules are large enough to support GO annotation at this sample
size.

**Two WGCNA modules survive BH correction**, both strongly correlated with
genotype: `blue` (1,819 genes, r = +0.993, padj = 2.5e-05) and `turquoise`
(2,234 genes, r = −0.985, padj = 1.1e-04). These two are themselves mutually
correlated and represent a single dominant axis of variation, not two
independent findings.

**Module membership is not stable across powers, even though the themes are.**
The median adjusted Rand index between any two powers in 12–18 is 0.31, and the
gene-level (size-weighted) mean module stability is 0.243. In practice the
biological *themes* reproduce across powers while individual gene-to-module
assignments do not, so module-level claims are safe and gene-level ones are not.
In particular, **C4b's module assignment is not reproducible**: across powers
12–18 it lands in `blue` four times, `turquoise` twice and `pink` once, with a
gene-level stability of 0.38. No claim should rest on which module C4b falls
into. Its assignments at each power are in
`wgcna_power_sensitivity_assignments.csv`. `04b_wgcna_power_sensitivity.R`
regenerates all of these diagnostics; it is deliberately excluded from
`run_all.R` because it rebuilds one TOM per power (~30 min).

**Sex is confounded with genotype by design.** All four C4-OE animals are male;
Control is 2F/2M, so genotype and sex correlate at r = 0.577. The design
`~ sex + genotype` adjusts for sex using the within-Control sex contrast, but it
cannot fully separate a genotype effect from a male-specific one. Sex was
*inferred* from Xist/Ddx3y counts, not recorded.

**Cross-species pathway selection is post-hoc.** The RNA-localization result is
selected in `R/08` from a hardcoded list of pathway names derived from the same
GSEA that identified it, and `R/09` describes those modules as "prespecified".
They were not. Either genuinely pre-register the module list or describe the
analysis as exploratory.

**Gene-level cross-species overlap is weak, and directionally random.** The
primary comparison (human FDR < 0.10, mouse padj < 0.05) gives 43 observed vs
36.2 expected genes, enrichment 1.19, Fisher p = 0.13 — not significant. The
stricter comparison (human FDR < 0.05) does reach p = 0.0023, but it is the more
discordant of the two: 16 of its 23 genes move in opposite directions across
species. Across the full 43-gene overlap only 19 agree in direction (44%,
binomial p = 0.54), which is what chance would give. The cross-species claim
rests on pathway-level convergence, not on shared genes.

**Pathway-level NES correlations are anticonservative.** `R/09` and `R/12`
Spearman-correlate human against mouse NES across 2,739 GO:BP terms. GO terms
share large fractions of their gene membership, so those observations are
heavily non-independent and the p value will be tiny regardless of the biology.
Report rho descriptively; drop the p value or use a permutation null.

---

## 6. Verifying a run matched

Every run writes to `Outputs/_provenance/`:

- `sessionInfo_<timestamp>.txt` — full package versions
- `run_summary_<timestamp>.txt` — R version, platform, strands, timings
- `output_manifest_<timestamp>.csv` — every file produced, with size

To compare two runs, diff the manifests. Identical sizes across every table is
a strong signal; if they differ, `check_environment.R` output from both runs
usually identifies the cause immediately.

---

## 7. Current status

All three strands run end-to-end from `run_full_analysis.R` and complete without
manual intervention, and `Figures_paper.R` then builds the manuscript figures and
supplementary tables from the cached results. Strand 03 requires the single-cell reference
(`sc_data_C4OE_PCA.rds`); when it is absent, `run_full_analysis.R` skips strand 03
with an explanation and strands 01–02 still complete.
