# Input files

`raw_counts.csv` is tracked in the repository. The other inputs are not — they
are either too large for git or belong to other publications — so obtain them
from the sources below. `run_full_analysis.R` checks every one of them before
doing any work and names exactly what is missing.

See `REPRODUCIBILITY.md` at the repo root for sources and citations.

```
Input files/
  raw_counts.csv                        strand 01
  mh.all.v2023.2.Mm.symbols.gmt         strand 01   Hallmark
  m2.all.v2023.2.Mm.symbols.gmt         strand 01   curated, all
  m2.cp.v2023.2.Mm.symbols.gmt          strand 01   canonical pathways
  m3.all.v2023.2.Mm.symbols.gmt         strand 01   regulatory targets
  m5.go.bp.v2023.2.Mm.symbols.gmt       strand 01   GO Biological Process
  m8.all.v2023.2.Mm.symbols.gmt         strand 01   cell type signatures
  cross_species/
    Table S2. Results of MS_MS analysis of DLPFC synapse proteomes.xlsx
    Table S3. Results of GSEA of changes in SCZ and BP synapse proteomes.xlsx
    Table S4. Results of GO analyses of module proteins.xlsx
  sc_data_C4OE_PCA.rds                  strand 03   Allen sc reference
```

## Notes

- **`raw_counts.csv`** — line endings were normalised from the original
  CR-only export to standard LF. Gene identifiers are Ensembl IDs.
- **MSigDB `.gmt` files** — mouse collections, **v2023.2**, symbol versions.
  The version matters: a different release changes gene set membership and
  therefore every enrichment result. Download from the
  [MSigDB mouse collections page](https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp).
- **`cross_species/`** — the three Aryal et al. supplementary workbooks
  (human DLPFC synapse proteome). Filenames are matched by pattern
  (`table s2`, `table s3`, `table s4`), so keep the "Table SN" prefix.
- **`sc_data_C4OE_PCA.rds`** — the Allen single-cell cortical reference used by
  strand 03. Treated as optional: if it is absent, strand 03 is skipped and the
  rest of the run still completes.

## What is NOT an input

The mouse DEG table consumed by the cross-species comparison is **not** placed
here. It is strand 01's output
(`Outputs/01_bulk_rnaseq_DE/02_differential_expression/tables/S1_Table_DESeq.csv`),
copied into place automatically by `run_full_analysis.R`, so the cross-species
analysis always uses the DEG table produced by this pipeline.
