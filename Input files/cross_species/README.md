# Cross-species inputs (strand 02)

Put the three Aryal et al. supplementary workbooks here:

```
Table S2. Results of MS_MS analysis of DLPFC synapse proteomes.xlsx
Table S3. Results of GSEA of changes in SCZ and BP synapse proteomes.xlsx
Table S4. Results of GO analyses of module proteins.xlsx
```

Keep the `Table SN` prefix. `01_validate_inputs.R` finds each file by matching
the lowercased patterns `table s2`, `table s3`, `table s4` against the
filename, so renaming them will break the run.

These are not tracked in git — they belong to another publication. See
`REPRODUCIBILITY.md` for the citation.

## Only the Aryal workbooks go here

The mouse DEG table used by the cross-species comparison is **not** an input —
it is an output of strand 01:

```
Outputs/01_bulk_rnaseq_DE/02_differential_expression/tables/S1_Table_DESeq.csv
```

`run_full_analysis.R` copies it into the strand's working `data_raw/` folder
automatically, so the cross-species comparison always uses the DEG table
produced by this pipeline. There is no need to place any mouse table here.
