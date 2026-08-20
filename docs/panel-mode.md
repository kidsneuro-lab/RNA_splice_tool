# Panel mode

Panel mode uses the same per-test-sample comparison and reporting code as
default mode, with three material differences.

1. Genes or RefSeq transcripts come from `genelist`, not the sample table.
2. Every gene in that resolved list is retained in each test sample's report.
3. Controls are all samples with a different family identifier; there is no
   same-gene exclusion and no coverage filter.

Source: [`cortar()`](../R/0_run_cortar.R#L140-L159),
[`compareSplicing()`](../R/4_compare_splicing.R#L116-L238) and
[`generateReport()`](../R/5_generate_report.R#L18-L68).

## Control rule

For test sample `p`:

```text
family(p)   = every sample with familyID equal to familyID(p), including p
controls(p) = every sample with a different familyID
```

Control selection ignores `sampletype`, tissue, batch, sex, phenotype, library
size and gene-level coverage. Other test samples therefore act as controls.

## Statistics

The calculations are identical to default mode:

- row-wise mean sample proportion across controls;
- row-wise mean raw count across controls;
- row-wise sample standard deviation of control proportions;
- test-minus-control-mean difference;
- strict two-, three- and four-SD flags;
- family `unique` ratio.

See [default mode](default-mode.md#reported-comparison-statistics) for the exact
formulae and edge cases.

These values are still calculated independently for each event row. There is no
cross-gene normalisation, panel-wide statistical model, event ranking, gene
burden, pathway analysis, composite panel score or multiple-testing correction.

## Output

For every row marked `sampletype == "test"`, panel mode writes:

- `<prefix><sampleID>_panel_combined_full.tsv`;
- `<prefix><sampleID>_panel_combined_dt_.xlsx`.

The workbook has one worksheet named `panel`, not one sheet per gene. The TSV
and worksheet combine all selected genes. PDF generation is disabled by the
hard-coded report flag.
