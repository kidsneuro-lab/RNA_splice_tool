# Research mode

Research mode aggregates every sample in the sample table into one descriptive
cohort. It does not use `sampletype`, `familyID`, the `genes` column or coverage
thresholds when comparing samples. Target genes or RefSeq transcripts come from
`genelist`.

Source: [`compareSplicing()`](../R/4_compare_splicing.R#L239-L282) and
[`generateReport()`](../R/5_generate_report.R#L69-L89).

## Cohort statistics

For each event row `r` and all sample rows `S`:

```text
controlavg[r]      = mean(pct[s,r] for s in S)
controlavgreads[r] = mean(count[s,r] for s in S)
controlsd[r]       = sample standard deviation(pct[s,r] for s in S)
controln[r]        = number of samples
```

The `control` prefix is retained even though all samples are peers. The mean and
SD do not remove missing values. With one sample, `controlsd` is `NA`.

No per-sample difference, outlier flag, z-score, p-value, confidence interval,
minimum prevalence or multiple-testing correction is calculated. There is also
no modelling of population structure, repeated families, batch or other sample
metadata.

## Sorting and data retention

Rows are sorted by ascending `gene` and then ascending `controlavg`. The final
table retains these columns:

- assembly and genomic coordinates;
- gene, event, exact annotation and approximate frame annotation;
- cohort mean proportion, SD, sample count and mean raw read count;
- intron number and `SJ_IR` evidence class.

Individual `pct_<sampleID>` and `count_<sampleID>` columns are removed from the
returned comparison table and final TSV. They are present only in the earlier
event table, which can currently be captured through internal calls; the debug
pipeline is affected by a separate defect described in
[assumptions and limitations](assumptions-and-limitations.md).

## Output

Research mode writes:

- `<prefix>splicing_analysis_combined_full.tsv`;
- `<prefix>splicing_analysis_<gene>_normalSpliceMap.pdf`;
- `<prefix>splicing_analysis_<gene>_normalSpliceMap_bar.pdf`.

The first plot shows canonical SJ cohort means with a mean ± two-SD ribbon. It
does not show individual samples. In the second plot, research mode sets
`difference` to zero, so its bars carry no cohort contrast; the plot is not an
outlier display.

## Interpretation

`controlavg` is the average of per-sample, within-intron event shares. It is not
the pooled count proportion:

```text
mean(count_event / sample_intron_total)
```

is generally different from:

```text
sum(count_event) / sum(sample_intron_total)
```

Every sample therefore receives equal weight in `controlavg` regardless of its
read depth, while `controlavgreads` separately reports the unnormalised mean
count.
