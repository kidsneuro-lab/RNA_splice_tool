# One-versus-many analysis: default mode

Default mode performs one comparison for every sample table row whose
`sampletype` is exactly `test`. “Many” means the eligible sample columns used as
controls; there is no matched-pair model or group-level case analysis.

Source: [`compareSplicing()`](../R/4_compare_splicing.R#L86-L238).

## Target selection

All non-blank values from the sample table's `genes` column are passed together
to transcript selection and read extraction. The `transcript` column is ignored.
If different test rows name different genes, the combined event table initially
contains all of them, but each test report is filtered back to that row's gene.

## Family and control sets

For test sample `p`:

```text
family(p)  = every sample with familyID equal to familyID(p), including p
controls(p) = every sample with a different familyID and a different genes value
```

The code accepts either `family` or `familyID`, and either `gene` or `genes`,
within `compareSplicing()`. Control selection does not inspect `sampletype`.
Consequently, an unrelated row marked `test` is used as a control when its gene
value differs. Samples with the same gene label are excluded even when they are
otherwise valid controls.

## Coverage filter

Eligible controls are filtered using the test row's `coverage` value. For the
test gene, cortar takes each control's median raw count across rows satisfying:

```text
gene == test_gene AND annotated == "canonical" AND SJ_IR == "SJ"
```

A control is kept only when its median is **strictly greater than** the
threshold:

| `coverage` value | Threshold |
|---|---:|
| `het` | 60 |
| `hom` or `hemi` | 30 |
| empty string | 0 |
| anything else | `as.numeric(value)` |

The filter is sample-wide for the gene, not event-specific. A value exactly at
the threshold fails. Invalid non-numeric values become `NA` and normally remove
all controls. The full `cortar()` entry point does not validate that the
`coverage` column exists.

## Reported comparison statistics

For every event row `r` and retained control set `C`:

```text
controlavg[r]      = mean(pct[c,r] for c in C)
controlavgreads[r] = mean(count[c,r] for c in C)
controlsd[r]       = sample standard deviation(pct[c,r] for c in C)
controln[r]        = number of retained controls
difference[r]      = pct[p,r] - controlavg[r]
```

`rowMeans()` and `sd()` are called without `na.rm = TRUE`. One control produces
`controlsd = NA`; no controls produce `NA` means and SD with `controln = 0`.

The Boolean flags are strict comparisons:

```text
two_sd   = abs(difference) > 2 * controlsd
three_sd = abs(difference) > 3 * controlsd
four_sd  = abs(difference) > 4 * controlsd
```

These are descriptive thresholds, not hypothesis tests. When all controls have
exactly the same percentage, `controlsd` is zero and any non-zero difference
sets all three flags to `TRUE`.

## `unique` column

For an event with `controlavg == 0`, cortar counts how many family percentage
columns, including the proband, are non-zero. It reports `<count>/<number of
family members>` when at least one is non-zero, otherwise an empty string.

This means `unique` does not mean unique to the proband, nor does it require a
minimum read count. It means observed in at least one family member and absent
from the **mean** of retained controls.

## Sorting and output

Rows are sorted by descending `abs(difference)`. The TSV includes family
percentages and raw counts, control summaries, SD flags, annotation and intron
assignment. An Excel workbook contains the same comparison table with styling.
Default-mode plots are disabled in the current implementation.
