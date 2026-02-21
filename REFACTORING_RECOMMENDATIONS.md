# Refactoring Recommendations

A series of small, independently testable refactoring changes to incrementally improve the cortar codebase. Each recommendation is self-contained and can be implemented as a single commit/PR.

Recommendations are grouped into tiers from lowest to highest risk. Within each tier, they are ordered by priority.

---

## Tier 1 — Safe Fixes (No Behaviour Change)

These changes fix latent bugs or style violations with zero risk to existing functionality. All can be verified by confirming the existing test suite still passes.

### 1.1 Replace `T`/`F` with `TRUE`/`FALSE` in `5_generate_report.R`

**File:** `R/5_generate_report.R`, lines 32–35

**Problem:** `T` and `F` are not reserved words in R. A user or loaded package could reassign them (e.g. `T <- 0`), causing silent breakage. `TRUE` and `FALSE` are reserved and cannot be overwritten.

**Current code:**
```r
report <- T
splicing_diagnostics_report <- F
full_all_genes_report <- T
figure <- F
```

**Change to:**
```r
report <- TRUE
splicing_diagnostics_report <- FALSE
full_all_genes_report <- TRUE
figure <- FALSE
```

**Also in `4_compare_splicing.R`, line 115:**
```r
decreasing = T
```
Change to:
```r
decreasing = TRUE
```

**And in `5_generate_report.R`, line 167 (`generateReport` default mode column ordering):**
```r
with = F
```
Change to:
```r
with = FALSE
```

Do a project-wide search for ` = T)`, ` = T,`, ` = F)`, ` = F,`, `<- T`, and `<- F` to catch any remaining instances. Exclude test files and documentation.

**Verification:** Run the full test suite (`devtools::test()`). No behaviour change expected.

---

### 1.2 Fix `debug = debug` Default in `cortar_batch()`

**File:** `R/0_run_cortar.R`, line 246

**Problem:** The `cortar_batch()` function signature has `debug = debug` as a default parameter value. This evaluates `debug` at call time in the calling environment, not at definition time. If no `debug` variable exists in the caller's scope, this produces an opaque error. The `cortar()` function uses `debug = FALSE`.

**Current code (line 246):**
```r
cortar_batch <- function(folder,
                           pattern = "*.tsv",
                           ...
                           debug = debug,
                           ria = TRUE){
```

**Change to:**
```r
cortar_batch <- function(folder,
                           pattern = "*.tsv",
                           ...
                           debug = FALSE,
                           ria = TRUE){
```

**Verification:** Write a test that calls `cortar_batch()` in a clean environment without a `debug` argument. It should not error on the function signature. Existing tests continue to pass.

---

### 1.3 Replace Unsafe `seq(1, n)` and `1:length(x)` Patterns

**Problem:** `seq(1, 0)` produces `c(1, 0)` (a two-element descending vector), not an empty sequence. Similarly, `1:0` produces `c(1, 0)`. When these iterate over zero-row data or zero-length vectors, they cause incorrect loop execution and potential errors. `seq_len(n)` and `seq_along(x)` correctly return empty vectors when `n == 0` or `length(x) == 0`.

**Locations to change (each file, with line numbers):**

| File | Line(s) | Current | Replacement |
|------|---------|---------|-------------|
| `R/0_run_cortar.R` | 250 | `seq(1,length(batches_in))` | `seq_along(batches_in)` |
| `R/0_run_cortar.R` | 347 | `seq(1, length(genes))` | `seq_along(genes)` |
| `R/1_select_gene_and_transcript.R` | 78 | `seq(1, length(genes))` | `seq_along(genes)` |
| `R/2_extract_count_reads.R` | 85 | `1:length(sample_names)` | `seq_along(sample_names)` |
| `R/2_extract_count_reads.R` | 125 | `1:length(sample_names)` | `seq_along(sample_names)` |
| `R/2_extract_count_reads.R` | 206 | `1:length(sample_names)` | `seq_along(sample_names)` |
| `R/3_annotate_quantify_events.R` | 28 | `seq(1, nrow(introns))` | `seq_len(nrow(introns))` |
| `R/3_annotate_quantify_events.R` | 132 | `seq(1,nrow(query_intron.dt))` | `seq_len(nrow(query_intron.dt))` |
| `R/3_annotate_quantify_events.R` | 209 | `seq(1,nrow(query_intron.dt))` | `seq_len(nrow(query_intron.dt))` |
| `R/4_compare_splicing.R` | 21 | `seq(1, nrow(Sample_File))` | `seq_len(nrow(Sample_File))` |
| `R/4_compare_splicing.R` | 94 | `seq(1, nrow(all_splicing_events_sample))` | `seq_len(nrow(all_splicing_events_sample))` |
| `R/5_generate_report.R` | 19 | `seq(1, nrow(Sample_File))` | `seq_len(nrow(Sample_File))` |

**Verification:** Run the full test suite. Additionally, write a targeted test that passes an empty character vector to `tx_extraction()` to confirm it returns a zero-row data.table without error (this test exists at `test-1_select_gene_and_transcript.R:46` — confirm it passes after the change).

---

### 1.4 Simplify Mode Validation in `cortar()`

**File:** `R/0_run_cortar.R`, lines 57–62

**Problem:** The current mode validation uses a confusing double-negative nested if/else with an empty `if` body:

```r
if (mode != "default") {
    if (mode == "panel" | mode == "research") {
    } else {
      stop("'", mode, "' is not an available cortar mode ('default','panel','research')")
    }
  }
```

**Change to:**
```r
if (mode %nin% c("default", "panel", "research")) {
    stop("'", mode, "' is not an available cortar mode ('default','panel','research')")
}
```

**Verification:** Run existing tests. The `test-4_compare_splicing.R:125–135` test that passes `mode = "invalid_mode"` validates this path. Optionally add a direct test for `cortar()` with an invalid mode.

---

### 1.5 Remove Commented-Out and Dead Code

**Locations:**

1. **`R/0_run_cortar.R`, line 332:** Remove `# file <- data.table::fread(file)` inside `subsetBamfiles()`. This is leftover from a previous refactoring.

2. **`R/2_extract_count_reads.R`, line 81:** Remove `# isPaired = T` inside the `ScanBamParam` call.

3. **`R/5_generate_report.R`, line 33:** Remove `splicing_diagnostics_report <- F`. This variable is assigned but never read.

4. **`R/5_generate_report.R`, line 44:** Remove `# will need a for loop`. This is a stale TODO.

**Verification:** Run the full test suite. Grep the codebase for `splicing_diagnostics_report` to confirm no other references.

---

## Tier 2 — Structural Improvements (Low Risk)

These changes improve code clarity and reduce fragility. They change function signatures or return types but preserve semantics.

### 2.1 Name the List Returned by `selectGenesTranscripts()`

**File:** `R/1_select_gene_and_transcript.R`, lines 45–50

**Problem:** `selectGenesTranscripts()` returns an unnamed list accessed by numeric index throughout `cortar()`:
```r
genes_tx[[1]]           # genes GRanges
genes_tx[[2]][[1]]      # introns GRanges
genes_tx[[2]][[2]]      # introns data.table
genes_tx[[3]]           # introns other tx GRanges
genes_tx[[4]][[1]]      # intron starts GRanges
genes_tx[[4]][[2]]      # intron ends GRanges
```
If the function's return structure ever changes, all call sites silently break.

**Step 1 — Name the return value in `selectGenesTranscripts()` (line 45–50):**

Change:
```r
return(list(
    genes.GRanges,
    introns.GRanges,
    introns_other_tx.GRanges,
    unlist(introns_jx.GRanges)
  ))
```

To:
```r
return(list(
    genes = genes.GRanges,
    introns = introns.GRanges,
    introns_other_tx = introns_other_tx.GRanges,
    junctions = unlist(introns_jx.GRanges)
  ))
```

**Also name the sub-lists returned by `introns_to_GRange()` (line 198):**

Change:
```r
return(list(introns.GRanges, introns))
```
To:
```r
return(list(granges = introns.GRanges, metadata = introns))
```

**And `introns_jx_to_GRange()` (line 290):**

Change:
```r
return(list(intron_starts.GRanges, intron_ends.GRanges))
```
To:
```r
return(list(starts = intron_starts.GRanges, ends = intron_ends.GRanges))
```

**Step 2 — Update all call sites in `cortar()` (`R/0_run_cortar.R`, lines 165–192):**

Change:
```r
reads <- extractCountReads(
    genes.GRanges = genes_tx[[1]],
    introns.GRanges = genes_tx[[2]][[1]],
    intron_starts.GRanges = genes_tx[[4]][[1]],
    intron_ends.GRanges = genes_tx[[4]][[2]],
    ...
```
To:
```r
reads <- extractCountReads(
    genes.GRanges = genes_tx$genes,
    introns.GRanges = genes_tx$introns$granges,
    intron_starts.GRanges = genes_tx$junctions$starts,
    intron_ends.GRanges = genes_tx$junctions$ends,
    ...
```

And similarly for `annotateQuantifyEvents()`:
```r
events <- annotateQuantifyEvents(
    ...
    introns.GRanges = genes_tx$introns$granges,
    introns_other_tx.GRanges = genes_tx$introns_other_tx,
    introns = genes_tx$introns$metadata,
    ...
```

**Verification:** Run the full integration test (`test-0_cortar.R`). Named list access is backward-compatible with positional access, so nothing should break.

---

### 2.2 Fix `ria == TRUE` Comparison to Idiomatic R

**File:** `R/3_annotate_quantify_events.R`, line 42

**Problem:** `if(ria == TRUE)` is fragile if `ria` is ever `NA` (NA == TRUE returns NA, not FALSE). Use `isTRUE(ria)` or simply `if (ria)`.

**Change:**
```r
if(ria == TRUE){
```
To:
```r
if (isTRUE(ria)) {
```

**Same pattern in other files — search for `== TRUE` and `== FALSE`:**
- `R/0_run_cortar.R`, line 123: `if(debug == TRUE)` → `if (isTRUE(debug))`
- `R/2_extract_count_reads.R`, line 89: `if (paired == FALSE)` → `if (!paired)`
- `R/2_extract_count_reads.R`, line 94: `else if (paired == TRUE)` → `else`
- `R/5_generate_report.R`, line 45: `if(figure == TRUE)` → `if (isTRUE(figure))`
- `R/5_generate_report.R`, line 52: `if (report == TRUE)` → `if (isTRUE(report))`
- `R/5_generate_report.R`, line 61: `if (full_all_genes_report == TRUE)` → `if (isTRUE(full_all_genes_report))`
- `R/99_generate_samplefile.R`, line 104: `if(sex_matched == TRUE)` → `if (isTRUE(sex_matched))`
- `R/99_generate_samplefile.R`, line 127: `if(use.db.paths == TRUE)` → `if (isTRUE(use.db.paths))`

**Verification:** Run the full test suite.

---

### 2.3 Use `|` vs `||` Correctly for Scalar Logical Tests

**Problem:** The `|` operator is vectorised and evaluates both sides. For scalar if-conditions, `||` is correct (short-circuits and is guaranteed scalar).

**Locations:**
- `R/0_run_cortar.R`, line 58: `mode == "panel" | mode == "research"` → `mode == "panel" || mode == "research"` (moot after 1.4 but relevant if 1.4 is deferred)
- `R/0_run_cortar.R`, line 145: `mode == "panel" | mode == "research"` → `mode == "panel" || mode == "research"`
- `R/4_compare_splicing.R`, line 19: `mode == "default" | mode == "panel"` → `mode == "default" || mode == "panel"`
- `R/5_generate_report.R`, line 18: `mode == "default" | mode == "panel"` → `mode == "default" || mode == "panel"`
- `R/5_generate_report.R`, line 114: `mode == "default" | mode == "panel"` → `mode == "default" || mode == "panel"`

Alternatively, consider introducing a helper like `is_comparison_mode <- function(mode) mode %in% c("default", "panel")` in `utils.R` and using it throughout, since this pattern repeats 5 times.

**Verification:** Run the full test suite.

---

## Tier 3 — Extract Testable Helper Functions

These changes break large functions into smaller, independently testable units. They improve test coverage without changing behaviour.

### 3.1 Extract `calculate_unique_events()` from `compareSplicing()`

**File:** `R/4_compare_splicing.R`, lines 94–110

**Problem:** A 17-line nested loop is embedded inside a 200-line function. It identifies splicing events unique to the proband family (present in family, absent in controls). This logic is independently testable but currently can only be tested through the full `compareSplicing()` call.

**Extract into a new function in `R/4_compare_splicing.R`:**

```r
#' Identify Unique Splicing Events
#'
#' Determines which splicing events are present in family members but absent in
#' controls (control average == 0 and family member value != 0).
#'
#' @param events_dt A data.table with splicing event data
#' @param family_pct_cols Character vector of column names for family pct values
#'
#' @return Character vector with entries like "1/2" (unique in 1 of 2 family members)
#'   or "" (not unique)
#' @keywords internal
calculate_unique_events <- function(events_dt, family_pct_cols) {
  unique_col <- character(nrow(events_dt))
  for (i in seq_len(nrow(events_dt))) {
    event_unique_count <- 0
    for (member in family_pct_cols) {
      if (events_dt$controlavg[i] == 0 &
        events_dt[[member]][i] != 0) {
        event_unique_count <- event_unique_count + 1
      }
    }
    if (event_unique_count >= 1) {
      unique_col[i] <- paste0(event_unique_count, "/", length(family_pct_cols))
    } else {
      unique_col[i] <- ""
    }
  }
  return(unique_col)
}
```

**Then replace lines 93–110 in `compareSplicing()` with:**
```r
all_splicing_events_sample$unique <- calculate_unique_events(
    all_splicing_events_sample, familycols
)
```

**New tests to write (`tests/testthat/test-4_compare_splicing.R`):**

```r
test_that("calculate_unique_events identifies unique events correctly", {
  dt <- data.table(
    controlavg = c(0, 0.5, 0),
    pct_sample1 = c(0.3, 0.4, 0),
    pct_sample2 = c(0.2, 0.1, 0)
  )
  result <- calculate_unique_events(dt, c("pct_sample1", "pct_sample2"))
  expect_equal(result, c("2/2", "", ""))
})

test_that("calculate_unique_events handles empty data.table", {
  dt <- data.table(controlavg = numeric(0), pct_s1 = numeric(0))
  result <- calculate_unique_events(dt, "pct_s1")
  expect_equal(result, character(0))
})

test_that("calculate_unique_events handles partial family uniqueness", {
  dt <- data.table(
    controlavg = c(0),
    pct_sample1 = c(0.3),
    pct_sample2 = c(0)
  )
  result <- calculate_unique_events(dt, c("pct_sample1", "pct_sample2"))
  expect_equal(result, c("1/2"))
})
```

**Verification:** Existing `test-4_compare_splicing.R` tests pass. New unit tests pass.

---

### 3.2 Extract `filter_controls_by_coverage()` from `compareSplicing()`

**File:** `R/4_compare_splicing.R`, lines 50–59

**Problem:** The logic for filtering controls by median coverage threshold is inlined within the large `compareSplicing()` function. This is a critical filtering step in `default` mode that determines which controls are used. It cannot be tested in isolation.

**Extract into a new function in `R/4_compare_splicing.R`:**

```r
#' Filter Controls by Coverage Threshold
#'
#' Removes control samples whose median read count at canonical splice sites
#' falls below the coverage threshold for the proband's gene.
#'
#' @param events_dt A data.table with all splicing events
#' @param gene Character string of the proband's gene name
#' @param ctrl_pct_cols Character vector of control pct column names
#' @param ctrl_read_cols Character vector of control read count column names
#' @param coverage_type Character string for coverage type (e.g. "het", "hom", "50")
#'
#' @return A list with two elements:
#'   - `ctrl_pct_cols`: filtered character vector of control pct column names
#'   - `ctrl_read_cols`: filtered character vector of control read count column names
#' @keywords internal
filter_controls_by_coverage <- function(events_dt, gene, ctrl_pct_cols, ctrl_read_cols, coverage_type) {
  canon_splicing_counts <- events_dt[
    gene == gene & annotated == "canonical" & SJ_IR == "SJ",
    ..ctrl_read_cols
  ]
  threshold <- get_coverage_threshold(coverage_type)
  keep <- sapply(canon_splicing_counts, median) > threshold

  list(
    ctrl_pct_cols = ctrl_pct_cols[keep],
    ctrl_read_cols = ctrl_read_cols[keep]
  )
}
```

**Then replace lines 50–59 in `compareSplicing()` with:**
```r
if (mode == "default") {
  filtered <- filter_controls_by_coverage(
    all_splicing_events_sample,
    Sample_File$gene[sample_number],
    ctrlscols,
    ctrlsreadcols,
    Sample_File$coverage[sample_number]
  )
  ctrlscols <- filtered$ctrl_pct_cols
  ctrlsreadcols <- filtered$ctrl_read_cols
}
```

**Note:** The current code on line 51 uses `gene == Sample_File$gene[sample_number]` where `gene` is a column in `events_dt` — verify the scoping is correct when extracting. In data.table, unquoted `gene` in `[` refers to the column. In the extracted function, use `.SD` or pass the gene name as a separate argument and filter like `events_dt[gene_col == gene_name, ...]`.

**New tests to write (`tests/testthat/test-4_compare_splicing.R`):**

```r
test_that("filter_controls_by_coverage removes low-coverage controls", {
  dt <- data.table(
    gene = rep("EMD", 3),
    annotated = rep("canonical", 3),
    SJ_IR = rep("SJ", 3),
    count_ctrl1 = c(100, 120, 110),  # median = 110
    count_ctrl2 = c(10, 5, 8)        # median = 8
  )
  result <- filter_controls_by_coverage(
    dt, "EMD",
    c("pct_ctrl1", "pct_ctrl2"),
    c("count_ctrl1", "count_ctrl2"),
    "het"  # threshold = 60
  )
  expect_equal(result$ctrl_pct_cols, c("pct_ctrl1"))
  expect_equal(result$ctrl_read_cols, c("count_ctrl1"))
})

test_that("filter_controls_by_coverage keeps all controls when threshold is 0", {
  dt <- data.table(
    gene = rep("EMD", 2),
    annotated = rep("canonical", 2),
    SJ_IR = rep("SJ", 2),
    count_ctrl1 = c(5, 3),
    count_ctrl2 = c(2, 1)
  )
  result <- filter_controls_by_coverage(
    dt, "EMD",
    c("pct_ctrl1", "pct_ctrl2"),
    c("count_ctrl1", "count_ctrl2"),
    ""  # threshold = 0
  )
  expect_equal(length(result$ctrl_pct_cols), 2)
})
```

**Verification:** Existing `test-4_compare_splicing.R` tests pass. New unit tests pass.

---

### 3.3 Extract `read_sj_sample()` and `read_bam_sample()` from `extractCountReads()`

**File:** `R/2_extract_count_reads.R`, lines 76–193

**Problem:** `extractCountReads()` is 231 lines with two large branches (`input == "bamfile"` at line 76 and `input == "sj"` at line 124). Each branch handles per-sample reading and processing in a long loop body. These branches should be separate functions.

**Extract the SJ-file reading logic (lines 124–191) into a new function:**

```r
#' Read Splice Junction File for a Single Sample
#'
#' Internal helper that reads a STAR SJ.out.tab file and an intron retention
#' BED file, converts them to GRanges, and filters to the genes of interest.
#'
#' @param sjfile Path to the SJ.out.tab file
#' @param irfile Path to the intron retention BED file
#' @param genes.GRanges GRanges of gene regions for filtering
#' @param introns.GRanges GRanges of intron regions for IR calculation
#'
#' @return A list with `sj` (GRanges of splice junctions) and `ir` (GRanges of intron retention)
#' @keywords internal
read_sj_sample <- function(sjfile, irfile, genes.GRanges, introns.GRanges) {
  # ... extracted loop body from lines 129-191
}
```

**Similarly extract the BAM reading logic (lines 85–122) into:**

```r
#' Read BAM File for a Single Sample
#'
#' @param bamfile Path to the BAM file
#' @param param ScanBamParam object
#' @param paired Logical for paired-end
#' @param stranded Strand mode
#' @param intron_starts.GRanges GRanges for intron start junctions
#' @param intron_ends.GRanges GRanges for intron end junctions
#' @param introns.GRanges GRanges for intron regions
#' @param Genome_Assembly BSgenome object
#'
#' @return A list with `sj` (GRanges) and `ir` (GRanges)
#' @keywords internal
read_bam_sample <- function(bamfile, param, paired, stranded,
                            intron_starts.GRanges, intron_ends.GRanges,
                            introns.GRanges, Genome_Assembly) {
  # ... extracted loop body from lines 86-122
}
```

**New tests for `read_sj_sample()`** (these can use small mock TSV files without needing BAM files):
- Test strand remapping: `"0"` → `"*"`, `"1"` → `"+"`, `"2"` → `"-"`
- Test that BED coordinate adjustment occurs (`start + 1` on the IR file)
- Test with empty SJ/IR files (zero rows)

**Verification:** Run existing integration test (`test-0_cortar.R`). New unit tests for `read_sj_sample()` pass.

---

### 3.4 Extract Genome Assembly Lookup Function

**File:** `R/2_extract_count_reads.R`, lines 59–71

**Problem:** The assembly/annotation lookup is a 13-line nested if/else block. A similar pattern exists in `1_select_gene_and_transcript.R` (lines 31–37) for reference data. The BSgenome lookup specifically could be a utility function.

**Create in `R/utils.R`:**

```r
#' Get BSgenome Object for Assembly and Annotation
#'
#' @param assembly Character: "hg38" or "hg19"
#' @param annotation Character: "UCSC", "1000genomes", or "NCBI"
#'
#' @return A BSgenome object
#' @keywords internal
get_genome_assembly <- function(assembly, annotation) {
  if (assembly == "hg19" && annotation == "UCSC") {
    return(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  } else if (assembly == "hg19" && annotation == "1000genomes") {
    return(BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)
  } else if (assembly == "hg38" && annotation == "UCSC") {
    return(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  } else if (assembly == "hg38" && annotation == "NCBI") {
    return(BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38)
  } else {
    stop("Unsupported assembly/annotation combination: ", assembly, "/", annotation)
  }
}
```

**Verification:** Write tests for valid and invalid combinations. Run existing tests.

---

## Tier 4 — Vectorise Loops (Performance + Clarity)

These changes replace R anti-patterns (row-by-row loops building vectors) with vectorised equivalents. Each is testable by comparing output against the original loop-based version.

### 4.1 Vectorise `eventAnnotation()`

**File:** `R/3_annotate_quantify_events.R`, lines 129–187

**Problem:** The function builds a character vector one element at a time inside a for-loop. This is the classic R anti-pattern: `events[event] <- value` causes repeated vector reallocation. The function has 58 lines of nested if/else that map cleanly to a vectorised approach.

**Approach:** Replace the loop with `data.table::fcase()` or `dplyr::case_when()`. Since the package already depends on `data.table`, use `fcase()`:

```r
eventAnnotation <- function(query_intron.dt) {
  start_jxn <- query_intron.dt$intron_jxn_start
  end_jxn <- query_intron.dt$intron_jxn_end
  strand <- as.character(query_intron.dt$strand)
  sj_ir <- query_intron.dt$SJ_IR

  # Pre-compute range values
  min_range <- pmin(start_jxn, end_jxn, na.rm = FALSE)
  max_range <- pmax(start_jxn, end_jxn, na.rm = FALSE)
  range_diff <- max_range - min_range

  # Build exon skipping labels
  skip_labels <- vapply(seq_len(nrow(query_intron.dt)), function(i) {
    if (!is.na(min_range[i]) && !is.na(max_range[i]) && range_diff[i] > 0) {
      paste0("exon ", paste0(seq(min_range[i] + 1, max_range[i]), collapse = "-"), " skipping")
    } else {
      ""
    }
  }, character(1))

  # Vectorised event classification
  events <- fcase(
    # Both ends annotated
    !is.na(start_jxn) & !is.na(end_jxn) & range_diff == 0 & sj_ir == "SJ",
      paste0("canonical exon ", min_range, "-", max_range + 1, " splicing"),
    !is.na(start_jxn) & !is.na(end_jxn) & range_diff == 0 & sj_ir != "SJ",
      paste0("intron ", min_range, " retention"),
    !is.na(start_jxn) & !is.na(end_jxn) & range_diff > 0,
      skip_labels,

    # Only start annotated
    !is.na(start_jxn) & is.na(end_jxn) & strand == "-",
      paste0("cryptic donor ~ exon ", start_jxn + 1),
    !is.na(start_jxn) & is.na(end_jxn) & strand == "+",
      paste0("exon ", start_jxn, " ~ cryptic acceptor"),
    !is.na(start_jxn) & is.na(end_jxn) & strand == "*",
      "cryptic (strand unknown)",

    # Only end annotated
    is.na(start_jxn) & !is.na(end_jxn) & strand == "-",
      paste0("exon ", end_jxn, " ~ cryptic acceptor"),
    is.na(start_jxn) & !is.na(end_jxn) & strand == "+",
      paste0("cryptic donor ~ exon ", end_jxn + 1),
    is.na(start_jxn) & !is.na(end_jxn) & strand == "*",
      "cryptic (strand unknown)",

    # Neither annotated
    default = "unannotated junctions"
  )

  return(events)
}
```

**Testing strategy:**
1. Keep the old function temporarily renamed as `eventAnnotation_loop()`.
2. Write a test that runs both versions on the mock data from `test-4_compare_splicing.R` and asserts `expect_equal(eventAnnotation(dt), eventAnnotation_loop(dt))`.
3. Once confirmed, remove the loop version.

**New unit tests to add (`tests/testthat/test-3_annotate_quantify_events.R` — new file):**

```r
test_that("eventAnnotation labels canonical splicing correctly", {
  dt <- data.table(
    intron_jxn_start = 5, intron_jxn_end = 5,
    strand = "+", SJ_IR = "SJ"
  )
  expect_equal(eventAnnotation(dt), "canonical exon 5-6 splicing")
})

test_that("eventAnnotation labels intron retention correctly", {
  dt <- data.table(
    intron_jxn_start = 3, intron_jxn_end = 3,
    strand = "+", SJ_IR = "IR"
  )
  expect_equal(eventAnnotation(dt), "intron 3 retention")
})

test_that("eventAnnotation labels exon skipping correctly", {
  dt <- data.table(
    intron_jxn_start = 2, intron_jxn_end = 5,
    strand = "+", SJ_IR = "SJ"
  )
  expect_equal(eventAnnotation(dt), "exon 3-4-5 skipping")
})

test_that("eventAnnotation labels cryptic splice sites correctly", {
  dt <- data.table(
    intron_jxn_start = c(3, 3, NA, NA),
    intron_jxn_end = c(NA, NA, 4, 4),
    strand = c("+", "-", "+", "-"),
    SJ_IR = rep("SJ", 4)
  )
  expected <- c(
    "exon 3 ~ cryptic acceptor",
    "cryptic donor ~ exon 4",
    "cryptic donor ~ exon 5",
    "exon 4 ~ cryptic acceptor"
  )
  expect_equal(eventAnnotation(dt), expected)
})

test_that("eventAnnotation handles empty data.table", {
  dt <- data.table(
    intron_jxn_start = integer(0), intron_jxn_end = integer(0),
    strand = character(0), SJ_IR = character(0)
  )
  expect_equal(eventAnnotation(dt), character(0))
})
```

**Verification:** All new unit tests pass. Existing integration test (`test-0_cortar.R`) produces identical output.

---

### 4.2 Vectorise `framed()`

**File:** `R/3_annotate_quantify_events.R`, lines 200–247

**Problem:** Same row-by-row loop pattern as `eventAnnotation()`. The function checks whether splice junction coordinates are annotated in the RefSeq lookup table and determines if the reading frame is conserved.

**Approach:** Use vectorised `%in%` and `match()` for the lookup operations:

```r
framed <- function(query_intron.dt, assembly) {
  if (assembly == "hg38") {
    rfsq <- refseq_introns_exons_hg38
  } else if (assembly == "hg19") {
    rfsq <- refseq_introns_exons_hg19
  }

  n <- nrow(query_intron.dt)
  frame <- rep("", n)

  is_sj <- query_intron.dt$SJ_IR == "SJ"
  starts <- query_intron.dt$start
  ends <- query_intron.dt$end

  start_annotated <- starts %in% rfsq$region_start
  end_annotated <- ends %in% rfsq$region_end

  for (event in seq_len(n)) {
    if (!is_sj[event]) next

    if (start_annotated[event] && end_annotated[event]) {
      pairend <- unique(rfsq$region_end[rfsq$region_start == starts[event]])
      if (ends[event] %in% pairend) {
        frame[event] <- TRUE
      } else {
        dist2authentic <- abs(ends[event] - pairend)
        frame[event] <- any(dist2authentic %% 3 == 0)
      }
    } else if (start_annotated[event]) {
      pairend <- unique(rfsq$region_end[rfsq$region_start == starts[event]])
      dist2authentic <- abs(ends[event] - pairend)
      frame[event] <- any(dist2authentic %% 3 == 0)
    } else if (end_annotated[event]) {
      pairstart <- unique(rfsq$region_start[rfsq$region_end == ends[event]])
      dist2authentic <- abs(starts[event] - pairstart)
      frame[event] <- any(dist2authentic %% 3 == 0)
    } else {
      frame[event] <- NA
    }
  }
  return(frame)
}
```

**Note:** Full vectorisation of this function is harder because each row's lookup against `rfsq` can return multiple matches. A partial improvement is to precompute the `%in%` checks outside the loop (as shown above) and later convert the inner lookups to keyed data.table joins if profiling shows this is a bottleneck.

**New tests to add (`tests/testthat/test-3_annotate_quantify_events.R`):**

```r
test_that("framed returns TRUE for in-frame canonical junctions", {
  # Use a known in-frame junction from the RefSeq data
  # (Specific coordinates depend on test data available)
})

test_that("framed returns empty string for intron retention events", {
  dt <- data.table(SJ_IR = "IR", start = 100, end = 200)
  result <- framed(dt, "hg38")
  expect_equal(result, "")
})

test_that("framed returns NA for unannotated junctions", {
  dt <- data.table(SJ_IR = "SJ", start = 1, end = 2)
  result <- framed(dt, "hg38")
  expect_equal(result, NA)
})
```

**Verification:** Compare output against the original loop version on mock data. Existing integration test passes.

---

### 4.3 Vectorise Unique Events Detection in `compareSplicing()`

**File:** `R/4_compare_splicing.R`, lines 94–110 (or the extracted `calculate_unique_events()` from 3.1)

**Problem:** If recommendation 3.1 has been implemented, `calculate_unique_events()` still uses a nested loop. This can be vectorised:

```r
calculate_unique_events <- function(events_dt, family_pct_cols) {
  if (nrow(events_dt) == 0 || length(family_pct_cols) == 0) {
    return(character(nrow(events_dt)))
  }

  family_matrix <- as.matrix(events_dt[, ..family_pct_cols])
  is_unique <- events_dt$controlavg == 0 & family_matrix != 0
  counts <- rowSums(is_unique)

  ifelse(counts >= 1,
         paste0(counts, "/", length(family_pct_cols)),
         "")
}
```

**Verification:** Same tests as 3.1 should still pass. Additionally test with larger mock data to verify performance improvement.

---

## Tier 5 — Fix Latent Bugs

### 5.1 Fix `familycols` Reference in Research Mode of `generateReport()`

**File:** `R/5_generate_report.R`, lines 75–83

**Problem:** In the `research` mode branch, `normalSpliceMap()` is called with `familycols[1]` (line 77), but `familycols` is only defined inside the `for` loop in the `default`/`panel` branch (line 28). If research mode reaches this code, it will error with `object 'familycols' not found`.

Currently this is masked because `normalSpliceMap` is only called in research mode from this branch (there's no `figure` gate here unlike default/panel). This is a real bug that will produce an error for any research-mode run.

**Fix:** In research mode, `normalSpliceMap` expects a `familycols` parameter that represents sample pct columns. For research mode, all samples are treated equally:

```r
}else if(mode == "research"){
    all_splicing_events_sample <- comparisons
    testgenes <- unique(all_splicing_events_sample$gene)
    proband <- "splicing_analysis"
    sample_pct_cols <- grep("^pct_", names(all_splicing_events_sample), value = TRUE)
    for(gene in testgenes){
      normalSpliceMap(all_splicing_events_sample,
                      sample_pct_cols[1],
                      proband,
                      ...
```

**Verification:** Write a test or manually run `cortar()` in research mode and verify `normalSpliceMap` is called without error. Check that the generated PDF is correct.

---

### 5.2 Fix Magic Number `minoverlap = 8` in `extractCountReads()`

**File:** `R/2_extract_count_reads.R`, lines 105, 113

**Problem:** The value `8` is a hard-coded minimum overlap for junction read counting. This should be a named constant in `utils.R` with documentation explaining its biological significance.

**Add to `R/utils.R`:**
```r
# Minimum overlap (bp) required for a read to be counted at an exon-intron junction.
# This filters out reads that only marginally overlap the junction boundary,
# ensuring counted reads genuinely span the splice site.
MIN_JUNCTION_OVERLAP <- 8
```

**Replace in `R/2_extract_count_reads.R`:**
```r
minoverlap = MIN_JUNCTION_OVERLAP,
```

**Verification:** Run integration test. No behaviour change.

---

## Tier 6 — Testing Gaps to Fill

These are not code changes but new test files to improve coverage of currently untested functions.

### 6.1 Add Tests for `eventAnnotation()` and `framed()`

**New file:** `tests/testthat/test-3_annotate_quantify_events.R`

See the test cases described in recommendations 4.1 and 4.2. These should be implemented regardless of whether the vectorisation changes are made.

### 6.2 Add Tests for `extractCountReads()` with SJ Input

**New file:** `tests/testthat/test-2_extract_count_reads.R`

Create small mock SJ.out.tab and intron BED files in `tests/testthat/data/` and test:
- Strand remapping logic
- BED coordinate adjustment
- Empty file handling

### 6.3 Add Tests for `generateReport()` Output Files

**New file:** `tests/testthat/test-5_generate_report.R`

Test with mock comparison data in a temporary directory (`withr::with_tempdir()`):
- Excel file is created with expected name
- TSV full report is created
- Column names in output match expectations

---

## Suggested Implementation Order

For a gradual approach with minimal risk:

1. **PR 1:** Tier 1 items (1.1–1.5) — pure cleanup, zero behaviour change
2. **PR 2:** Tier 2 items (2.1–2.3) — naming and idiom improvements
3. **PR 3:** Tier 5 items (5.1–5.2) — fix latent bugs
4. **PR 4:** Tier 3 items (3.1–3.2) + Tier 6.1 — extract helpers and add tests for `compareSplicing`
5. **PR 5:** Tier 3 items (3.3–3.4) + Tier 6.2 — extract helpers and add tests for `extractCountReads`
6. **PR 6:** Tier 4 items (4.1–4.3) — vectorise loops (use new tests from PR 4 to validate)
7. **PR 7:** Tier 6.3 — report generation tests
