# Developer Documentation for cortar

## Overview

This guide provides information for developers who want to contribute to or extend the cortar package. It covers the package architecture, coding standards, testing practices, and contribution guidelines.

## Package Architecture

### Main Components

The cortar package is organized into several key modules:

1. **Pipeline Entry Point** (`0_run_cortar.R`)
   - Main `cortar()` function
   - Batch processing with `cortar_batch()`
   - Parameter validation and workflow orchestration

2. **Gene Selection** (`1_select_gene_and_transcript.R`)
   - Gene and transcript extraction
   - GRanges object creation for genes and introns
   - Junction coordinate calculation

3. **Read Extraction** (`2_extract_count_reads.R`)
   - BAM file reading
   - Split-read and intron retention counting
   - Read aggregation across samples

4. **Event Annotation** (`3_annotate_quantify_events.R`)
   - Splicing event identification
   - Event quantification
   - Annotation with transcript context

5. **Comparison Analysis** (`4_compare_splicing.R`)
   - Sample comparison logic
   - Mode-specific filtering (default, panel, research)
   - Statistical analysis

6. **Report Generation** (`5_generate_report.R`)
   - Excel report creation
   - PDF visualization
   - Formatting and highlighting

7. **Utilities** (`utils.R`)
   - Helper functions
   - Constants and thresholds
   - Input validation

### Data Flow

```
Input: Sample File + BAM Files
    ↓
1. Gene Selection (selectGenesTranscripts)
    ↓
2. Read Extraction (extractCountReads)
    ↓
3. Event Annotation (annotateQuantifyEvents)
    ↓
4. Comparison Analysis (compareSplicing)
    ↓
5. Report Generation (generateReport)
    ↓
Output: Excel + PDF Reports
```

## Coding Standards

### Style Guide

We follow the [Tidyverse Style Guide](https://style.tidyverse.org/) with these key principles:

#### Naming Conventions

- **Functions**: Use `camelCase` or `snake_case` consistently
- **Variables**: Use descriptive names, prefer `snake_case`
- **Constants**: Use `UPPER_SNAKE_CASE`

```r
# Good
DEFAULT_COVERAGE_HET <- 60
sample_count <- length(samples)
extractCountReads <- function(...) { }

# Bad
default.coverage.het <- 60
x <- length(samples)
Extract_Count_Reads <- function(...) { }
```

#### Boolean Values

Always use explicit `TRUE` and `FALSE`, never `T` or `F`:

```r
# Good
if (paired == TRUE) { }
is_valid <- FALSE

# Bad
if (paired == T) { }
is_valid <- F
```

#### Function Documentation

All functions must have roxygen2 documentation:

```r
#' Extract and Count Reads
#'
#' This function extracts split-reads from BAM files and counts them.
#'
#' @param bamfiles Character vector of BAM file paths
#' @param sample_names Character vector of sample names
#' @param assembly Genome assembly ("hg38" or "hg19")
#' @param paired Logical, is the data paired-end?
#'
#' @return A GRanges object with read counts
#' @export
#'
#' @examples
#' \dontrun{
#' reads <- extractCountReads(
#'   bamfiles = c("sample1.bam", "sample2.bam"),
#'   sample_names = c("sample1", "sample2"),
#'   assembly = "hg38",
#'   paired = TRUE
#' )
#' }
extractCountReads <- function(bamfiles, sample_names, assembly, paired) {
  # Function implementation
}
```

### Code Organization

#### Constants

Define all magic numbers as named constants in `utils.R`:

```r
# In utils.R
DEFAULT_COVERAGE_HET <- 60
DEFAULT_COVERAGE_HOM_HEMI <- 30
DEFAULT_COVERAGE_NONE <- 0
INTRON_JUNCTION_UPSTREAM <- 4
INTRON_JUNCTION_DOWNSTREAM <- 3
```

#### Helper Functions

Extract repeated operations into helper functions:

```r
# Instead of repeating GRanges creation:
# Bad
gr1 <- GenomicRanges::GRanges(
  seqnames = chrom1,
  IRanges::IRanges(start = start1, end = end1),
  strand = strand1
)
gr2 <- GenomicRanges::GRanges(
  seqnames = chrom2,
  IRanges::IRanges(start = start2, end = end2),
  strand = strand2
)

# Good - use helper
gr1 <- create_granges(chrom1, start1, end1, strand1, annotation)
gr2 <- create_granges(chrom2, start2, end2, strand2, annotation)
```

#### Error Handling

Use informative error messages:

```r
# Bad
if (!file.exists(file)) {
  stop("File not found")
}

# Good
if (!file.exists(file)) {
  stop(sprintf(
    "File '%s' does not exist or is not readable. Current working directory: '%s'",
    file,
    getwd()
  ))
}
```

Use validation helpers:

```r
# Good
validate_parameter("assembly", assembly, allowed_values = c("hg19", "hg38"))
validate_bam_files(bamfiles)
```

## Testing

### Test Organization

Tests are located in `tests/testthat/` with the naming convention `test-<module>.R`:

- `test-utils.R`: Tests for utility functions
- `test-0_cortar.R`: Tests for main cortar function
- `test-1_select_gene_and_transcript.R`: Tests for gene selection
- `test-4_compare_splicing.R`: Tests for comparison logic

### Writing Tests

Use the `testthat` framework:

```r
library(testthat)
library(cortar)

test_that("function handles valid input correctly", {
  result <- my_function(valid_input)
  expect_true(!is.null(result))
  expect_equal(length(result), 5)
})

test_that("function rejects invalid input", {
  expect_error(
    my_function(invalid_input),
    "expected error message"
  )
})
```

### Test Coverage Goals

- **Utility functions**: 100% coverage
- **Helper functions**: >90% coverage
- **Main pipeline functions**: Integration tests covering common paths
- **Edge cases**: Specific tests for boundary conditions

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-utils.R")

# Run with coverage report
covr::package_coverage()
```

## Adding New Features

### Process

1. **Create an issue** describing the feature
2. **Write tests** for the new functionality (TDD approach)
3. **Implement the feature** following coding standards
4. **Update documentation** (roxygen2 and vignettes if applicable)
5. **Run tests and checks** (`devtools::check()`)
6. **Submit a pull request**

### Example: Adding a New Helper Function

```r
# 1. Add to utils.R
#' Calculate Mean Coverage
#'
#' Internal helper to calculate mean coverage from a vector of read counts.
#'
#' @param counts Numeric vector of read counts
#' @param na_rm Logical, remove NA values?
#'
#' @return Numeric mean coverage
#' @keywords internal
calculate_mean_coverage <- function(counts, na_rm = TRUE) {
  if (length(counts) == 0) {
    stop("Cannot calculate mean of empty vector")
  }
  return(mean(counts, na.rm = na_rm))
}

# 2. Add tests in test-utils.R
test_that("calculate_mean_coverage works correctly", {
  expect_equal(calculate_mean_coverage(c(10, 20, 30)), 20)
  expect_equal(calculate_mean_coverage(c(10, NA, 30), na_rm = TRUE), 20)
  expect_error(
    calculate_mean_coverage(numeric(0)),
    "Cannot calculate mean of empty vector"
  )
})

# 3. Document and export if needed
```

## CI/CD Integration

### GitHub Actions Workflows

The package uses GitHub Actions for continuous integration:

#### Pull Request Workflow (`.github/workflows/pull_request.yml`)

Runs on every pull request:
- Builds Docker image
- Runs test suite
- Checks CLI functionality

#### Release Workflow (`.github/workflows/create_release.yml`)

Runs on releases to main branch:
- Creates tagged releases
- Generates artifacts

### Adding lintr

To add automated style checking:

1. Add to `DESCRIPTION`:
```
Suggests:
    testthat (>= 3.0.0),
    withr,
    lintr
```

2. Create `.lintr` configuration:
```
linters: with_defaults(
  line_length_linter(120),
  object_name_linter = NULL
)
```

3. Add to CI workflow:
```yaml
- name: Run lintr
  run: |
    docker run --entrypoint Rscript --rm rna_splice_tool \
      -e "lintr::lint_package()"
```

### Adding Code Coverage

To add code coverage reporting:

1. Add to `DESCRIPTION`:
```
Suggests:
    testthat (>= 3.0.0),
    withr,
    covr
```

2. Add to CI workflow:
```yaml
- name: Test coverage
  run: |
    docker run --entrypoint Rscript --rm rna_splice_tool \
      -e "covr::package_coverage()"
```

## Performance Considerations

### Profiling

Use `profvis` to identify bottlenecks:

```r
library(profvis)

profvis({
  cortar(
    file = "test_samples.tsv",
    mode = "default",
    assembly = "hg38",
    annotation = "UCSC",
    paired = TRUE,
    stranded = 2,
    output_dir = "output"
  )
})
```

### Optimization Tips

1. **Use data.table** for large data operations
2. **Vectorize** instead of looping when possible
3. **Cache** frequently accessed data structures
4. **Use appropriate** GRanges operations (overlap queries are optimized)
5. **Consider parallel processing** for independent operations

Example of caching:

```r
# Cache genome assembly
.genome_cache <- new.env(parent = emptyenv())

get_genome_assembly <- function(assembly, annotation) {
  cache_key <- paste(assembly, annotation, sep = "_")
  
  if (exists(cache_key, envir = .genome_cache)) {
    return(get(cache_key, envir = .genome_cache))
  }
  
  # Load genome (expensive operation)
  if (assembly == "hg38" && annotation == "UCSC") {
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
  # ... other cases
  
  # Cache and return
  assign(cache_key, genome, envir = .genome_cache)
  return(genome)
}
```

## Debugging

### Debug Mode

Enable debug mode to save intermediate files:

```r
cortar(
  file = "samples.tsv",
  mode = "default",
  assembly = "hg38",
  annotation = "UCSC",
  paired = TRUE,
  stranded = 2,
  output_dir = "output",
  debug = TRUE
)
```

### Using browser()

Add `browser()` to pause execution:

```r
my_function <- function(x) {
  browser()  # Execution pauses here
  result <- process(x)
  return(result)
}
```

### Checking GRanges Objects

```r
# Print structure
str(my_granges)

# View as data.table
as.data.table(my_granges)

# Check overlap operations
findOverlaps(gr1, gr2)
```

## Common Pitfalls

### GRanges Seqlevels

Always set seqlevels style consistently:

```r
# Set UCSC style (chr1, chr2, ...)
GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"

# Or use helper
gr <- create_granges(seqnames, starts, ends, strands, annotation = "UCSC")
```

### File Path Issues

Use absolute paths or check working directory:

```r
# Bad
file.exists("samples/test.bam")

# Good
file.exists("/full/path/to/samples/test.bam")
# Or
file.exists(file.path(getwd(), "samples", "test.bam"))
```

### Memory Management

Clean up large objects when done:

```r
# After processing large alignment object
alignment <- GenomicAlignments::readGAlignmentPairs(bamfile)
# ... use alignment ...
rm(alignment)
gc()  # Force garbage collection
```

## Contributing

### Pull Request Guidelines

1. **Create a feature branch** from `main`
2. **Write or update tests** for your changes
3. **Update documentation** as needed
4. **Run checks**: `devtools::check()`
5. **Ensure tests pass**: `devtools::test()`
6. **Submit PR** with clear description

### Code Review Process

All PRs require review from a maintainer. The reviewer will check:

- Code follows style guide
- Tests are comprehensive
- Documentation is complete
- Changes are minimal and focused
- No breaking changes to public API

### Issue Reporting

When reporting issues, include:

- **Description** of the problem
- **Steps to reproduce**
- **Expected vs actual behavior**
- **Session info** (`sessionInfo()`)
- **Sample data** if possible (anonymized)

## Resources

- [Bioconductor Guidelines](https://bioconductor.org/developers/)
- [Tidyverse Style Guide](https://style.tidyverse.org/)
- [R Packages Book](https://r-pkgs.org/)
- [testthat Documentation](https://testthat.r-lib.org/)
- [roxygen2 Documentation](https://roxygen2.r-lib.org/)

## Contact

For questions or discussions:
- Open an issue on GitHub
- Contact the maintainer (see DESCRIPTION file)
