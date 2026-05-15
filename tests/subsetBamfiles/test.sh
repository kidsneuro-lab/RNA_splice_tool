#!/bin/bash

set -euo pipefail

echo $PWD

echo "🏃 Running subsetBamfiles test script"

echo "📂 Output directory: $PWD/tests/subsetBamfiles/output"
echo "📂 Temp directory: $PWD/tests/subsetBamfiles/temp"
echo "📂 Reference fasta: $PWD/tests/subsetBamfiles/data/chr21.fasta.gz"
echo "📂 Alignment files: $PWD/tests/subsetBamfiles/data/samples.tsv"
echo "🧬 Genes: ETS2"

rm -rf $PWD/tests/subsetBamfiles/output
rm -rf $PWD/tests/subsetBamfiles/temp
mkdir -p $PWD/tests/subsetBamfiles/output
mkdir -p $PWD/tests/subsetBamfiles/temp

./inst/subsetBamfiles.sh \
  --alignmentfiles $PWD/tests/subsetBamfiles/data/samples.tsv \
  --genes ETS2 \
  --assembly "38" \
  --ref-fasta $PWD/tests/subsetBamfiles/data/chr21.fasta.gz \
  --output $PWD/tests/subsetBamfiles/output \
  --temp-dir $PWD/tests/subsetBamfiles/temp

echo "Obtain number of lines in the input bam file corresponding to ETS2 gene"
ORIGINAL_LINE_COUNT=$(samtools view $PWD/tests/subsetBamfiles/data/SAMPLE1.bam chr21:38804182-38825955 | wc -l)

SAMPLE1_SUBSET_LINE_COUNT=$(samtools view $PWD/tests/subsetBamfiles/output/SAMPLE1_subset.bam | wc -l)
SAMPLE2_SUBSET_LINE_COUNT=$(samtools view $PWD/tests/subsetBamfiles/output/SAMPLE2_subset.bam | wc -l)

echo "Original line count for SAMPLE1.bam: $ORIGINAL_LINE_COUNT"
echo "Subset line count for SAMPLE1.bam: $SAMPLE1_SUBSET_LINE_COUNT"
echo "Subset line count for SAMPLE2.bam: $SAMPLE2_SUBSET_LINE_COUNT"

# Check that both sample 1 and sample 2 have the same number of lines as the original bam file
if [[ "$SAMPLE1_SUBSET_LINE_COUNT" -eq "$ORIGINAL_LINE_COUNT" && "$SAMPLE2_SUBSET_LINE_COUNT" -eq "$ORIGINAL_LINE_COUNT" ]]; then
  echo "✅ Test passed: Both subset bam files have the same number of lines as the original bam file."
else    
  echo "❌ Test failed: Subset bam files do not have the same number of lines as the original bam file."
  exit 1
fi