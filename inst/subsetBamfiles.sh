#!/usr/bin/env bash
set -e

usage() {
  echo "Usage: $0 --genes <gene1,gene2,...> --assembly <37|38> --cramfiles <path/to/cramfiles.txt> --ref-fasta <path/to/reference.fa> --output <path/to/output_directory> [--temp-dir <path/to/temp_directory>]"
  exit 1
}

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

GENES=""
ASSEMBLY=""
CRAMFILES=""
REF_FASTA=""
OUTPUT_DIR=""
TEMP_DIR=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    -g|--genes)        GENES="$2"; shift 2;;
    -a|--assembly)     ASSEMBLY="$2"; shift 2;;
    -c|--cramfiles)    CRAMFILES="$2"; shift 2;;
    -r|--ref-fasta)    REF_FASTA="$2"; shift 2;;
    -o|--output)       OUTPUT_DIR="$2"; shift 2;;
    -t|--temp-dir)     TEMP_DIR="$2"; shift 2;;
    -h|--help)         usage;;
    --)                shift; break;;
    -*)                echo "Unknown option: $1"; usage;;
    *)                 break;;
  esac
done

# (Optional) basic required-args checks
[ -z "$GENES" ] && echo "Missing --genes" && usage
[ -z "$ASSEMBLY" ] && echo "Missing --assembly" && usage
[ -z "$CRAMFILES" ] && echo "Missing --cramfiles" && usage
[ -z "$REF_FASTA" ] && echo "Missing --ref-fasta" && usage
[ -z "$OUTPUT_DIR" ] && echo "Missing --output" && usage

echo "GENES.     : $GENES"
echo "ASSEMBLY   : $ASSEMBLY"
echo "CRAMFILES  : $CRAMFILES"
echo "REF_FASTA  : $REF_FASTA"
echo "OUTPUT_DIR : $OUTPUT_DIR"

# Check for mandatory arguments
if [[ -z "$GENES" || -z "$ASSEMBLY" || -z "$CRAMFILES" || -z "$REF_FASTA" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: All arguments --genes, --assembly, --cramfiles, --ref-fasta, and --output are required."
    usage
fi

# Validate assembly
if [[ "$ASSEMBLY" != "37" && "$ASSEMBLY" != "38" ]]; then
    echo "Error: Assembly must be either 37 or 38."
    exit 1
fi

# Validate cramfiles.txt
if [[ ! -f "$CRAMFILES" ]]; then
    echo "Error: CRAM files list '$CRAMFILES' does not exist or is not a file."
    exit 1
fi

# Validate REF_FASTA
if [[ ! -f "$REF_FASTA" ]]; then
    echo "Error: Reference FASTA file '$REF_FASTA' does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Output directory '$OUTPUT_DIR' does not exist. Creating it..."
    mkdir -p "$OUTPUT_DIR"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to create output directory '$OUTPUT_DIR'."
        exit 1
    fi
fi

### MODIFICATION ###
# Set temporary directory, if not specified default to ./temp
if [[ -z "$TEMP_DIR" ]]; then
    TEMP_DIR="./temp"
fi

# Create temporary directory if it doesn't exist
mkdir -p "$TEMP_DIR"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to create temporary directory '$TEMP_DIR'."
    exit 1
fi
### END MODIFICATION ###

# Obtain gene coordinates by executing get_genes_coords.R
echo "Obtaining gene coordinates..."
echo "Executing: get_genes_coords.R --genes="$GENES" --hg="$ASSEMBLY" --overhang=1000"
$BASE_DIR/get_genes_coords.R --genes="$GENES" --hg="$ASSEMBLY" --overhang=1000 > $TEMP_DIR/regions.bed
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to obtain gene coordinates."
    exit 1
fi

echo "Gene coordinates obtained:"
echo "--------------------"
cat $TEMP_DIR/regions.bed
echo ""
echo "--------------------"

# Preserve original IFS
OIFS="$IFS"
IFS=$'\n'

# Loop through each CRAM file listed in cramfiles.txt
echo "Starting subsetting of CRAM files..."
for CRAM in $(cat "$CRAMFILES"); do
    echo "Subsetting $CRAM"

    # Check if CRAM file exists
    if [[ ! -f "$CRAM" ]]; then
        echo "Error: CRAM file '$CRAM' does not exist. Exiting."
        IFS="$OIFS"
        exit 1
    fi

    PREFIX=$(basename "$CRAM" .cram)
    echo "Processing prefix: $PREFIX"

    if [[ -f "$OUTPUT_DIR/${PREFIX}_subset.bam" ]]; then
        echo "$OUTPUT_DIR/${PREFIX}_subset.bam already exists. Skipping..."
    else
        # Define temporary and final BAM file paths
        TMP_BAM="$TEMP_DIR/tmp_${PREFIX}_subset.bam"
        FINAL_BAM="$TEMP_DIR/${PREFIX}_subset.bam"

        echo Executing: samtools view for regions "$TEMP_DIR/regions.bed"

        # Subset the CRAM file using samtools view
        samtools view \
            -b \
            --region-file "$TEMP_DIR/regions.bed" \
            -o "$TMP_BAM" \
            -T "$REF_FASTA" \
            "$CRAM"

        # Check if samtools view was successful
        if [[ $? -ne 0 ]]; then
            echo "Error: samtools view failed for '$CRAM'."
            IFS="$OIFS"
            exit 1
        fi

        # Sort the subset BAM file
        samtools sort \
            -@ 4 \
            -O bam \
            -o "$FINAL_BAM" "$TMP_BAM"

        # Check if samtools sort was successful
        if [[ $? -ne 0 ]]; then
            echo "Error: samtools sort failed for '$TMP_BAM'."
            IFS="$OIFS"
            exit 1
        fi

        # Index the sorted BAM file
        samtools index "$FINAL_BAM"

        # Check if samtools index was successful
        if [[ $? -ne 0 ]]; then
            echo "Error: samtools index failed for '$FINAL_BAM'."
            IFS="$OIFS"
            exit 1
        fi

        # Remove temporary BAM file
        rm "$TMP_BAM"

        echo "Successfully subsetted '$CRAM' to '$FINAL_BAM'."
    fi
done

IFS="$OIFS"

# Copy subsetted BAM files to the output directory
echo "Copying subsetted BAM files to '$OUTPUT_DIR'..."
find "$TEMP_DIR" -name '*subset.bam*' -exec cp -v {} "$OUTPUT_DIR" \;

echo "All subsetted BAM files have been copied to '$OUTPUT_DIR'."

# Optional: Clean up temporary directory
# Uncomment the following line if you wish to remove the temp directory after copying
rm -rf "$TEMP_DIR"

echo "Subsetting process completed successfully."

