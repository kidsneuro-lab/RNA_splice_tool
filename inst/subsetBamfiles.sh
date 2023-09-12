#!/bin/sh

rm ./temp/*

OIFS="$IFS"
IFS=$'\n'

for CRAM in `cat cramfiles.txt`; do
    echo "Subsetting $CRAM"

    PREFIX=`basename $CRAM .cram`
    echo "$PREFIX"

    samtools view \
        -@ 4 \
        -b \
        -o "temp/tmp_${PREFIX}_subset.bam" \
        -T #replace this tag with the reference .fasta# \
        $CRAM #replace this tag with gene coordinates#

    samtools sort \
        -@ 4 \
        -O bam \
        -o "temp/${PREFIX}_subset.bam" "temp/tmp_${PREFIX}_subset.bam"

    samtools index "temp/${PREFIX}_subset.bam"
    rm "temp/tmp_${PREFIX}_subset.bam"

    if [ $? -ne 0 ]; then
        echo "Encountered error while subsetting $CRAM"
        IFS="$OIFS"
        exit 1
    fi
done

IFS="$OIFS"

find ./temp -name '*subset.bam*' -exec cp -v {} "#replace this tag with /path/to/destination/directory#" \;