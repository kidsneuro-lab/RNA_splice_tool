# Scripts testing

## Test data preparation for subsetBamfiles.sh

### Preparation of test data

1. Download ENCODE bam file

```bash
wget -qc https://www.encodeproject.org/files/ENCFF433OOU/@@download/ENCFF433OOU.bam
samtools sort -@ 8 -o ENCFF433OOU.sorted.bam ENCFF433OOU.bam
samtools index ENCFF433OOU.sorted.bam
```

2. Download Chromosome 21
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr2.fa.gz
gunzip chr21.fa.gz
mv chr21.fa chr21.fasta
bgzip chr21.fasta
samtools faidx chr21.fasta
```

3. Subset original BAM file to chr21:37900354-41084017 and call it SAMPLE1.bam

```bash
samtools view -o SAMPLE1.bam -T chr21.fasta.gz ENCFF433OOU.sorted.bam chr21:37900354-41084017
```

4. Subset original BAM file to chr21:37900354-41084017 and call it SAMPLE2.cram

> [!NOTE]
>
> The awk script portion below is required to ensure reads are corrected 'formatted'
> to ensure they are compliant with a CRAM file containing a single chromosome in the header

```bash
samtools view -h \
  ENCFF433OOU.sorted.bam \
  chr21:37900354-41084017 |
awk 'BEGIN{OFS="\t"} /^@SQ/ { if ($0 ~ /SN:chr21(\t|$)/) print; next } /^@/ { print; next } { print }' |
samtools view -b -o SAMPLE2.chr21.onlyheader.bam

samtools view \
  -O CRAM \
  -o SAMPLE2.cram \
  -T chr21.fasta.gz \
  SAMPLE2.chr21.onlyheader.bam

rm SAMPLE2.chr21.onlyheader.bam
```

5. samples.tsv

```csv
tests/subsetBamfiles/data/SAMPLE1.bam
tests/subsetBamfiles/data/SAMPLE2.cram
```