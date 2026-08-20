# Scientific assumptions and limitations

## Scope

This document identifies the assumptions required to interpret each event type
reported by Cortar. A reported category describes how an observed alignment
relates to the selected RefSeq transcript; it does not establish the biological
mechanism, pathogenicity or abundance of a complete RNA isoform.

The exact event definitions and equations are given in
[event detection and quantification](event-detection-and-quantification.md).

## Assumptions shared by all event types

### The selected transcript is an appropriate reference

All intron numbers, canonical boundaries and event descriptions are defined by
the selected RefSeq transcript. A junction can be canonical relative to one
transcript and alternative, cryptic or exon-skipping relative to another.
Results can therefore change when transcript choice or annotation release
changes.

Gene-level selection uses transcript records marked canonical in the bundled
annotation. Cortar does not verify that exactly one such transcript exists or
that it is the biologically dominant transcript in the studied tissue.

### Exact coordinate matching is biologically meaningful

Cortar requires exact equality between observed and annotated junction
boundaries. It does not cluster nearby splice sites or allow coordinate
tolerance. Distinct one-base shifts are separate events. Conversely, an exact
coordinate match is accepted without a minimum supporting-read count or an
event-specific quality model.

### The observed junction set is sufficiently complete

An event found in any sample is represented in every sample, with zero assigned
where no supporting read was counted. Absence of support can reflect biological
absence, low coverage, mapping failure or sample-specific technical effects.
Zero is not evidence that the site was adequately assayed.

### The within-intron denominator is an appropriate normalisation

All event types are divided by the sum of SJ counts and the averaged IR score
assigned to the same reference intron. This treats heterogeneous quantities as
components of one local evidence pool. The reported values are not conventional
PSI, transcript fractions or expression-normalised abundance estimates.

Because the values are compositional, adding or removing one event changes the
proportions of all other events in that intron even when their raw counts remain
constant. Differences in junction discovery between cohorts can therefore
change denominators as well as numerators.

### Alignment counts are comparable

The BAM route counts loaded alignments supporting each junction, whereas the
STAR route uses uniquely mapping junction reads. Cortar does not impose a common
minimum mapping quality, anchor length or multi-mapping policy across the two
routes. Results from the two input types should not be assumed interchangeable
without empirical validation.

Duplicate-marked and secondary BAM alignments are excluded, but no UMI-based
deduplication is performed. Primary multi-mapping and supplementary alignments
are not explicitly excluded. Counts can consequently depend on the upstream
aligner, alignment settings and duplicate-marking strategy.

## Canonical splicing

### Assumptions

- Exact agreement with both boundaries of a reference intron is treated as
  canonical junction use.
- The number of supporting alignments is assumed to be proportional to use of
  that junction.
- All competing evidence assigned to the intron is assumed suitable for the
  shared denominator.

### Limitations

- “Canonical” is coordinate-based and does not imply normal function or benign
  clinical interpretation.
- The canonical proportion is not the fraction of all transcripts containing
  the corresponding exon–exon junction. It is only the canonical count divided
  by Cortar's local event-pool total.
- A lower canonical proportion may result from increased cryptic, skipping or
  IR evidence even if the absolute canonical count is unchanged.
- Junction-specific mapping bias, sequence uniqueness and anchor length are not
  modelled.

## Intron retention

### Assumptions

- Coverage spanning four exonic and four intronic bases at a canonical boundary
  is treated as evidence that the boundary was not spliced at that position.
- Averaging the two boundary counts is treated as a suitable single IR score.
- The two boundary measurements are assumed comparable even though they may be
  supported by different molecules and have different sequence mappability.
- The averaged IR score is treated as commensurate with whole-junction read
  counts in the common denominator.

### Limitations of the boundary definition

- Cortar does not measure coverage through the intron body and does not require
  a molecule to connect both boundaries. Local transcription, pre-mRNA,
  incomplete RNA processing or mapping artefacts can produce boundary coverage
  without stable full-intron retention.
- A read can support only one boundary and still increase the IR score by half
  a count. Strongly asymmetric boundary coverage is concealed by the mean unless
  the underlying left and right counts are examined separately; those counts
  are not retained in the final report.
- Short introns, overlapping features and nearby splice sites can cause one
  alignment or read pair to cover both boundary windows.
- The score can be fractional, while SJ evidence is generally integer-valued.
  Combining these quantities assumes that their scales are comparable.

### Split-read contamination of IR evidence

Cortar does not restrict IR boundary counting to alignments without splice
junctions. It tests whether any continuously aligned block covers the complete
eight-base window.

A conventional canonical split read does not cover the boundary window because
its skipped interval begins or ends there. However:

- a cryptic donor inside an intron can leave a continuously aligned block across
  the canonical intron-start boundary;
- a cryptic acceptor inside an intron can leave a continuously aligned block
  across the canonical intron-end boundary;
- a split alignment with a junction elsewhere can also cover a canonical
  boundary continuously.

Such an alignment can contribute to a cryptic SJ count and to one side of the
IR score. The categories are therefore not guaranteed to be supported by
disjoint molecules. This can inflate both a cryptic-junction numerator and the
IR component of the denominator.

### Precomputed IR input

For external IR files, Cortar imports boundary coverage and averages it. It does
not confirm that the external values were calculated over the same eight-base
windows, from non-split reads, or with equivalent alignment filters. The
scientific definition of IR then partly depends on the upstream file-generation
method and must be reported separately.

### Interpretation boundary

The IR result should be described as an **averaged exon–intron boundary-coverage
proportion**, not as full-intron PSI or the percentage of mature transcripts
retaining the intron. Confirmation should use intron-body coverage, sashimi/read
inspection, transcript-aware quantification or an orthogonal assay.

## Cryptic donor and acceptor junctions

### Assumptions

- A junction sharing exactly one selected reference boundary is treated as use
  of a cryptic donor or acceptor, with donor/acceptor identity determined by
  transcript strand.
- Each exact novel coordinate represents a separate event.
- The inferred junction strand is correct.

### Limitations

- “Cryptic” means novel relative to the selected transcript, not necessarily
  absent from all biological annotations. A junction can be present in another
  transcript and simultaneously labelled as a cryptic event relative to the
  selected reference.
- Cortar does not require a canonical splice motif, minimum distance from the
  reference site, minimum read count or minimum junction anchor.
- Nearby cryptic sites are not clustered, so alignment jitter or ambiguous
  placement may fragment one biological signal across several rows.
- Non-canonical motifs can produce unknown strand and prevent donor-versus-
  acceptor classification.
- As described above, a cryptic junction read can also contribute to one IR
  boundary score.
- The proportion is relative to every event assigned to the intron, not merely
  cryptic versus canonical junction counts.

## Exon-skipping junctions

### Assumptions

- A junction connecting boundaries belonging to different reference introns is
  interpreted as skipping the intervening reference exon or exons.
- Intron numbering in the selected transcript accurately represents the exon
  sequence relevant to the sample.
- A long junction is taken as direct evidence of the corresponding skip even if
  alternative transcript structures could produce the same coordinates.

### Limitations

- The event describes one junction, not a complete transcript isoform. It does
  not establish that the flanking exons occur together with the rest of the
  selected transcript.
- The same skipping junction is emitted once for each boundary intron and, when
  `ria` is enabled, for every fully spanned intermediate intron. These rows are
  repeated contexts for one junction, not independent events.
- The raw count is repeated across those rows, but each row has a different
  intron-specific denominator. Their proportions can differ and must not be
  summed or averaged as independent measurements.
- Disabling `ria` does not remove the skipping junction from its two endpoint
  introns; it only removes intermediate-intron assignments.
- No maximum junction length, minimum anchor, motif-quality or read-count filter
  is imposed by Cortar.
- A skipping junction annotated in another transcript may be normal alternative
  splicing rather than aberrant exon loss.

## Alternative-transcript annotation

### Assumptions

- Exact matching to an intron in a non-selected RefSeq transcript of the same
  gene is sufficient to identify known alternative transcript structure.
- The bundled transcript catalogue is appropriate for the tissue and assembly.

### Limitations

- `alternative` is an annotation status, not a mutually exclusive event type.
  The same row can be described as cryptic or exon skipping relative to the
  selected transcript.
- If a junction matches both the selected and another transcript, canonical
  status takes precedence and the alternative match is not evident.
- The final report does not retain the alternative transcript identifier, the
  number of supporting annotations or tissue-specific expression evidence.
- Absence from the bundled alternative transcripts does not establish novelty.

## Junctions with neither reference boundary annotated

### Assumptions

- With `ria` enabled, a long junction fully containing a reference intron is
  considered relevant to that intron even if neither junction boundary matches
  the selected transcript.

### Limitations

- A fully novel junction is not reported solely because it lies inside the
  selected gene. It must share a selected boundary or, with `ria` enabled,
  contain a complete selected intron.
- Short novel intronic or exonic junctions that meet neither condition are
  omitted.
- Qualifying long junctions can be repeated across multiple contained introns,
  again producing non-independent rows with different denominators.
- The label `unannotated junctions` is relative to selected-intron boundaries;
  it is not proof that the event is absent from external annotation databases.

## Approximate frame conservation

### Assumptions

The frame field assumes that a splice-boundary displacement divisible by three
is indicative of frame conservation.

### Limitations

The assessment is a genomic coordinate modulo-three heuristic. It is not
reliably constrained by gene, chromosome, strand, transcript, coding region or
CDS phase. It does not model untranslated regions, start/stop codons, premature
termination or nonsense-mediated decay. IR events receive no frame assessment.

The frame field should not be reported as a predicted protein consequence or
used to infer reading-frame rescue.

## Cross-sample comparisons

### Control means and SD flags

For each event, Cortar calculates the arithmetic mean and sample standard
deviation of control proportions. The test-sample difference is the test
proportion minus this mean. Two-, three- and four-SD flags use strict inequalities
against multiples of the control SD.

These flags are not hypothesis tests. There are no p-values, confidence
intervals, variance stabilisation, minimum event-level coverage requirements or
multiple-testing corrections. With one control the SD is undefined; with zero
control variance any non-zero difference exceeds all SD thresholds. Missing
sample proportions propagate into the mean and SD rather than being removed.

### Control selection

- Default mode excludes samples with the same family or gene label, then applies
  a gene-level median canonical-junction coverage threshold.
- Panel mode excludes the same family only and applies no coverage filter.
- Control selection ignores whether another row is itself labelled as a test,
  so unrelated test samples can act as controls.
- A control can pass the gene-level threshold while having no coverage at the
  specific intron being compared.

The analysis does not model tissue, age, sex, ancestry, batch, RNA integrity,
read length or sequencing protocol. These factors must be controlled in cohort
design.

### Research mode

Research mode calculates the unweighted arithmetic mean and sample SD of the
per-sample event proportions. Each sample has equal weight regardless of read
depth. This is not the same as pooling raw counts across samples. Family
relationships and other covariates are not modelled, and individual sample
values are not retained in the final research TSV.

## Reference and operational limitations affecting scientific interpretation

- Only the bundled hg19 and hg38 RefSeq-derived annotations are used; runtime
  GTF/GFF annotation is not supported.
- The sample table's separate `transcript` column is ignored by the main
  analysis. A transcript must be supplied through the identifier actually used
  for target selection. This can cause the analysed transcript to differ from
  the one recorded in sample metadata.
- `hg38` with the `1000genomes` sequence convention passes initial option
  checking but is not supported by the genome-selection stage.
- Gene-symbol recognition is restricted and can reject valid symbols containing
  punctuation or outside the accepted length range.
- Unversioned transcript identifiers and multiple canonical records can resolve
  ambiguously without explicit rejection.
- The precomputed IR route assumes complete and consistently ordered boundary
  rows; missing or duplicated rows can lead to incorrect boundary pairing.
- Debug intermediate output is currently unreliable, reducing traceability of a
  production run unless lower-level results are captured separately.

## Recommended interpretation boundary

Until event-specific synthetic validation and orthogonal biological benchmarking
are complete, Cortar should be treated as an exploratory junction-prioritisation
tool. Results can identify alignments for manual review and follow-up, but should
not alone support claims of abnormal splicing, transcript abundance, statistical
significance, reading-frame preservation or clinical pathogenicity.
