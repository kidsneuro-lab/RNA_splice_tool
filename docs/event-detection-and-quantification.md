# Event detection and quantification

## Purpose of this document

This document describes the scientific rules used by Cortar to detect, classify
and quantify RNA-splicing events. It is written for scientists interpreting
Cortar results. The central point is that Cortar reports an event's **share of
all evidence assigned to a particular reference intron**. This quantity is not
conventional percent-spliced-in (PSI).

## Reference model and terminology

Cortar analyses selected assembly-specific RefSeq transcript records. A
gene-level request uses record(s) marked canonical in the bundled annotation,
whereas a transcript accession selects the corresponding record. Their introns
are the reference introns used to define canonical boundaries, intron numbers
and event categories. Introns from non-selected RefSeq transcripts of the same
gene are retained only to identify exact matches to annotated alternative
transcript structures.

For reference intron \(i\), let:

- \(s_i\) be its genomic start coordinate;
- \(e_i\) be its genomic end coordinate;
- \(E_i\) be the set of events assigned to that intron;
- \(C_{ijr}\) be the raw evidence count for event \(r\), sample \(j\) and
  intron \(i\).

Splice-junction coordinates describe the genomic interval removed by the split
alignment. Event classification depends on exact equality between the observed
junction boundaries and reference intron boundaries.

## Evidence used by Cortar

### Splice-junction evidence

For BAM input, a splice junction is detected when an RNA-seq alignment contains
a skipped genomic interval. The raw count for a junction is the number of loaded
alignments containing that exact junction. Duplicate-marked and secondary
alignments are excluded. Cortar applies no additional minimum mapping-quality,
uniqueness or junction-anchor threshold; primary multi-mapping alignments and
supplementary alignments are not explicitly removed.

For precomputed STAR junction input, the raw count is STAR's uniquely mapping
read count for that exact junction. STAR's multi-mapping count and maximum
overhang fields are not used by Cortar.

Counts produced by the BAM and STAR-input routes are therefore not guaranteed to
be numerically equivalent.

### Intron-retention evidence

Intron retention is not measured from coverage across the entire intron.
Instead, Cortar estimates retention separately at the two exon–intron
boundaries of each reference intron.

Each boundary window contains eight genomic bases:

```text
Intron-start window = [intron start − 4, intron start + 3]
Intron-end window   = [intron end − 3, intron end + 4]
```

Each window therefore contains four bases on the exonic side and four bases on
the intronic side of its boundary. A read or read-pair alignment contributes to
a boundary count only if a continuously aligned block covers all eight bases.
Bases skipped by a splice junction do not satisfy the overlap requirement.

Let \(L_{ij}\) and \(R_{ij}\) be the counts at the two boundary windows for
intron \(i\) in sample \(j\). The raw intron-retention evidence is:

\[
C^{IR}_{ij} = \frac{L_{ij} + R_{ij}}{2}.
\]

This value can be fractional. The reads counted at the two boundaries need not
be the same molecules, and no read is required to traverse the intron body.

All loaded alignments are eligible for boundary counting; Cortar does not first
select only alignments without splice junctions. A conventional split read whose
gap begins or ends at the canonical boundary will not span that eight-base
window and will not count there. However, a split read with a cryptic splice
site inside the intron can have a continuously aligned block that spans the
opposite canonical boundary. Such a read can support both a cryptic junction and
one side of the IR score.

For example, if 12 alignments cover the intron-start window and 8 cover the
intron-end window:

```text
IR count = (12 + 8) / 2 = 10
```

For precomputed IR input, Cortar imports the supplied coverage value at each
boundary and applies the same arithmetic mean. The eight-base definition then
depends on how the external IR file was produced; Cortar does not reconstruct
or verify it.

## Assignment of events to reference introns

An observed event with genomic interval \([s,e]\) is assigned to reference
intron \(i\) when at least one of the following applies:

1. \(s=s_i\): the event shares the reference intron's genomic start boundary;
2. \(e=e_i\): the event shares the reference intron's genomic end boundary;
3. when reads-in-absentia (`ria`) is enabled, \([s_i,e_i]\) lies completely
   inside \([s,e]\).

The comparison is strand-aware. An unknown event strand can match either
strand.

The third rule assigns a long exon-skipping junction to intermediate introns
that it spans without sharing either boundary. Disabling `ria` removes these
intermediate assignments, but the junction remains assigned to the first and
last introns whose boundaries it shares.

This intron-centred assignment has two consequences:

- the same genomic junction may appear in several output rows, once for each
  reference intron to which it is assigned;
- its percentage can differ between those rows because each intron has a
  different event pool and denominator.

## Common quantification rule

Every event category uses the same within-intron denominator. For sample \(j\)
and reference intron \(i\), the total evidence is:

\[
T_{ij} = \sum_{q \in E_i} C_{ijq}.
\]

The reported proportion for event \(r\) is:

\[
P_{ijr} = \frac{C_{ijr}}{T_{ij}}.
\]

The event pool can contain:

- the canonical splice junction;
- cryptic donor and acceptor junctions;
- exon-skipping junctions assigned to the intron;
- exact junctions from alternative transcripts;
- qualifying junctions with neither boundary annotated in the selected
  transcript;
- the averaged IR boundary score.

When \(T_{ij}=0\), Cortar reports all event proportions for that intron and
sample as zero. When \(T_{ij}>0\), the proportions within that sample–intron
pool sum to one, subject to numerical precision.

This calculation is not adjusted for library size, gene expression, read
length, insert size, mappability, GC content or the number of possible junctions.
It is also compositional: adding evidence for one event decreases the reported
proportion of every other event in the same pool even when their raw counts do
not change.

## Event-specific detection and calculation

### Canonical splicing

#### Detection

An observed splice junction is labelled canonical when both of its boundaries
exactly match the start and end of the same reference intron.

#### Raw count

The canonical raw count is the number of alignments, or STAR uniquely mapping
reads, supporting that exact junction.

#### Reported proportion

\[
P^{canonical}_{ij} = \frac{C^{canonical}_{ij}}{T_{ij}}.
\]

The denominator is not limited to canonical versus retained-intron evidence; it
contains all qualifying event types listed above.

#### Interpretation

The result is the canonical junction's share of Cortar's local event evidence.
It is not the fraction of all gene transcripts using that junction, nor a direct
measure of exon inclusion.

### Intron retention

#### Detection

An IR row is generated for every selected reference intron. Its evidence is the
mean count across the two eight-base exon–intron boundary windows, including
four exonic and four intronic bases at each boundary.

#### Raw count

\[
C^{IR}_{ij} = \frac{L_{ij} + R_{ij}}{2}.
\]

#### Reported proportion

\[
P^{IR}_{ij} = \frac{(L_{ij}+R_{ij})/2}{T_{ij}}.
\]

#### Interpretation

This is the averaged boundary-coverage proxy's share of the local event pool.
It is not a conventional IR ratio based on intron-body coverage, junction
exclusion, transcript abundance or a requirement that the same molecule support
both boundaries.

#### Cryptic-junction overlap

A cryptic donor inside the intron can leave an aligned block spanning the
canonical intron-start boundary; a cryptic acceptor can leave a block spanning
the canonical intron-end boundary. If that block covers all eight bases, the
read contributes to one boundary count while its split portion also contributes
to the cryptic-junction count. A single cryptically spliced molecule can
therefore influence both numerator categories in the same overall analysis.

### Cryptic acceptor use

#### Detection

A junction is labelled as cryptic acceptor use when it retains the annotated
donor boundary but terminates at a genomic coordinate that is not the matching
reference acceptor boundary. Strand determines which genomic end represents the
acceptor.

For a plus-strand transcript, the observed junction shares a reference intron
start but has a novel end. For a minus-strand transcript, it shares a reference
intron end but has a novel start.

#### Raw count and proportion

The raw count is the number of reads supporting the exact cryptic junction:

\[
P^{cryptic\ acceptor}_{ijr} =
\frac{C^{cryptic\ acceptor}_{ijr}}{T_{ij}}.
\]

Each distinct cryptic coordinate is a separate event; nearby sites are not
clustered.

### Cryptic donor use

#### Detection

A junction is labelled as cryptic donor use when it begins at a genomic
coordinate that is not the matching reference donor boundary but retains the
annotated acceptor boundary. Strand again determines which genomic end
represents the donor.

For a plus-strand transcript, the observed junction has a novel start and shares
a reference intron end. For a minus-strand transcript, it shares a reference
intron start and has a novel end.

#### Raw count and proportion

\[
P^{cryptic\ donor}_{ijr} =
\frac{C^{cryptic\ donor}_{ijr}}{T_{ij}}.
\]

Each exact coordinate is quantified independently. No minimum distance from the
canonical site or minimum read count is imposed by Cortar.

### Exon skipping

#### Detection

A junction is labelled as exon skipping when its two boundaries match
boundaries belonging to different reference introns. If the matched intron
numbers are \(a\) and \(b\), with \(a<b\), the skipped exons are labelled
\(a+1\) through \(b\).

For example, a junction connecting the start boundary of intron 2 to the end
boundary of intron 5 is labelled as skipping exons 3–5.

#### Raw count and proportion

The raw count is the count of the exact long junction. The junction is assigned
at least to the two boundary introns. With `ria` enabled, it is also assigned to
every intermediate reference intron completely enclosed by the junction.

For each assigned intron \(i\):

\[
P^{skip}_{ijr} = \frac{C^{skip}_{jr}}{T_{ij}}.
\]

The numerator can be identical in several output rows, while the denominators
and proportions differ by intron. These repeated rows are contextual views of
one junction and must not be treated as independent observations or summed.

### Junctions annotated in alternative transcripts

“Alternative” is an annotation status rather than a separate mutually exclusive
event label. A junction receives this status when its two boundaries exactly
match an intron in another RefSeq transcript of the selected gene. If the same
junction also matches the selected reference transcript, canonical status takes
precedence.

The human-readable event label is still determined relative to the selected
reference transcript. Consequently, an alternative-transcript junction can be
described as cryptic donor use, cryptic acceptor use, exon skipping or an
unannotated junction relative to that reference.

Its raw count is the exact junction count and its proportion remains:

\[
P^{alternative}_{ijr} = \frac{C^{alternative}_{ijr}}{T_{ij}}.
\]

The final output does not identify which alternative transcript or transcripts
contain the junction.

### Junctions with neither selected boundary annotated

A junction with neither boundary matching a selected reference intron is
labelled `unannotated junctions`. Such a junction is not automatically included
merely because it lies inside the gene. It enters an intron's event pool only if
`ria` is enabled and the junction interval completely contains that reference
intron.

Short fully novel junctions that share no selected boundary and contain no
selected intron are therefore absent from the results.

For every intron to which a qualifying junction is assigned:

\[
P^{unannotated}_{ijr} = \frac{C^{unannotated}_{jr}}{T_{ij}}.
\]

### Unknown-strand cryptic junctions

When exactly one boundary matches but the junction strand is unknown, Cortar
labels the event as `cryptic (strand unknown)` rather than donor or acceptor use.
It is counted and normalised in the same way as other junctions. Non-canonical
splice motifs can cause BAM-derived junction strand to be unknown.

## Worked multi-event example

For one sample at one reference intron:

| Evidence category | Raw evidence |
|---|---:|
| Canonical junction | 80 |
| Cryptic acceptor junction | 20 |
| Exon-skipping junction | 5 |
| Intron-start boundary coverage | 12 |
| Intron-end boundary coverage | 8 |

First calculate the IR count:

```text
IR count = (12 + 8) / 2 = 10
```

The event-pool total is:

```text
Total = 80 + 20 + 5 + 10 = 115
```

The reported proportions are:

| Event | Calculation | Proportion |
|---|---:|---:|
| Canonical junction | 80 / 115 | 0.6957 |
| Cryptic acceptor | 20 / 115 | 0.1739 |
| Exon skipping | 5 / 115 | 0.0435 |
| Intron retention | 10 / 115 | 0.0870 |

These proportions describe evidence composition at this intron only. The
exon-skipping junction may have another percentage in the event pool of another
intron that it spans.

## Cross-sample representation

Cortar forms the union of event coordinates observed across all samples. If a
junction is observed in one sample but not another, it is represented with a
raw count and proportion of zero in the latter sample. Zero therefore means “no
supporting alignment was counted”; it does not establish adequate coverage for
confident absence.

In one-versus-many analyses, Cortar calculates the arithmetic mean and sample
standard deviation of control proportions for each event. The test-sample
difference is its proportion minus the control mean. Two-, three- and four-SD
flags indicate whether the absolute difference is strictly greater than the
corresponding multiple of the control SD. These are descriptive rules, not
hypothesis tests.

In research mode, the arithmetic mean and sample standard deviation are
calculated across all samples. Samples are equally weighted regardless of read
depth; the result is the mean of sample proportions, not a pooled-count
proportion.

## Approximate frame annotation

Splice-junction rows also receive a frame-conservation field. Cortar compares
the displacement between observed and annotated genomic boundaries modulo
three. IR rows receive no frame assessment.

This is a coordinate-based heuristic, not a coding-sequence consequence model.
It does not account reliably for transcript-specific coding start, coding phase,
untranslated regions or premature termination. It should not be used as evidence
that an RNA or protein product remains in frame.

## Minimum information to report

A scientific report using Cortar should state:

- Cortar version or commit;
- genome assembly and reference annotation;
- how the reference transcript was selected;
- BAM or STAR/IR input route;
- paired-end and strandedness settings;
- whether `ria` was enabled;
- control-selection mode and control count;
- coverage threshold used in default mode;
- that reported values are within-intron evidence proportions rather than
  conventional PSI;
- that IR is an averaged two-boundary proxy and can include boundary-spanning
  blocks from cryptically split reads.

## Implementation audit provenance

The scientific rules above were traced to the extraction, annotation,
quantification and comparison stages at commit `9f693bd`. Package metadata at
that commit reports version 1.0.0, while the repository tag is 0.1.12.
The corresponding source audit is summarised in the
[developer guide](developer-codebase-guide.md), and outstanding validation work
is listed in the [synthetic BAM validation plan](synthetic-bam-test-plan.md).
