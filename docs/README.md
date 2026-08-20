# Cortar technical documentation

These documents describe the behaviour of cortar version 1.0.0 at commit
`9f693bd`. They are an implementation audit, not a description of an idealised
splice-analysis method. Where comments, older documentation and executable code
disagree, the executable code is treated as authoritative.

## Start here

- [Developer guide to the codebase](developer-codebase-guide.md): newcomer
  orientation, repository map, call graph, every `R/` module, tests and safe
  change workflow.
- [Manuscript-ready methods brief](cortar-overview.md): scientific calculation
  method, essential limitations and study-specific details to report.
- [Event detection, annotation and quantification](event-detection-and-quantification.md):
  exact counting rules and formulae.
- [One-versus-many/default mode](default-mode.md): test, family and control
  selection and comparison statistics.
- [Panel mode](panel-mode.md): how panel mode differs from default mode.
- [Research mode](research-mode.md): cohort aggregation and outputs.
- [Assumptions, limitations and known defects](assumptions-and-limitations.md):
  biological assumptions, statistical caveats and implementation hazards.
- [Synthetic BAM validation plan](synthetic-bam-test-plan.md): recommended tools,
  fixtures, expected results and test architecture.

## Confidence labels

The documents use three labels:

- **Implemented**: directly established by executable source at the audited
  commit.
- **Inferred**: a consequence of library semantics or control flow, but not
  asserted by a dedicated test.
- **Unverified intent**: a name or comment suggests an intention that the code
  does not reliably implement.

The existing mode documents previously contained unsupported features such as
Ensembl identifier support, panel burden scores, pathway analysis, multiple
testing and batch correction. Those features are not implemented and have been
removed from this documentation.

## Current automated evidence

Running `Rscript -e 'devtools::test(stop_on_failure = FALSE)'` on 18 August 2026
produced 113 passes, with no failures, warnings or skips. The suite includes one
snapshot-style BAM end-to-end test for EMD. It does not yet contain orthogonal
ground-truth BAMs for every event class, strand mode, flag filter, control rule or
zero-coverage edge case. See the [validation plan](synthetic-bam-test-plan.md).
