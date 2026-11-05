# Known issues

- Laskowski: the 1,345 imported nonempty conductivity cells are from a damaged secondary transcription, not a primary-verified recovery of the reported 1,346 entries. These rows are excluded from the default primary-clean subset.
- DDSE: no bulk data obtained. The current public interface restricts bulk copy to collaborators.
- Source overlap: 3,335 source records are not 3,335 independent samples; use the cross-source candidate links and leakage groups.
- Room-temperature candidates: numeric temperature, bulk/total definition, synthesis, and direct-versus-extrapolated status are often incomplete.
- Released CIFs include coordinate perturbations. Structural matching is not independently verified for every record.
- Additional structural-field audit: OBELiX `r5y` has main composition `Li2CdCl4`, but `True Composition` is `Li1.6Cd8.0Cl32.0`. Their largest elemental atomic-fraction discrepancy is about0.247. The source note says “composition corrected”, but the two fields still disagree. The original fields and conductivity were preserved. Do not use this row to calculate cell-normalized atom volume or lithium density without resolving the discrepancy; see `data/processed/unified/structural_metadata_review.csv`. Its main composition-conductivity record is not silently deleted.

No new transport law or physical mechanism has been established by this dataset integration.

## Original close-match flag parsing correction

The original table uses `Yes`; the previous unified adapter only recognized `true`. This has been corrected:61 CIF records in all data and55 in the primary-clean subset are approximate matches. Labels and raw source files are unchanged. See `reports/effective_dataset_coverage.md` for composition deduplication, conservative structure evidence, processing coverage and external-test limitations.
