# Source-specific data terms

This repository is an integration of separately licensed third-party sources. The publisher of this integration does not claim ownership of the original data or grant permissions beyond the original terms.

## OBELiX and bundled secondary transcriptions

CC-BY-4.0. Preserve attribution to the original authors and cite DOI 10.1039/D5DD00441A and https://github.com/NRC-Mila/OBELiX. The original license is retained at `data/raw/obelix/LICENSE.txt`. The pinned upstream commit is `4eaac889809dd489e3bae468c45d5c0e02b9147a`. Released CIFs include perturbed coordinates and do not grant access to the ICSD database.

The Laskowski CSV in this release is a damaged secondary transcription distributed in that repository. Its provenance and uncertainty are retained; this is not a verified reissue of the complete original publisher SI. Cite the underlying Laskowski paper, DOI 10.1039/D2EE03499A. No additional rights to the original SI are asserted.

## Liverpool Ionics Dataset

The official University of Liverpool page states:

> Where a user wishes to use this Database for the purposes of any extraction or re-utilisation for commercial use, then such use may only be granted under licence by the University. To request such consent and the licence terms under which such consent is granted, please contact impactscience@liverpool.ac.uk and ipadmin@liverpool.ac.uk. Any use of the Database, whether for the purpose of illustration for teaching or academic research or under a commercial licence, should include an appropriate citation to acknowledge the source of the Database as the University of Liverpool.

Source: https://pcwww.liverpool.ac.uk/~msd30/lmds/LiIonDatabase.html (snapshot retained at `data/raw/unified_sources/liion_original.html`). Cite DOI 10.1038/s41524-022-00951-z. The integrated tables include Liverpool-derived rows and are not offered as unrestricted commercial-use data. Contact the University for commercial licensing questions.

## McHaffie / Caltech

The CaltechDATA record specifies CC0-1.0 for the experimental conductivity dataset: https://data.caltech.edu/records/23mvv-6gk43 (dataset DOI 10.22002/23mvv-6gk43). The saved metadata is at `data/raw/new_sources/caltech_record_metadata.json`. Cite the associated paper DOI 10.1039/D5DD00052A.

## DDSE / DigBat

No bulk conductivity dataset from this service is included. The saved public access-status text is evidence of availability, not a license for the underlying database.

## Integrated files and scripts

Per-row source attribution, original terms, raw values, and quality flags remain in the tables and source manifest. No single blanket open-data license is applied to the combined dataset. The accompanying original cleaning scripts have no separately selected software license in this release; do not infer one from an upstream data license.
