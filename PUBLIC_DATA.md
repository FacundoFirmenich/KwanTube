# KwanTube Public Data Repository

## Overview

The primary computational artefacts and code are hosted on GitHub:
**https://github.com/FacundoFirmenich/KwanTube**

However, the complete companion public-data repository (~1.52 GB) cannot be
included in the GitHub repository due to file size constraints.
It is publicly available, with free access, at the following link:

## 🔗 Access

**Google Drive (free, public, no login required):**

> **https://drive.google.com/drive/folders/1DmXBlcdwP7gY-k56RKdNGRvGSRwd5ZzU?usp=sharing**

Folder name: `KwanTube_public_data_repo`

---

## Contents

The repository is a curated collection of public-domain scientific data
used to construct the L4 empirical validation layer of the KwanTube pipeline.
It contains:

### Structural Data (RCSB Protein Data Bank)
- Full PDB/mmCIF records for tubulin structures (e.g., 1JFF, 6DPU, 8X9P, 4H6U, and others)
- Used by: `src/qmc_mt/pdb_tubulin_analysis.py` → `raw_json/metrics/pdb_tubulin_analysis.json`

### Compound Information (PubChem)
- JSON records for tubulin inhibitors and modulators:
  Epothilone, Methotrexate, Hemiasterlin, Vinorelbine, Colchicine, Taxol, and others
- Used by: `src/scripts/public_data/fetch_public_data.py`

### Literature Search Summaries (OpenAlex, CrossRef, Europe PMC, RCSB)
- `fetch_summary.json`: aggregated search results from 4 scientific databases
- Covers topics: Tau-MT NMR spectroscopy, cancer drug resistance,
  tubulin vibrational modes, UV superradiance in tryptophan networks
- Used by: `src/scripts/public_data/curate_compact.py` → `raw_csv/compact/`

### Spectroscopic & Structural Audit
- 362 tubulin PDB structures with usable B-factor data
- 93 literature studies reporting vibrational modes for microtubules
- Source for the empirical η-proxy audit reported in the manuscript (Section 2.1.1)

---

## Usage Instructions

### Option 1 — Automatic download (pipeline integration)

The pipeline's `fetch_public_data.py` script re-fetches all data
from the original public APIs (RCSB, PubChem, OpenAlex, CrossRef, Europe PMC).
This is the reproducible-from-scratch path.

```bash
python src/scripts/public_data/fetch_public_data.py
```

### Option 2 — Manual download from Drive

1. Open the Drive link above.
2. Download the entire `KwanTube_public_data_repo` folder.
3. Place the contents in:
   ```
   KwanTube/data_downloaded_public_repos/
   ```
4. The pipeline scripts will detect the local data and skip re-fetching.

---

## Data Provenance and Licences

All data in this repository originates from open-access public databases:

| Source | Licence | URL |
|:-------|:--------|:----|
| RCSB Protein Data Bank | CC0 1.0 | https://www.rcsb.org |
| PubChem | Public domain | https://pubchem.ncbi.nlm.nih.gov |
| OpenAlex | CC0 | https://openalex.org |
| CrossRef | CC0 | https://www.crossref.org |
| Europe PMC | CC BY | https://europepmc.org |

No proprietary or restricted data is included. All records are reproducible
by re-running `fetch_public_data.py` against the original APIs.

---

## Code Repository — Zenodo DOI

The KwanTube software package is permanently archived on Zenodo with a citable DOI:

> **DOI: [10.5281/zenodo.19744599](https://doi.org/10.5281/zenodo.19744599)**

To cite this software:

```bibtex
@software{qmcmt_zenodo_2026,
  author    = {Firmenich, Facundo and Firmenich, Pau and Firmenich, León},
  title     = {KwanTube: reproducible open-system and convergence-validation
               framework for tubulin quantum dynamics},
  version   = {3.5.1.1},
  year      = {2026},
  doi       = {10.5281/zenodo.19744599},
  url       = {https://doi.org/10.5281/zenodo.19744599}
}
```

The public data companion repository (this Drive folder) is hosted separately
due to its size (~1.52 GB). It contains exclusively public-domain data
re-fetchable from the original APIs via `fetch_public_data.py`.

---

*For questions, contact: f.firmenich@cedesur.org*
