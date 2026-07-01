# HalC8 Detection Pipeline v2.2

Public software-only release of the HalC8 v2.2 domain-scoring script.

## Contents

- `halc8_score_domains.py` — standalone HalC8-family domain scoring script.
- `environment.yml` — conda environment specification.
- `requirements.txt` — Python package requirements.
- `LICENSE` — software license.
- `CITATION.cff` — citation metadata for this software repository.

## Installation

Using conda:

    conda env create -f environment.yml
    conda activate halc8-pipeline

Using pip:

    python3 -m pip install -r requirements.txt

## Basic usage

    python3 halc8_score_domains.py --input input_domains.csv --output scored_domains.csv --input-type domain --ss-mode none

Input CSV files require at least two columns: `accession` and `sequence`.

## License

See `LICENSE`.

## Data provenance

The 2021 prior-label table used for calibration-frame provenance is archived at Zenodo: https://doi.org/10.5281/zenodo.21116990. This table is used as a historical prior-label framework, not as an experimentally validated benchmark or independent external validation set.

