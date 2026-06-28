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
