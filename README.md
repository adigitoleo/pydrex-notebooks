# PyDRex notebooks

Jupyter notebooks for PyDRex and related CPO calculations.

## Dependencies

This is not a python package, so there is no `pyproject.toml`.
Direct dependencies should instead be listed in `requirements.in`.
To compile a `requirements.txt` manifest for the GitHub codespace,
and install the necessary Python dependencies, use the following command:

    pip install --upgrade pip pip-tools && pip-compile --resolver=backtracking && pip-sync

It can take a few minutes.
