"""Portable locations for rebuilding displays from saved measurements."""
import csv
import hashlib
import os
from pathlib import Path
PROJECT = Path(__file__).resolve().parents[2]
ROOT = PROJECT / 'evidence'
RUN = Path(os.environ['ENTEROBACTER_REBUILD_OUTPUT']).resolve()
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rows(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))
