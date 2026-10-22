"""Select a checked local or relocated development package, without new policy."""
from __future__ import annotations

import copy
import json
from pathlib import Path

from .evidence import sha256
from .reference_package import validate_package


def load_runtime_package(path: Path, *, expected_sha256: str):
    path = Path(path).absolute()
    if not path.is_file() or sha256(path) != expected_sha256:
        raise ValueError("reference package SHA256 mismatch")
    data = json.loads(path.read_text())
    from .portable_reference import SCHEMA as PORTABLE_SCHEMA, load_portable
    if data.get("schema_version") == PORTABLE_SCHEMA:
        return load_portable(path, expected_sha256=expected_sha256)
    return LocalReference(path, expected_sha256, data)


class LocalReference:
    def __init__(self, path: Path, digest: str, data: dict):
        self.validation = validate_package(data)
        if sha256(path) != digest:
            raise ValueError("local reference package changed during validation")
        self._manifest = copy.deepcopy(data)
        self.package_path = str(path)
        self.package_sha256 = digest

    @property
    def manifest(self) -> dict:
        return copy.deepcopy(self._manifest)

    def engine_record(self, component_id: str, binary_path: Path | None = None) -> dict:
        record = self._manifest["recipe"]["components"][component_id]["skani_binary"]
        path = Path(binary_path).absolute() if binary_path else Path(record["path"])
        if not path.is_file() or sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise ValueError("skani executable differs from package pin")
        return {**record, "path": str(path)}
