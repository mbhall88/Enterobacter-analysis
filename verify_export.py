#!/usr/bin/env python3
"""Verify the exact publication payload, excluding this manifest and .git metadata."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys


def verify(root):
    root = Path(root)
    rows = json.loads((root / 'MANIFEST.json').read_text())['files']
    expected = {}
    for row in rows:
        name = row['path']
        parts = PurePosixPath(name).parts
        if not parts or name.startswith('/') or '..' in parts or '\\' in name or name in expected:
            raise ValueError('Unsafe or duplicate manifest path: ' + name)
        expected[name] = row
    actual = set()
    for p in root.rglob('*'):
        rel = p.relative_to(root)
        if rel.parts[0] == '.git':
            continue
        if p.is_symlink():
            raise ValueError('Unexpected symlink: ' + str(rel))
        if p.is_dir() or str(rel) == 'MANIFEST.json':
            continue
        if not p.is_file():
            raise ValueError('Nonregular file: ' + str(rel))
        name = rel.as_posix()
        actual.add(name)
        row = expected.get(name)
        if row is None:
            raise ValueError('Unexpected file: ' + name)
        if p.stat().st_size != row['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('Changed file: ' + name)
    if actual != set(expected):
        raise ValueError('Missing files: ' + ', '.join(sorted(set(expected) - actual)))
    return {'files_verified': len(actual), 'all_hashes_match': True}


if __name__ == '__main__':
    try:
        print(json.dumps(verify(Path(__file__).resolve().parent)))
    except (ValueError, KeyError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
