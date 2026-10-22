#!/usr/bin/env python3
"""Restore hash-bound historic packages from explicit archives, without running analyses."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile


def safe_relative(value):
    if not isinstance(value, str) or not value or '\\' in value:
        raise ValueError(f'invalid relative path: {value!r}')
    path = PurePosixPath(value)
    if path.is_absolute() or any(p in ('', '.', '..') for p in value.split('/')):
        raise ValueError(f'unsafe relative path: {value!r}')
    return path


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def validate_manifest(manifest):
    if manifest.get('schema_version') != 1:
        raise ValueError('unsupported manifest schema')
    archives = manifest['archives']
    targets = manifest['targets']
    seen = set()
    for name, entry in archives.items():
        if len(safe_relative(name).parts) != 1:
            raise ValueError('archive names must be basenames')
        if not re.fullmatch('[0-9a-f]{64}', entry['sha256']):
            raise ValueError('invalid archive hash')
        if type(entry['bytes']) is not int or entry['bytes'] < 0:
            raise ValueError('invalid archive size')
    for row in targets:
        safe_relative(row['target'])
        safe_relative(row['member'])
        if row['target'] in seen:
            raise ValueError('duplicate restoration target')
        seen.add(row['target'])
        if row['archive'] not in archives:
            raise ValueError('unknown source archive')
        if not re.fullmatch('[0-9a-f]{64}', row['sha256']):
            raise ValueError('invalid member hash')
        if type(row['bytes']) is not int or row['bytes'] < 0:
            raise ValueError('invalid member size')
    for path in seen:
        if any(str(parent) in seen for parent in PurePosixPath(path).parents):
            raise ValueError('file/directory target collision')
    return archives, targets


def restore(manifest_path, archive_dir, output):
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    archives, targets = validate_manifest(manifest)
    archive_dir, output = Path(archive_dir), Path(output)
    if output.exists() or output.is_symlink():
        raise ValueError('output must be a new directory; existing paths are never overwritten')
    for parent in [output.parent, *output.parents]:
        if parent.is_symlink():
            raise ValueError('output ancestors must not be symbolic links')
    if not output.parent.is_dir():
        raise ValueError('output parent must already exist')
    # Hash every whole archive before creating output or reading member contents.
    for name, expected in archives.items():
        source = archive_dir / name
        if source.is_symlink() or not source.is_file():
            raise ValueError(f'archive is missing or is a symlink: {name}')
        if source.stat().st_size != expected['bytes'] or digest(source) != expected['sha256']:
            raise ValueError(f'archive checksum/size mismatch: {name}')
    output.mkdir()
    (output / 'INCOMPLETE').write_text('Restoration has not completed; do not use these packages.\n')
    completed = []
    for name in archives:
        wanted = {}
        for row in targets:
            if row['archive'] == name:
                wanted.setdefault(row['member'], []).append(row)
        seen = set()
        with tarfile.open(archive_dir / name, 'r|gz') as tar:
            for member in tar:
                # Inspect every header, including unused public index assets. Never extractall.
                safe_relative(member.name.rstrip('/') if member.isdir() else member.name)
                if member.name in seen:
                    raise ValueError('duplicate archive member')
                seen.add(member.name)
                if member.isdir():
                    continue
                if not member.isfile():
                    raise ValueError('archive links and special members are refused')
                rows = wanted.get(member.name)
                if not rows:
                    continue
                if any(member.size != row['bytes'] for row in rows):
                    raise ValueError(f'member size mismatch: {member.name}')
                handles = []
                try:
                    for row in rows:
                        destination = output / row['target']
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        handles.append(destination.open('xb'))
                    state = hashlib.sha256()
                    source = tar.extractfile(member)
                    for block in iter(lambda: source.read(1024 * 1024), b''):
                        state.update(block)
                        for handle in handles:
                            handle.write(block)
                    observed = state.hexdigest()
                    if any(observed != row['sha256'] for row in rows):
                        raise ValueError(f'member checksum mismatch: {member.name}')
                    completed.extend(row['target'] for row in rows)
                finally:
                    for handle in handles:
                        handle.close()
        if any(member not in seen for member in wanted):
            raise ValueError(f'required member missing from archive: {name}')
    if len(completed) != len(targets):
        raise ValueError('not every target was restored')
    # Verify actual destination bytes independently of the streaming copy.
    for row in targets:
        if digest(output / row['target']) != row['sha256']:
            raise ValueError('restored-file readback mismatch')
    for relative, expected in manifest.get('original_bound_manifests', {}).items():
        safe_relative(relative)
        if digest(output / relative) != expected:
            raise ValueError('original bound manifest differs after restoration')
    receipt = {
        'manifest_sha256': digest(manifest_path),
        'restored_files': len(targets),
        'restored_bytes': sum(row['bytes'] for row in targets),
        'all_archive_and_target_hashes_verified': True,
        'original_bound_manifests_verified': len(manifest.get('original_bound_manifests', {})),
        'classification_executed': False,
        'portable_full_rerun_verified': False,
    }
    (output / 'restoration-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (output / 'INCOMPLETE').unlink()
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=Path(__file__).with_name('restoration-manifest.json'))
    parser.add_argument('--archive-dir', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    try:
        receipt = restore(args.manifest, args.archive_dir, args.output)
    except (OSError, ValueError, KeyError, TypeError, tarfile.TarError) as error:
        parser.exit(1, f'Restoration failed: {error}\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
