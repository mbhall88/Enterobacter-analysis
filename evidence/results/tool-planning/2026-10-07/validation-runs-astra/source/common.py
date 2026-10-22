"""Execution-only utilities for frozen validation design revision 6."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path('${PROJECT_ROOT}')
RUN = ROOT / 'results/tool-planning/2026-10-07/validation-runs-astra'
DESIGN = ROOT / 'results/tool-planning/2026-10-06/validation-design-claude'
ART = DESIGN / 'outputs/design-artifacts-attempt08'
FREEZE = DESIGN / 'outputs/freeze-hashes-attempt01'
DOWNLOAD = ROOT / 'results/ncbi/2026-10-07/validation-reserve-claude/manifest.tsv'
EXPERT = ROOT / 'results/tool-planning/2026-10-07/expert-labelling-claude/seal-attempt01/expert-truth.tsv'
PYTHON = ROOT / '.work/envs/verticall/bin/python'
PINS = {
    'package': (ROOT / 'results/tool-planning/2026-10-06/portable-reference-package-astra/relocated-package/portable.json',
                '12aa0705b774de8ea262869fdc4c85d513c3aa97d39c674f08e85d1c9055b88f'),
    'rank': (ROOT / 'results/tool-planning/2026-10-06/frozen-eligibility-claude/inputs/rank-rule.json',
             '5849c5db6649678963215ca9363a5dfae2ef2b0282d96946aa843814fafa4188'),
    'config': (ROOT / 'results/tool-planning/2026-10-06/near-boundary-claude/inputs/frozen-eligibility-config-2.1.json',
               '3188518eb7bdaad287e7f335ffb18cadfbb41b46308f018f5b69ff33cab7c626'),
    'mlst': (ROOT / 'results/tool-planning/2026-10-06/mlst-runner-integration-astra/inputs/mlst-config.json',
             '1fce9a33cf423967e0c089830ee3e6549930b5d45761ab95ac662ccd1430a937'),
    'skani': (Path('/sw/local/rocky9/noarch/rcc/software/gtdb-tk/gtdbtk_2.7.2/bin/skani'),
              'b3000043feb4ce8a5200b3716bae155f34f73e2a123b8b055d574a7897a14151'),
}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, fields=None):
    rows = list(rows)
    with Path(path).open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def json_write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_object():
    """Verify actual package sources, not just the frozen table listing them."""
    if sha(FREEZE / 'freeze-hashes.tsv') != '60e10c5fd87fc41a1a0ab2ed2989ea81b1cfa12ac5ddfc5cdbbf6f1c8fb799aa':
        raise ValueError('posted freeze inventory changed')
    sources = read(FREEZE / 'package-source.tsv')
    digest = hashlib.sha256(''.join(f"{r['sha256']}  {r['path']}\n"
                            for r in sorted(sources, key=lambda r: r['path'])).encode()).hexdigest()
    if digest != '52ec296c27795b35399712183210d96ba888eabac4edcf8384838ea866c30089':
        raise ValueError('posted package source inventory changed')
    actual = {str(p.relative_to(ROOT)) for p in (ROOT / 'src/enterobacter_classifier').rglob('*')
              if p.is_file() and '__pycache__' not in p.parts}
    expected = {r['path'] for r in sources}
    if actual != expected or len(sources) != 19:
        raise ValueError('frozen package source inventory changed')
    # The full 40-file preflight also checked the large, labelling-only files.
    # Recheck all design/code/artifact/object inputs used by this execution worker.
    small = [r for r in read(FREEZE / 'freeze-hashes.tsv')
             if r['group'] in {'design', 'artifacts', 'code', 'pins'}]
    checks = [(ROOT / r['path'], r['sha256']) for r in sources + small] + list(PINS.values())
    for path, digest in checks:
        if sha(path) != digest:
            raise ValueError(f'frozen object changed: {path}')
    return {'package_source_files': len(sources), 'object_pins': len(PINS), 'changed': 0}


def command(query, output):
    """Exactly the frozen section 1 analyze command: no extra evidence or QC."""
    return [str(PYTHON), '-m', 'enterobacter_classifier', 'analyze',
            '--reference-package', str(PINS['package'][0]), '--expected-package-sha256', PINS['package'][1],
            '--experimental-rank-rule', str(PINS['rank'][0]),
            '--experimental-eligibility-config', str(PINS['config'][0]),
            '--expected-eligibility-config-sha256', PINS['config'][1], '--experimental-genomovars',
            '--mlst-config', str(PINS['mlst'][0]), '--expected-mlst-config-sha256', PINS['mlst'][1],
            '--skani-binary', str(PINS['skani'][0]), '--threads', '2', '--query', str(query), '--output', str(output)]


def batch(rows, index, size):
    if size < 1 or index < 0 or index * size >= len(rows):
        raise ValueError('array index outside manifest')
    return rows[index * size:(index + 1) * size]
