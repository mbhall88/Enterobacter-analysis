import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

from restore_packages import restore


class RestoreTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def fixture(self, name='asset', content=b'ACGT\n', kind=tarfile.REGTYPE):
        archive = self.root / 'input.tar.gz'
        with tarfile.open(archive, 'w:gz') as tf:
            member = tarfile.TarInfo(name)
            member.type = kind
            if kind == tarfile.REGTYPE:
                member.size = len(content)
                tf.addfile(member, io.BytesIO(content))
            else:
                member.linkname = 'outside'
                tf.addfile(member)
        manifest = {
            'schema_version': 1,
            'archives': {'input.tar.gz': {'bytes': archive.stat().st_size,
                                        'sha256': hashlib.sha256(archive.read_bytes()).hexdigest()}},
            'targets': [{'target': 'package/asset', 'archive': 'input.tar.gz', 'member': name,
                         'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}],
        }
        return manifest

    def run_restore(self, manifest):
        path = self.root / 'manifest.json'
        path.write_text(json.dumps(manifest))
        return restore(path, self.root, self.root / 'output')

    def test_exact_restoration_and_receipt(self):
        receipt = self.run_restore(self.fixture())
        self.assertEqual((self.root / 'output/package/asset').read_bytes(), b'ACGT\n')
        self.assertEqual(receipt['restored_files'], 1)
        self.assertFalse(receipt['classification_executed'])
        self.assertFalse((self.root / 'output/INCOMPLETE').exists())

    def test_existing_output_is_preserved(self):
        manifest = self.fixture()
        out = self.root / 'output'
        out.mkdir()
        (out / 'sentinel').write_text('keep')
        with self.assertRaisesRegex(ValueError, 'new directory'):
            self.run_restore(manifest)
        self.assertEqual((out / 'sentinel').read_text(), 'keep')

    def test_target_traversal_rejected_before_writing(self):
        manifest = self.fixture()
        manifest['targets'][0]['target'] = '../escape'
        with self.assertRaisesRegex(ValueError, 'unsafe relative'):
            self.run_restore(manifest)
        self.assertFalse((self.root / 'output').exists())

    def test_tar_traversal_rejected(self):
        manifest = self.fixture('../escape')
        manifest['targets'][0]['member'] = 'other'
        with self.assertRaisesRegex(ValueError, 'unsafe relative'):
            self.run_restore(manifest)
        self.assertFalse((self.root / 'escape').exists())

    def test_absolute_tar_header_rejected(self):
        manifest = self.fixture('/outside')
        manifest['targets'][0]['member'] = 'other'
        with self.assertRaisesRegex(ValueError, 'unsafe relative'):
            self.run_restore(manifest)

    def test_symlink_member_rejected(self):
        with self.assertRaisesRegex(ValueError, 'links and special'):
            self.run_restore(self.fixture(kind=tarfile.SYMTYPE))

    def test_corrupt_archive_rejected_before_writing(self):
        manifest = self.fixture()
        (self.root / 'input.tar.gz').write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'archive checksum/size'):
            self.run_restore(manifest)
        self.assertFalse((self.root / 'output').exists())

    def test_member_checksum_and_incomplete_marker(self):
        manifest = self.fixture()
        manifest['targets'][0]['sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'member checksum'):
            self.run_restore(manifest)
        self.assertTrue((self.root / 'output/INCOMPLETE').exists())

    def test_missing_member_rejected(self):
        manifest = self.fixture()
        manifest['targets'][0]['member'] = 'absent'
        with self.assertRaisesRegex(ValueError, 'required member missing'):
            self.run_restore(manifest)

    def test_duplicate_target_rejected(self):
        manifest = self.fixture()
        manifest['targets'].append(manifest['targets'][0].copy())
        with self.assertRaisesRegex(ValueError, 'duplicate restoration'):
            self.run_restore(manifest)

    def test_original_manifest_identity_checked(self):
        manifest = self.fixture()
        manifest['original_bound_manifests'] = {'package/asset': '0' * 64}
        with self.assertRaisesRegex(ValueError, 'original bound manifest'):
            self.run_restore(manifest)


if __name__ == '__main__':
    unittest.main()
