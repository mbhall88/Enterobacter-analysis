"""End-to-end table reproduction and provenance/refusal tests; no external data."""
import json,tempfile,unittest
from pathlib import Path
import rebuild_tables as r
BASE=Path(__file__).resolve().parent
class ReproductionTests(unittest.TestCase):
    def test_all_23_tables_match_frozen_hashes(self):
        with tempfile.TemporaryDirectory(dir=BASE.parent) as tmp:
            out=Path(tmp)/'rebuilt';r.execute(BASE,out)
            result=json.loads((out/'reconstruction-checks.json').read_text())
            self.assertEqual(len(result['tables']),23)
            self.assertTrue(all(x['exact'] for x in result['tables']))
            self.assertEqual((out/'TableS13.tsv').read_text(),'coded_id\tdataset\tsetting\tcategory\tforeign_challenge_references\n')
    def test_compressed_input_corruption_refused(self):
        item=json.loads((BASE/'input-inventory.json').read_text())[0];item['export_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'Compressed input checksum mismatch'):r.Input(item,BASE).read_bytes()
    def test_decoded_input_corruption_refused(self):
        item=json.loads((BASE/'input-inventory.json').read_text())[0];item['decoded_export_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'Decoded input checksum mismatch'):r.Input(item,BASE).read_bytes()
    def test_existing_output_preserved(self):
        with tempfile.TemporaryDirectory(dir=BASE.parent) as tmp:
            out=Path(tmp);sentinel=out/'keep.txt';sentinel.write_text('unchanged')
            with self.assertRaises(FileExistsError):r.execute(BASE,out)
            self.assertEqual(sentinel.read_text(),'unchanged')
if __name__=='__main__':unittest.main()
