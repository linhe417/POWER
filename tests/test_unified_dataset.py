"""Scientific integrity checks for source integration and label interpretation."""
import unittest, json, math
from pathlib import Path
import pandas as pd
from src.unify_datasets import DEST, ROOT, parse_sigma, parse_formula, formula_syntax, short_hash

class UnifiedDatasetIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=pd.read_csv(DEST/'measurement_records.csv',keep_default_na=False)

    def test_censoring_is_not_an_exact_regression_target(self):
        for text,relation in [('<1E-10','<'),('≤ 1e-8','<='),('>2e-5','>'),('~1e-3','~')]:
            val,rel,log,_=parse_sigma(text)
            self.assertGreater(val,0);self.assertEqual(rel,relation);self.assertTrue(math.isnan(log))
        c=self.d[self.d.sigma_relation.ne('=')]
        self.assertTrue(c.log10_sigma_S_cm.eq('').all())
        self.assertFalse(c.eligible_primary_clean_numeric.any())

    def test_formula_parser_rejects_detached_pdf_subscripts(self):
        self.assertEqual(parse_formula('Li La Zr O 7 3 2 12')[2],'ambiguous_detached_subscripts')
        self.assertFalse(parse_formula('Li6La3XxO12')[0])
        a=parse_formula('"(Li2S)0.8 (P2S5)0.2"')
        self.assertEqual(a[2],'parsed')
        self.assertAlmostEqual(a[1]['Li'],1.6)
        self.assertEqual(parse_formula('Li4SnS3Se')[0],parse_formula('Li8Sn2S6Se2')[0])
        self.assertNotEqual(parse_formula('Li4SnS3Se')[0],parse_formula('Li4SnSSe3')[0])

    def test_source_count_conservation_and_unique_ids(self):
        self.assertEqual(self.d.dataset.value_counts().to_dict(),{'laskowski':1345,'liion':820,'obelix':599,'mchaffie':571})
        self.assertEqual(len(self.d),3335)
        self.assertTrue(self.d.record_id.is_unique)
        self.assertFalse(self.d.source_payload_json.eq('').any())
        for value in self.d.source_payload_json:json.loads(value)

    def test_temperature_units_and_missingness(self):
        li=self.d[self.d.dataset.eq('liion')]
        self.assertTrue(((pd.to_numeric(li.temperature_raw)+273.15-pd.to_numeric(li.measurement_temperature_K)).abs()<1e-8).all())
        other=self.d[~self.d.dataset.eq('liion')]
        self.assertTrue(other.measurement_temperature_K.eq('').all())
        ex=self.d[self.d.measurement_method_status.eq('rt_extrapolated_from_higher_T')]
        self.assertEqual(len(ex),135)
        self.assertFalse(ex.eligible_reported_rt_candidate.any())

    def test_official_liion_correction_is_used(self):
        r=self.d.set_index('record_id').loc['liion:703']
        self.assertEqual(r.formula_raw,'Li6NBr3')
        self.assertEqual(r.doi_normalized,'10.1002/anie.198000741')
        self.assertAlmostEqual(float(r.sigma_value_S_cm),3.1e-7)

    def test_unverified_transcriptions_cannot_enter_primary_clean_view(self):
        la=self.d[self.d.dataset.eq('laskowski')]
        self.assertTrue(la.provenance_status.eq('secondary_damaged_transcription_unverified').all())
        self.assertFalse(la.eligible_primary_clean_numeric.any())
        self.assertEqual(len(pd.read_csv(DEST/'laskowski_all_transcription_fragments.csv')),1928)

    def test_excel_preserves_relation_strings_without_formulas(self):
        import sys
        sys.path.insert(0,str(ROOT/'experiments/vendor'))
        import openpyxl
        w=openpyxl.load_workbook(DEST/'unified_electrolytes.xlsx',read_only=True,data_only=False)
        rows=w['all_source_records'].iter_rows()
        fields=[c.value for c in next(rows)]
        idx=fields.index('sigma_relation')
        self.assertEqual(fields,list(self.d.columns))
        count=0
        for cells in rows:
            for cell in cells:self.assertNotEqual(cell.data_type,'f')
            self.assertEqual(cells[idx].value,self.d.iloc[count].sigma_relation)
            count+=1
        self.assertEqual(count,len(self.d));w.close()

    def test_structures_are_hash_verified_and_not_fabricated_for_icsd_ids(self):
        import hashlib
        m=pd.read_csv(DEST/'structure_manifest.csv')
        self.assertEqual(len(m),321)
        for r in m.itertuples():
            self.assertEqual(hashlib.sha256((ROOT/r.cif_path).read_bytes()).hexdigest(),r.cif_sha256)
        self.assertTrue(self.d[self.d.dataset.eq('mchaffie')].cif_path.eq('').all())

    def test_cross_source_links_do_not_assert_sample_identity(self):
        l=pd.read_csv(DEST/'cross_source_links.csv')
        self.assertFalse(l.same_physical_sample_confirmed.any())
        lookup=self.d.set_index('record_id')
        for r in l.itertuples():
            a,b=lookup.loc[r.record_id_a],lookup.loc[r.record_id_b]
            self.assertEqual(a.composition_key,b.composition_key)
            self.assertIn(r.shared_doi,a.doi_normalized.split(';'))
            self.assertIn(r.shared_doi,b.doi_normalized.split(';'))

    def test_no_composition_or_publication_crosses_leakage_group(self):
        groups={}
        for r in self.d.itertuples():
            keys=(['c:'+r.composition_key] if r.composition_key else [])+['d:'+d for d in r.doi_normalized.split(';') if d]
            for k in keys:
                if k in groups:self.assertEqual(groups[k],r.leakage_group_id)
                groups[k]=r.leakage_group_id

if __name__=='__main__':unittest.main()
