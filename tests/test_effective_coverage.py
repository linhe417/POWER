import unittest,json
import pandas as pd
from src.common import ROOT

class EffectiveCoverageIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=ROOT/'data/processed/effective_coverage'
        cls.s=json.loads((cls.p/'effective_coverage_summary.json').read_text())
        cls.d=pd.read_csv(ROOT/'data/processed/unified/measurement_records.csv',keep_default_na=False)

    def test_composition_counts_are_not_asserted_to_be_independent_samples(self):
        clean=self.d[self.d.eligible_primary_clean_numeric]
        rt=self.d[self.d.eligible_reported_rt_candidate]
        self.assertEqual(clean.composition_id.nunique(),925)
        self.assertEqual(rt.composition_id.nunique(),799)
        self.assertIsNone(self.s['independent_physical_sample_count'])
        self.assertEqual(len(pd.read_csv(self.p/'primary_clean_composition_index.csv')),925)

    def test_conservative_structure_candidates_remain_unverified(self):
        e=pd.read_csv(self.p/'clean_structure_evidence.csv')
        self.assertEqual(len(e),260)
        self.assertEqual(e.composition_id.nunique(),233)
        q=e[e.tier.eq('unperturbed_not_close_candidate_unverified')]
        self.assertEqual(len(q),23);self.assertEqual(q.composition_id.nunique(),23)
        self.assertFalse(e.same_measurement_sample_primary_verified.any())

    def test_quantitative_processing_is_partial_and_not_invented_for_core(self):
        p=pd.read_csv(self.p/'partial_quantitative_processing_inventory.csv')
        self.assertEqual(len(p),33)
        self.assertEqual(int(p.in_github_core_table.sum()),1)
        self.assertEqual(int(p.metadata_consistent.sum()),29)
        self.assertFalse(p.fully_standardized_processing_protocol.any())
        self.assertEqual(self.s['processing']['core_structured_synthesis_T_and_time_records'],0)

    def test_new_publication_does_not_imply_new_composition(self):
        a=pd.read_csv(self.p/'supplementary_composition_overlap.csv',keep_default_na=False)
        self.assertIn('mchaffie:0045',a.record_id.tolist())
        primary=a[a.eligible_primary_clean_numeric]
        self.assertFalse(primary.empty)
        self.assertEqual(self.s['testing']['supplementary_records_from_new_composition'],0)
        self.assertEqual(self.s['testing']['locked_untouched_external_samples'],0)

if __name__=='__main__':unittest.main()
