"""按截点结局推导的核对：修订选择、分类决胜、确认规则与固定分母。"""
import json
import unittest
from pathlib import Path
from research_data.io import load
from research_data.endpoint import derive

def material(name):
    return load(str(Path('data')/name))[0]

def rec(rid,day,rec_day,ver,targets,nt='stable',nl=False,sid='s1'):
    return {'record_id':rid,'subject_id':sid,'visit_day':day,'recorded_day':rec_day,'version':ver,
            'targets':[{'lesion_id':l,'mm':m} for l,m in targets],'non_target':nt,'new_lesion':nl}

def study(visits,cc=84,kc=86,baseline=None,extra_subjects=True):
    subjects=[{'subject_id':'s1','baseline':baseline or [{'lesion_id':'L1','mm':30},{'lesion_id':'L2','mm':20}]}]
    if extra_subjects:subjects.append({'subject_id':'s2','baseline':[{'lesion_id':'L1','mm':25}]})
    return {'subjects':subjects,'visits':visits,'clinical_cutoff':cc,'knowledge_cutoff':kc}

def s1(out):
    return next(s for s in out['subjects'] if s['subject_id']=='s1')

class MaterialDerivationTests(unittest.TestCase):
    def test_normal_boundary_pr_and_confirmation(self):
        out=derive(material('normal.json'))
        rows=s1(out)['visits']
        self.assertEqual([r['classification'] for r in rows],['PR','PR','NE'])
        self.assertEqual([r['version'] for r in rows],[1,1,1])
        self.assertEqual(rows[0]['pct_of_baseline'],70.0)
        self.assertEqual(rows[0]['nadir_mm'],35)
        self.assertEqual(rows[2]['complete'],False)
        self.assertIn('缺测靶病灶:L2',rows[2]['notes'])
        self.assertEqual(rows[2]['nadir_mm'],32)
        sub=s1(out)
        self.assertEqual(sub['confirmed_response'],{'start_day':28,'confirmed_day':56})
        self.assertIsNone(sub['confirmed_cr'])
        self.assertEqual(sub['best_confirmed'],'PR')
        self.assertEqual(sub['first_candidate_day'],28)
        excluded={(e['visit_day'],e['version']):e for e in out['excluded']}
        self.assertIn('knowledge_cutoff',excluded[(56,2)]['reason'])
        cohort=out['cohort']
        self.assertEqual(cohort['denominator'],2)
        self.assertEqual(cohort['confirmed_response'],{'count':1,'denominator':2,'ratio':'1/2','value':0.5})
        self.assertEqual(cohort['unconfirmed_response']['count'],0)
        self.assertEqual(cohort['no_response']['count'],1)
    def test_confirmation_waiting_across_clinical_cutoff(self):
        out=derive(material('before-confirmation.json'))
        sub=s1(out)
        self.assertEqual([r['visit_day'] for r in sub['visits']],[28])
        self.assertEqual(sub['best_confirmed'],'NONE')
        self.assertIsNone(sub['confirmed_response'])
        self.assertEqual(sub['open_run'],{'start_day':28,'earliest_confirm_day':56})
        reasons={(e['visit_day'],e['version']):e['reason'] for e in out['excluded']}
        self.assertIn('clinical_cutoff',reasons[(56,1)])
        self.assertIn('clinical_cutoff',reasons[(84,1)])
        self.assertIn('knowledge_cutoff',reasons[(56,2)])
        cohort=out['cohort']
        self.assertEqual(cohort['confirmed_response']['count'],0)
        self.assertEqual(cohort['unconfirmed_response'],{'count':1,'denominator':2,'ratio':'1/2','value':0.5})
    def test_late_revision_full_rederivation(self):
        out=derive(material('late-revision.json'))
        rows=s1(out)['visits']
        self.assertEqual([r['classification'] for r in rows],['PR','PD','POST_PD'])
        self.assertEqual(rows[1]['version'],2)
        self.assertEqual(rows[1]['sum_mm'],46)
        self.assertEqual(rows[1]['nadir_mm'],35)
        self.assertEqual(rows[1]['increase_from_nadir_mm'],11)
        self.assertEqual(s1(out)['first_pd_day'],56)
        self.assertEqual(s1(out)['best_confirmed'],'NONE')
        self.assertIsNone(s1(out)['confirmed_response'])
        excluded={(e['visit_day'],e['version']):e for e in out['excluded']}
        self.assertIn('取代',excluded[(56,1)]['reason'])
        normal=derive(material('normal.json'))
        self.assertEqual(s1(normal)['visits'][1]['classification'],'PR')
        self.assertEqual(s1(normal)['best_confirmed'],'PR')
    def test_missing_target_with_new_lesion_triggers_pd(self):
        out=derive(material('missing-target-new-lesion.json'))
        rows=s1(out)['visits']
        self.assertEqual([r['classification'] for r in rows],['PD','POST_PD','POST_PD'])
        self.assertEqual(rows[0]['complete'],False)
        self.assertIsNone(rows[0]['sum_mm'])
        self.assertTrue(any('独立触发PD:新病灶' in n for n in rows[0]['notes']))
        self.assertEqual(s1(out)['first_pd_day'],28)
        self.assertIsNone(s1(out)['first_candidate_day'])
        self.assertEqual(out['cohort']['no_response'],{'count':2,'denominator':2,'ratio':'1','value':1.0})
    def test_pd_wins_over_pr_at_same_visit(self):
        out=derive(material('pd-priority.json'))
        rows=s1(out)['visits']
        self.assertEqual([r['classification'] for r in rows],['PR','PD','POST_PD'])
        self.assertEqual(rows[0]['nadir_mm'],10)
        self.assertEqual(rows[1]['pct_of_baseline'],50.0)
        self.assertEqual(rows[1]['increase_from_nadir_mm'],40)
        self.assertTrue(any('较nadir增大' in n for n in rows[1]['notes']))
        self.assertEqual(s1(out)['best_confirmed'],'NONE')
        self.assertEqual(out['cohort']['unconfirmed_response']['count'],1)
    def test_confirmed_cr_waits_out_short_interval(self):
        out=derive(material('confirmed-cr.json'))
        rows=s1(out)['visits']
        self.assertEqual([r['classification'] for r in rows],['CR','CR','CR'])
        sub=s1(out)
        self.assertEqual(sub['confirmed_response'],{'start_day':28,'confirmed_day':60})
        self.assertEqual(sub['confirmed_cr'],{'start_day':28,'confirmed_day':60})
        self.assertEqual(sub['best_confirmed'],'CR')
        self.assertIsNone(sub['open_run'])
        self.assertEqual(out['cohort']['confirmed_cr'],{'count':1,'denominator':2,'ratio':'1/2','value':0.5})
    def test_invalid_material_rejected(self):
        with self.assertRaises(ValueError):derive(material('invalid-domain.json'))
    def test_derivation_json_serializable(self):
        for p in Path('data').glob('*.json'):
            if p.name.startswith('invalid-'):continue
            with self.subTest(file=p.name):
                json.dumps(derive(material(p.name)),ensure_ascii=False,allow_nan=False)

class RuleTests(unittest.TestCase):
    def test_pd_boundary_equality(self):
        out=derive(study([rec('a',10,11,1,[('L1',30)])],baseline=[{'lesion_id':'L1','mm':25}]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'PD')
        out=derive(study([rec('a',10,11,1,[('L1',29.9)])],baseline=[{'lesion_id':'L1','mm':25}]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'SD')
        out=derive(study([rec('a',10,11,1,[('L1',24)])],baseline=[{'lesion_id':'L1','mm':20}]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'SD')
        out=derive(study([rec('a',10,11,1,[('L1',105)])],baseline=[{'lesion_id':'L1','mm':100}]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'SD')
    def test_cr_beats_pr_and_requires_non_target_absent(self):
        out=derive(study([rec('a',10,11,1,[('L1',0),('L2',0)],nt='absent')]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'CR')
        out=derive(study([rec('a',10,11,1,[('L1',0),('L2',0)],nt='stable')]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'PR')
    def test_independent_pd_triggers_on_complete_visit(self):
        out=derive(study([rec('a',10,11,1,[('L1',0),('L2',0)],nt='absent',nl=True)]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'PD')
        out=derive(study([rec('a',10,11,1,[('L1',5),('L2',5)],nt='progression')]))
        self.assertEqual(s1(out)['visits'][0]['classification'],'PD')
    def test_unknown_non_target_is_ne(self):
        out=derive(study([rec('a',10,11,1,[('L1',5),('L2',5)],nt='unknown')]))
        row=s1(out)['visits'][0]
        self.assertEqual((row['classification'],row['complete']),('NE',False))
    def test_ne_and_sd_break_confirmation_run(self):
        visits=[rec('a',28,29,1,[('L1',20),('L2',15)]),
                rec('b',56,57,1,[('L1',18),('L2',None)]),
                rec('c',84,85,1,[('L1',20),('L2',15)]),
                rec('d',112,113,1,[('L1',20),('L2',15)])]
        sub=s1(derive(study(visits,cc=112,kc=113)))
        self.assertEqual([r['classification'] for r in sub['visits']],['PR','NE','PR','PR'])
        self.assertEqual(sub['first_candidate_day'],28)
        self.assertEqual(sub['confirmed_response'],{'start_day':84,'confirmed_day':112})
        visits[1]=rec('b',56,57,1,[('L1',25),('L2',15)])
        sub=s1(derive(study(visits,cc=112,kc=113)))
        self.assertEqual(sub['visits'][1]['classification'],'SD')
        self.assertEqual(sub['confirmed_response'],{'start_day':84,'confirmed_day':112})
    def test_confirmation_survives_later_pd(self):
        visits=[rec('a',28,29,1,[('L1',20),('L2',15)]),
                rec('b',56,57,1,[('L1',18),('L2',14)]),
                rec('c',84,85,1,[('L1',40),('L2',30)])]
        sub=s1(derive(study(visits)))
        self.assertEqual(sub['confirmed_response'],{'start_day':28,'confirmed_day':56})
        self.assertEqual(sub['first_pd_day'],84)
        self.assertEqual(sub['best_confirmed'],'PR')
    def test_confirmed_cr_requires_cr_only_between(self):
        visits=[rec('a',28,29,1,[('L1',0),('L2',0)],nt='absent'),
                rec('b',56,57,1,[('L1',3),('L2',0)],nt='stable'),
                rec('c',84,85,1,[('L1',0),('L2',0)],nt='absent'),
                rec('d',112,113,1,[('L1',0),('L2',0)],nt='absent')]
        sub=s1(derive(study(visits,cc=112,kc=113)))
        self.assertEqual(sub['confirmed_response'],{'start_day':28,'confirmed_day':56})
        self.assertEqual(sub['confirmed_cr'],{'start_day':84,'confirmed_day':112})
        self.assertEqual(sub['best_confirmed'],'CR')
    def test_confirmed_cr_cutoff_boundary(self):
        visits=[rec('a',28,29,1,[('L1',20),('L2',15)]),
                rec('b',56,57,1,[('L1',0),('L2',0)],nt='absent'),
                rec('c',84,85,1,[('L1',0),('L2',0)],nt='absent')]
        sub=s1(derive(study(visits,cc=83)))
        self.assertEqual(sub['confirmed_response'],{'start_day':28,'confirmed_day':56})
        self.assertIsNone(sub['confirmed_cr'])
        self.assertEqual(sub['open_cr_segment'],{'start_day':56,'earliest_confirm_day':84})
        self.assertEqual(sub['best_confirmed'],'PR')
        sub=s1(derive(study(visits,cc=84)))
        self.assertEqual(sub['confirmed_cr'],{'start_day':56,'confirmed_day':84})
        self.assertEqual(sub['best_confirmed'],'CR')
    def test_same_version_conflict_rejected(self):
        visits=[rec('a',28,29,1,[('L1',20),('L2',15)]),rec('b',28,29,1,[('L1',21),('L2',15)])]
        with self.assertRaises(ValueError):derive(study(visits))
    def test_identical_duplicates_deduped(self):
        visits=[rec('a',28,29,1,[('L1',20),('L2',15)]),rec('b',28,29,1,[('L1',20),('L2',15)])]
        row=s1(derive(study(visits)))['visits'][0]
        self.assertEqual(row['record_ids'],['a','b'])
        self.assertEqual(row['classification'],'PR')
    def test_cutoffs_inclusive(self):
        visits=[rec('a',28,86,1,[('L1',20),('L2',15)]),rec('b',28,90,2,[('L1',26),('L2',20)])]
        row=s1(derive(study(visits,cc=28,kc=86)))['visits'][0]
        self.assertEqual((row['version'],row['sum_mm']),(1,35))
        row=s1(derive(study(visits,cc=28,kc=90)))['visits'][0]
        self.assertEqual((row['version'],row['sum_mm']),(2,46))
    def test_no_visible_version_keeps_denominator(self):
        visits=[rec('a',28,90,1,[('L1',20),('L2',15)])]
        out=derive(study(visits))
        self.assertEqual(s1(out)['visits'],[])
        self.assertTrue(any(e['version'] is None and '无可见修订' in e['reason'] for e in out['excluded']))
        self.assertEqual(out['cohort']['denominator'],2)
        self.assertEqual(out['cohort']['no_response']['count'],2)
    def test_exact_ratio_rendering(self):
        subjects=[{'subject_id':f's{i}','baseline':[{'lesion_id':'L1','mm':25}]} for i in range(3)]
        v={'subjects':subjects,'clinical_cutoff':84,'knowledge_cutoff':86,
           'visits':[rec('a',28,29,1,[('L1',10)],sid='s0'),rec('b',56,57,1,[('L1',10)],sid='s0')]}
        cohort=derive(v)['cohort']
        self.assertEqual(cohort['confirmed_response']['ratio'],'1/3')
        self.assertAlmostEqual(cohort['confirmed_response']['value'],1/3)
        self.assertEqual(cohort['no_response']['ratio'],'2/3')

if __name__=='__main__':unittest.main()
