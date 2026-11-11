"""核对靶病灶身份并清点原始随访修订。"""
from .io import unique,finite,integer

def summarize(v):
    subjects=v['subjects'];known=unique(subjects,'subject_id');lesions={}
    for s in subjects:
        ids=unique(s['baseline'],'lesion_id')
        if not ids or sum(finite(x['mm']) for x in s['baseline'])<=0 or any(x['mm']<0 for x in s['baseline']):raise ValueError('基线直径错误')
        lesions[s['subject_id']]=ids
    unique(v['visits'],'record_id');counts={s:0 for s in sorted(known)}
    for r in v['visits']:
        sid=r['subject_id']
        if sid not in known:raise ValueError('未知受试者')
        integer(r['visit_day']);integer(r['recorded_day']);integer(r['version'],1)
        ids=unique(r['targets'],'lesion_id')
        if not ids<=lesions[sid]:raise ValueError('未知靶病灶')
        for x in r['targets']:
            if x['mm'] is not None and finite(x['mm'])<0:raise ValueError('负直径')
        if r['non_target'] not in ('absent','stable','progression','unknown') or type(r['new_lesion']) is not bool:raise ValueError('访视状态错误')
        counts[sid]+=1
    return {'subjects':len(known),'raw_visit_records':len(v['visits']),'records_per_subject':counts}
