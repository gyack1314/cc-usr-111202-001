"""按截点推导受试者结局：修订选择、逐访视分类、确认缓解与固定分母比例。

研究内简化协议（不声称覆盖完整临床指南）：
- 每人固定一组基线靶病灶，基线直径总和必须为正（由domain.summarize校验）。
- 同一(受试者,访视日)可有多个修订版本：取recorded_day不晚于knowledge_cutoff的
  最高版本；同版本不同载荷视为冲突并拒绝；载荷完全相同的多条记录去重。
- 只纳入visit_day不晚于clinical_cutoff的访视；同日靶测量与非靶状态来自同一
  版本记录，不跨版本拼接。
- 完整访视：全部基线病灶均有非负测量且非靶状态已知（非unknown）。
- CR：直径总和为零且非靶全部消失；PR：非靶无进展且较基线缩小至少30%；
  PD：较截至本次（含基线与本次）的完整访视最小总和增大至少20%且至少5毫米；
  新病灶或非靶明确进展独立触发PD，即使靶测量缺失；其余完整访视为SD，
  不完整为NE。决胜顺序PD、CR、PR、SD；NE不更新nadir；首次PD后不再判定。
- 确认缓解：连续CR/PR段首日起至少28天后仍出现CR/PR，中间不能有SD/NE/PD；
  中间仍为CR/PR但间隔不足可继续等待。确认CR另要求连续CR段内两次CR间隔
  至少28天。确认访视落在截点之后则当前未确认；已完成的确认不因后续PD抹除。
- 分母固定为入组名单全体，无可用访视者保留在分母中。
- 推导是纯函数：任何修订变化都触发完整重推，不存在只改当前行的增量路径。
"""
from fractions import Fraction
from .domain import summarize
from .io import integer

RANK={'CR':0,'PR':1,'SD':2,'PD':3,'NE':4}
CONFIRM_DAYS=28

def _frac(x):
    if type(x) is int:return Fraction(x)
    if type(x) is float:return Fraction(str(x))
    raise ValueError('测量值类型错误')

def _num(f):
    return f.numerator if f.denominator==1 else float(f)

def _payload(r):
    return (r['recorded_day'],tuple(sorted((t['lesion_id'],t['mm']) for t in r['targets'])),r['non_target'],r['new_lesion'])

def _select(visits,kc):
    """按(受试者,访视日)选不晚于knowledge_cutoff的最新修订，返回(选中映射,排除清单)。"""
    groups={}
    for r in visits:groups.setdefault((r['subject_id'],r['visit_day']),{}).setdefault(r['version'],[]).append(r)
    chosen={};excluded=[]
    for sid,day in sorted(groups):
        versions={}
        for ver,recs in sorted(groups[(sid,day)].items()):
            if len({_payload(x) for x in recs})>1:raise ValueError('同一修订号内容冲突')
            versions[ver]=recs
        visible={ver:recs for ver,recs in versions.items() if recs[0]['recorded_day']<=kc}
        for ver,recs in versions.items():
            if ver not in visible:
                excluded.append({'subject_id':sid,'visit_day':day,'version':ver,'record_ids':[x['record_id'] for x in recs],
                                 'reason':f"recorded_day {recs[0]['recorded_day']} 晚于 knowledge_cutoff {kc}"})
        if not visible:
            excluded.append({'subject_id':sid,'visit_day':day,'version':None,'record_ids':[],
                             'reason':'knowledge_cutoff 前无可见修订，整访视排除'})
            continue
        best=max(visible)
        for ver,recs in visible.items():
            if ver!=best:
                excluded.append({'subject_id':sid,'visit_day':day,'version':ver,'record_ids':[x['record_id'] for x in recs],
                                 'reason':f'被同访视更晚可见修订 version {best} 取代'})
        chosen[(sid,day)]=visible[best]
    return chosen,excluded

def _classify(sid,baseline,chosen,cc,excluded):
    """按visit_day顺序逐访视分类，返回(逐访视行,基线总和,首次PD日)。"""
    base={x['lesion_id']:_frac(x['mm']) for x in baseline}
    bsum=sum(base.values(),Fraction(0))
    rows=[];nadir=bsum;first_pd=None
    for day in sorted(d for s,d in chosen if s==sid):
        recs=chosen[(sid,day)];r=recs[0]
        if day>cc:
            excluded.append({'subject_id':sid,'visit_day':day,'version':r['version'],'record_ids':[x['record_id'] for x in recs],
                             'reason':f'visit_day 晚于 clinical_cutoff {cc}'})
            continue
        measured={t['lesion_id']:t['mm'] for t in r['targets']}
        missing=[lid for lid in base if measured.get(lid) is None]
        nt=r['non_target'];complete=not missing and nt!='unknown'
        total=sum((_frac(measured[l]) for l in base),Fraction(0)) if complete else None
        notes=[]
        if missing:notes.append('缺测靶病灶:'+','.join(missing))
        if nt=='unknown':notes.append('非靶状态未知')
        row={'visit_day':day,'version':r['version'],'record_ids':[x['record_id'] for x in recs],
             'complete':complete,'sum_mm':_num(total) if total is not None else None,
             'baseline_sum_mm':_num(bsum),'non_target':nt,'new_lesion':r['new_lesion']}
        if first_pd is not None:
            row['classification']='POST_PD';row['nadir_mm']=_num(nadir)
            row['notes']=notes+[f'首次PD于第{first_pd}天后，不再判定']
            rows.append(row);continue
        ref=min(nadir,total) if complete else nadir
        if r['new_lesion'] or nt=='progression':
            cls='PD';first_pd=day
            trig=[]
            if r['new_lesion']:trig.append('新病灶')
            if nt=='progression':trig.append('非靶明确进展')
            notes.append('独立触发PD:'+'+'.join(trig))
            if complete:nadir=ref
        elif not complete:cls='NE'
        else:
            if total>=ref*Fraction(6,5) and total-ref>=5:
                cls='PD';first_pd=day;notes.append('较nadir增大至少20%且至少5毫米')
            elif total==0 and nt=='absent':cls='CR'
            elif total<=bsum*Fraction(7,10):cls='PR'
            else:cls='SD'
            nadir=ref
        row['classification']=cls;row['nadir_mm']=_num(nadir);row['notes']=notes
        if complete:
            row['pct_of_baseline']=round(float(total*100/bsum),1)
            row['increase_from_nadir_mm']=_num(total-ref)
            row['increase_from_nadir_pct']=None if ref==0 else round(float((total-ref)*100/ref),1)
        rows.append(row)
    return rows,bsum,first_pd

def _confirm(rows):
    """在连续CR/PR段上推导候选与确认日期，返回(首次候选,确认缓解,确认CR,等待中段,等待中CR段,最佳确认等级)。"""
    cls=[r['classification'] for r in rows];days=[r['visit_day'] for r in rows]
    runs=[];i=0
    while i<len(rows):
        if cls[i] in ('CR','PR'):
            j=i
            while j+1<len(rows) and cls[j+1] in ('CR','PR'):j+=1
            runs.append((i,j));i=j+1
        else:i+=1
    first_candidate=days[runs[0][0]] if runs else None
    resp=None;cr=None;last_cr_seg=None
    for a,b in runs:
        if resp is None:
            for k in range(a,b+1):
                if days[k]-days[a]>=CONFIRM_DAYS:resp={'start_day':days[a],'confirmed_day':days[k]};break
        k=a
        while k<=b:
            if cls[k]=='CR':
                m=k
                while m+1<=b and cls[m+1]=='CR':m+=1
                last_cr_seg=(k,m)
                if cr is None:
                    for t in range(k,m+1):
                        if days[t]-days[k]>=CONFIRM_DAYS:cr={'start_day':days[k],'confirmed_day':days[t]};break
                k=m+1
            else:k+=1
    open_run=None
    if resp is None and runs:
        a,_=runs[-1];open_run={'start_day':days[a],'earliest_confirm_day':days[a]+CONFIRM_DAYS}
    open_cr=None
    if cr is None and last_cr_seg is not None:
        open_cr={'start_day':days[last_cr_seg[0]],'earliest_confirm_day':days[last_cr_seg[0]]+CONFIRM_DAYS}
    best='CR' if cr else 'PR' if resp else 'NONE'
    return first_candidate,resp,cr,open_run,open_cr,best

def derive(v):
    """从原始资料完整重推截点结局，返回可JSON序列化的推导结果。"""
    summarize(v)
    cc=integer(v['clinical_cutoff']);kc=integer(v['knowledge_cutoff'])
    chosen,excluded=_select(v['visits'],kc)
    subjects=[];counts={'CR':0,'PR':0,'UNCONFIRMED':0,'NONE':0}
    for s in v['subjects']:
        sid=s['subject_id']
        rows,bsum,first_pd=_classify(sid,s['baseline'],chosen,cc,excluded)
        first_candidate,resp,cr,open_run,open_cr,best=_confirm(rows)
        seen=[r['classification'] for r in rows if r['classification']!='POST_PD']
        subjects.append({'subject_id':sid,
                         'baseline':{'lesions':[{'lesion_id':x['lesion_id'],'mm':x['mm']} for x in s['baseline']],'sum_mm':_num(bsum)},
                         'visits':rows,'first_pd_day':first_pd,
                         'best_overall_unconfirmed':min(seen,key=lambda c:RANK[c]) if seen else None,
                         'first_candidate_day':first_candidate,
                         'confirmed_response':resp,'confirmed_cr':cr,
                         'open_run':open_run,'open_cr_segment':open_cr,
                         'best_confirmed':best})
        counts['CR' if best=='CR' else 'PR' if best=='PR' else 'UNCONFIRMED' if first_candidate is not None else 'NONE']+=1
    excluded.sort(key=lambda e:(e['subject_id'],e['visit_day'],-1 if e['version'] is None else e['version']))
    n=len(subjects)
    def prop(c):
        if n==0:return {'count':c,'denominator':0,'ratio':None,'value':None}
        f=Fraction(c,n);return {'count':c,'denominator':n,'ratio':str(f),'value':float(f)}
    cohort={'denominator':n,
            'confirmed_cr':prop(counts['CR']),
            'confirmed_pr':prop(counts['PR']),
            'confirmed_response':prop(counts['CR']+counts['PR']),
            'unconfirmed_response':prop(counts['UNCONFIRMED']),
            'no_response':prop(counts['NONE']),
            'note':'确认缓解率仅含confirmed_cr与confirmed_pr；unconfirmed_response为未确认候选，不得计入确认缓解率；分母固定为入组名单全体。'}
    return {'protocol':'study-simplified-v1',
            'cutoffs':{'clinical_cutoff':cc,'knowledge_cutoff':kc},
            'subjects':subjects,'excluded':excluded,'cohort':cohort}
