"""文献来源同文异写复核及报告级适应证描述；不修改任何四格表或排除规则。"""
from __future__ import annotations
import concurrent.futures as cf
import json
import re
import subprocess
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from run_class_compare import OUT, OLD, ROOT, DRUGS, ENDPOINTS, norm, save, sha

def source_key(text):
    s=norm(text)
    if not s: return '', 'missing'
    compact=re.sub(r'[^A-Z0-9]','',s)
    if '000527356' in s or ('ALSARHANI' in compact and 'CENTRALSEROUSCHORIORETINOPATHYASSOCIATEDWITHTADALAFIL' in compact):
        return 'Alsarhani_tadalafil_CSCR_10.1159_000527356','DOI_or_author_and_distinctive_title'
    if 'MOHAMMADPOUR' in compact and 'CENTRALSEROUSCHORIORETINOPATHYFOLLOWINGINGESTIONOFSILDENAFILCITRATE' in compact:
        return 'Mohammadpour_sildenafil_CSCR_2019','author_and_distinctive_title'
    if 'GARGALLO' in compact and 'BILATERAL' in compact and 'DETACHMENT' in compact and 'TADALAFIL' in compact:
        return 'Gargallo_Benedicto_tadalafil_bilateral_detachment','author_and_distinctive_title_concepts'
    if 'ETMINAN' in compact and 'RISKOF' in compact and 'OCULARADVERSEEVENTS' in compact and 'PHOSPHODIESTERASE5INHIBITORS' in compact:
        return 'Etminan_PDE5i_ocular_events_cohort_2022','author_and_distinctive_title'
    if 'TORKASHVAND' in compact and 'CENTRALRETINALVEINANDARTERYOCCLUSION' in compact and 'SILDENAFIL' in compact:
        return 'Torkashvand_sildenafil_vascular_occlusion_2023','author_and_distinctive_title'
    return 'UNMERGED_'+hashlib.sha256(s.encode()).hexdigest()[:16], 'exact_normalized_literal_only'

def scan_indi(job):
    quarter,path,keys=job
    with path.open('rb') as f: cols=f.readline().decode('latin1').rstrip('\r\n').lower().split('$')
    pat=OUT/'pattern_cache'/f'INDI_{quarter}.txt';pat.parent.mkdir(exist_ok=True)
    pat.write_text(r'^(?:'+'|'.join(sorted(keys))+r')\$'+'\n',encoding='ascii')
    r=subprocess.run(['rg','--text','--no-heading','--no-line-number','--color','never','-f',str(pat),str(path)],capture_output=True)
    if r.returncode not in (0,1):raise RuntimeError(r.stderr)
    records=[]
    for line in r.stdout.decode('latin1').splitlines():
        vals=line.split('$')
        if len(vals)!=len(cols):raise ValueError(f'Malformed INDI {path}')
        d=dict(zip(cols,vals))
        if d['primaryid'] in keys:records.append(dict(primaryid=d['primaryid'],source_quarter=quarter,indi_pt=norm(d.get('indi_pt'))))
    return records,dict(quarter=quarter,path=str(path),sha256=sha(path),matched_rows=len(records))

def main():
    ledger=pd.read_csv(OUT/'three_drug_eye_source_ledger.csv.gz',dtype={'caseid':str,'primaryid':str},keep_default_na=False)
    ledger[['literature_source_key','source_match_rule']]=ledger.lit_ref.apply(lambda z:pd.Series(source_key(z)))
    ledger.to_csv(OUT/'source_normalized_eye_ledger.csv.gz',index=False,compression='gzip')
    save('literature_mapping_dictionary.csv',ledger[['lit_ref','literature_source_key','source_match_rule']].drop_duplicates())
    rows=[]
    for role in ['PS','ALL']:
      flag='primary_PS' if role=='PS' else 'primary'
      for drug in DRUGS:
       for endpoint,col in ENDPOINTS.items():
        sub=ledger[ledger.drug.eq(drug)&ledger[flag]&ledger[col]]
        counts=sub[sub.literature_source_key.ne('')].literature_source_key.value_counts()
        for key,n in counts.items():
            rows.append(dict(role=role,drug=drug,endpoint=endpoint,literature_source_key=key,
                reports=int(n),total_endpoint_reports=len(sub),share_of_endpoint_reports=n/len(sub),
                scope='source_family_only_not_independent_patient_count'))
    save('normalized_literature_source_counts.csv',pd.DataFrame(rows))
    # 提供不会将零事件+0.5修正误显示成强阳性的展示表。
    res=pd.read_csv(OUT/'ROR_PRR_BCPNN_all_scenarios.csv')
    main=res.query("mapping=='primary' and role=='PS'").copy()
    main['display_ROR']=main.sample_ROR;main['display_OR025']=main.exact_OR025;main['display_OR975']=main.exact_OR975
    main['display_note']=np.where(main.a.eq(0),'zero_observed_events; no_positive_signal; ignore_HA_Wald_for_interpretation','sample_ROR_with_conditional_exact_OR_interval')
    main['exploratory_fisher_BH_q']=np.nan
    mask=main.scenario.eq('full')&main.endpoint.isin(['CSCR_exact','Broad_RD'])
    main.loc[mask,'exploratory_fisher_BH_q']=false_discovery_control(main.loc[mask,'fisher_p'].to_numpy())
    save('primary_PS_display_table.csv',main)
    # 报告级适应证：每报告可以包含多个适应证，不能相加得到唯一患者数。
    exp=pd.read_csv(OUT/'eligible_three_drug_flags.csv.gz',dtype={'caseid':str,'primaryid':str},keep_default_na=False)
    exp=exp[exp.primary].copy()
    keys=exp.groupby('source_quarter').primaryid.agg(set).to_dict()
    audit=json.loads((OLD/'run_audit.json').read_text('utf-8')); jobs=[]; missing=[]
    for ent in audit['files']:
        if ent['table']!='DRUG':continue
        p=Path(ent['path']); candidates=[z for z in p.parent.iterdir() if z.name.upper().startswith('INDI') and z.suffix.lower()=='.txt']
        if len(candidates)!=1:missing.append(ent['quarter']);continue
        jobs.append((ent['quarter'],candidates[0],keys.get(ent['quarter'],set())))
    rr=[]; manifests=[]
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        for i,(rows,m) in enumerate(pool.map(scan_indi,jobs),1):
            rr+=rows;manifests.append(m)
            if i%13==0:print(f'INDI source {i}/52',flush=True)
    ind=pd.DataFrame(rr).drop_duplicates()
    ind.to_csv(OUT/'selected_report_INDI_terms.csv.gz',index=False,compression='gzip')
    merged=exp.merge(ind,on=['primaryid','source_quarter'],how='left',validate='many_to_many')
    summaries=[]
    for role in ['PS','ALL']:
      flag='primary_PS' if role=='PS' else 'primary'
      for drug in DRUGS:
        d=merged[merged.drug.eq(drug)&merged[flag]]
        for term,g in d.dropna(subset=['indi_pt']).groupby('indi_pt'):
            summaries.append(dict(role=role,drug=drug,report_level_co_reported_indication=term,n_reports=g.caseid.nunique(),drug_reports=exp.loc[exp.drug.eq(drug)&exp[flag],'caseid'].nunique()))
    save('report_level_indication_distribution.csv',pd.DataFrame(summaries))
    (OUT/'source_followup_audit.json').write_text(json.dumps(dict(status='complete',INDI_missing_quarters=missing,INDI_manifest=manifests,
       source_key_scope='posthoc_bibliographic_family_descriptive_only',report_indication_scope='co_reported;not_drug_sequence_linked',
       rule_sha256=sha(OUT/'source_followup_addendum.md')),indent=2,ensure_ascii=False),encoding='utf-8')
    print('SOURCE FOLLOWUP COMPLETE',flush=True)

if __name__=='__main__':main()
