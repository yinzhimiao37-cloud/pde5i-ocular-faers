"""FAERS 三种 PDE5 抑制剂的报告不成比例与来源聚集比较（只读原数据库）。

运行：py review_20261013/mural_microenvironment_extension/FAERS_class_compare/run_class_compare.py
依赖：numpy pandas scipy duckdb statsmodels；系统应有 rg (ripgrep)。
中文说明：不能把报告数当独立患者数，不能把 ROR 比较当临床药物效应比较。
原始全库分母来自已审计 DuckDB；原始 DRUG 再扫描，恢复三种药物完整分母。
"""
from __future__ import annotations
import concurrent.futures as cf
import hashlib
import json
import math
import platform
import re
import subprocess
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import scipy
from scipy.special import digamma
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.contingency_tables import StratifiedTable

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
OLD = ROOT / 'data/raw/FAERS/tadalafil_local_results/20261006_210246'
sys.path.insert(0, str(ROOT / 'analysis'))
from faers_bayes import bcpnn
DRUGS = ['TADALAFIL', 'SILDENAFIL', 'VARDENAFIL']
BRANDS = {'TADALAFIL': {'CIALIS', 'ADCIRCA'}, 'SILDENAFIL': {'VIAGRA', 'REVATIO'},
          'VARDENAFIL': {'LEVITRA', 'STAXYN'}}
SALTS = {'TADALAFIL': {'TADALAFIL'}, 'SILDENAFIL': {'SILDENAFIL', 'SILDENAFIL CITRATE'},
         'VARDENAFIL': {'VARDENAFIL', 'VARDENAFIL HYDROCHLORIDE', 'VARDENAFIL HYDROCHLORIDE TRIHYDRATE'}}
ENDPOINTS = {'CSCR_exact': 'cscr_exact', 'Broad_RD': 'rd_contains', 'Serous_RD_exact': 'serous_rd_exact'}
SEED = 20261013

def norm(x):
    return ' '.join(str(x or '').upper().split())

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''): h.update(b)
    return h.hexdigest()

def scan_drug(ent):
    """rg 是完整文件的无损候选筛选；逐行按 FDA $ 解析后才判定暴露。"""
    p = Path(ent['path'])
    with p.open('rb') as f: cols = f.readline().decode('latin1').rstrip('\r\n').lower().split('$')
    pattern = '|'.join(DRUGS + sorted(set().union(*BRANDS.values())))
    r = subprocess.run(['rg', '--text', '-i', '--no-heading', '--no-line-number', '--color', 'never', pattern, str(p)], capture_output=True)
    if r.returncode not in (0, 1): raise RuntimeError(r.stderr.decode())
    rows = []
    for line in r.stdout.decode('latin1').splitlines():
        v = line.split('$')
        if len(v) != len(cols): raise ValueError(f'Malformed DRUG row {p} {v[:1]}')
        z = dict(zip(cols, v)); name, ai, role = norm(z.get('drugname')), norm(z.get('prod_ai')), norm(z.get('role_cod'))
        tokens = set(re.split(r'[;/\\,+]', ai))
        tokens = {norm(t) for t in tokens}
        for drug in DRUGS:
            # 字边界防止 CHLOROPRETADALAFIL 等类似物被当作 TADALAFIL。
            generic = bool(re.search(r'(?<![A-Z0-9])' + drug + r'(?![A-Z0-9])', name))
            ingredient = bool(tokens & SALTS[drug])
            brand = name in BRANDS[drug]
            if not (generic or ingredient or brand): continue
            rows.append(dict(primaryid=z['primaryid'], source_quarter=ent['quarter'], drug=drug,
                role=role, drugname=name, prod_ai=ai, generic=generic, ingredient=ingredient,
                brand=brand, primary=(generic or ingredient) and role != 'DN',
                brand_inclusive=(generic or ingredient or brand) and role != 'DN',
                dn=role == 'DN', nonempty_ai_without_target=bool(ai) and not ingredient))
    return rows, dict(quarter=ent['quarter'], path=str(p), bytes=p.stat().st_size,
                     sha256=sha(p), parsed_matches=len(rows), schema=cols)

def extract_demo(job):
    ent, ids = job
    p = Path(ent['path'])
    with p.open('rb') as f: cols = f.readline().decode('latin1').rstrip('\r\n').lower().split('$')
    pattern = r'^(?:' + '|'.join(sorted(ids)) + r')\$'
    r = subprocess.run(['rg', '--text', '--no-heading', '--no-line-number', '--color', 'never', pattern, str(p)], capture_output=True)
    if r.returncode not in (0, 1): raise RuntimeError(r.stderr.decode())
    rows=[]
    for line in r.stdout.decode('latin1').splitlines():
        v=line.split('$')
        if len(v)!=len(cols): raise ValueError(f'Malformed DEMO {p}')
        z=dict(zip(cols,v))
        if z['primaryid'] in ids: z['source_quarter']=ent['quarter']; rows.append(z)
    return rows

def statistics(a,b,c,d, seed=None):
    """BCPNN 后验区间与 ROR/PRR 置信区间分列；零格不藏起来。"""
    a,b,c,d=map(int,(a,b,c,d)); out=dict(a=a,b=b,c=c,d=d,N=a+b+c+d,
        exposed_reports=a+b,event_reports=a+c,zero_correction='none',status='estimable')
    if min(a+b,c+d,a+c,b+d)<=0:
        out['status']='zero_margin_not_estimable'; return out
    x=np.array([a,b,c,d],float)
    if min(x)==0: x+=.5; out['zero_correction']='Haldane_Anscombe_all_four_cells_for_Wald_only'
    A,B,C,D=x; ror=A*D/(B*C); se=np.sqrt(np.sum(1/x)); prr=(A/(A+B))/(C/(C+D)); sep=np.sqrt(1/A-1/(A+B)+1/C-1/(C+D))
    out.update(ROR_Wald=ror,ROR025_Wald=np.exp(np.log(ror)-1.95996398454*se),ROR975_Wald=np.exp(np.log(ror)+1.95996398454*se),
        PRR=prr,PRR025=np.exp(np.log(prr)-1.95996398454*sep),PRR975=np.exp(np.log(prr)+1.95996398454*sep),
        sample_ROR=a*d/(b*c) if b*c else np.inf)
    exact=odds_ratio([[a,b],[c,d]],kind='conditional'); ci=exact.confidence_interval()
    out.update(conditional_OR=exact.statistic,exact_OR025=ci.low,exact_OR975=ci.high,
               fisher_p=fisher_exact([[a,b],[c,d]]).pvalue)
    if seed is not None:
        bc=bcpnn(a,b,c,d,mc_samples=1_000_000,seed=seed)
        post=np.array(bc['dirichlet_posterior'])
        mean=(digamma(post[0])+digamma(post.sum())-digamma(post[0]+post[1])-digamma(post[0]+post[2]))/np.log(2)
        out.update(IC=mean,IC025=bc['IC025_dirichlet_MC'],IC975=bc['IC975_dirichlet_MC'],
            IC025_MCSE=bc['IC025_dirichlet_MCSE_approx'],MC_draws=1_000_000,MC_seed=seed,
            BCPNN_model='empirical-margin Dirichlet BCPNN; equal-tail 95% credible interval')
    return out

def counts(bg, exp, endpoint):
    N=int(bg.N.sum()); eventN=int(bg[endpoint].sum()); e=len(exp); a=int(exp[endpoint].sum())
    return a,e-a,eventN-a,N-e-eventN+a

def save(name, df):
    df.to_csv(OUT / name, index=False, encoding='utf-8-sig')

def main():
    audit=json.loads((OLD/'run_audit.json').read_text('utf-8')); assert audit['status']=='complete'
    entries=[e for e in audit['files'] if e['table']=='DRUG']; assert len(entries)==52
    raw_cache=OUT/'mapped_drug_rows.csv.gz'
    if raw_cache.exists():
        raw=pd.read_csv(raw_cache,dtype={'primaryid':str,'source_quarter':str},keep_default_na=False)
    else:
        rows=[]; manifests=[]
        with cf.ThreadPoolExecutor(max_workers=3) as pool:
            for i,(x,m) in enumerate(pool.map(scan_drug,entries),1):
                rows+=x; manifests.append(m); print(f'DRUG scan {i}/52 {m["quarter"]}',flush=True)
        raw=pd.DataFrame(rows); raw.to_csv(raw_cache,index=False,compression='gzip')
        (OUT/'raw_DRUG_manifest.json').write_text(json.dumps(manifests,ensure_ascii=False,indent=2),encoding='utf-8')
    # 一个版本多药物行先聚合，不按 DRUG 明细数扩大报告数。
    raw['primary_PS']=raw.primary & raw.role.eq('PS'); raw['brand_inclusive_PS']=raw.brand_inclusive & raw.role.eq('PS')
    fields=['primary','primary_PS','brand_inclusive','brand_inclusive_PS']
    flags=raw.groupby(['primaryid','source_quarter','drug'],as_index=False)[fields].max()
    con=duckdb.connect(str(OLD/'working.duckdb'),read_only=True)
    con.execute("SET memory_limit='2GB'");con.execute('SET threads=2');con.register('new_flags',flags)
    exposures=con.execute('SELECT c.*,f.drug,f.primary,f.primary_PS,f.brand_inclusive,f.brand_inclusive_PS FROM cohort c JOIN new_flags f USING(primaryid,source_quarter)').fetchdf()
    assert not exposures.duplicated(['caseid','drug']).any()
    exposures['year']=exposures.init_fda_dt.str[:4]; exposures['reporter_country']=exposures.reporter_country.fillna('').map(norm)
    exposures.to_csv(OUT/'eligible_three_drug_flags.csv.gz',index=False,compression='gzip')
    # 背景仅分组降维，不丢病例。之后各限制均从同一个全库分母重算。
    bg=con.execute('''SELECT upper(trim(reporter_country)) reporter_country,substring(init_fda_dt,1,4) AS "year",
      count(*) N,sum(cscr_exact::INT) cscr_exact,sum(rd_contains::INT) rd_contains,sum(serous_rd_exact::INT) serous_rd_exact
      FROM cohort GROUP BY 1,2''').fetchdf()
    bg.reporter_country=bg.reporter_country.fillna('').map(norm)
    assert int(bg.N.sum())==16687482
    save('full_background_country_year_counts.csv',bg)
    # 所有眼部终点恢复完整来源字段；未取非终点文献字段，故文献只做 case-only 描述。
    ev=con.execute('SELECT primaryid,caseid,source_quarter,cscr_exact,rd_contains,serous_rd_exact FROM cohort WHERE cscr_exact OR rd_contains').fetchdf()
    source_cache=OUT/'all_eye_endpoint_DEMO.csv.gz'
    if source_cache.exists(): source=pd.read_csv(source_cache,dtype=str,keep_default_na=False)
    else:
        idmap=ev.groupby('source_quarter').primaryid.agg(set).to_dict()
        jobs=[(e,idmap[e['quarter']]) for e in audit['files'] if e['table']=='DEMO' and e['quarter'] in idmap]
        source_rows=[]
        with cf.ThreadPoolExecutor(max_workers=3) as pool:
            for i,rows in enumerate(pool.map(extract_demo,jobs),1):
                source_rows+=rows
                if i%13==0: print(f'DEMO recovery {i}/52',flush=True)
        source=pd.DataFrame(source_rows).fillna('');source.to_csv(source_cache,index=False,compression='gzip')
    assert len(source)==len(ev) and not source.duplicated(['primaryid','source_quarter']).any()
    srcfields=[k for k in ['primaryid','source_quarter','lit_ref','auth_num','mfr_num','mfr_sndr','occr_country','rept_cod'] if k in source]
    er=exposures[exposures.cscr_exact|exposures.rd_contains].merge(source[srcfields],on=['primaryid','source_quarter'],validate='many_to_one')
    er.to_csv(OUT/'three_drug_eye_source_ledger.csv.gz',index=False,compression='gzip')
    # 暴露与映射审计，包括品牌敏感性新增数、DN行、同时暴露。
    save('mapping_summary.csv', exposures.groupby('drug')[fields].sum().reset_index())
    dictionary=raw.groupby(['drug','drugname','prod_ai','generic','ingredient','brand','dn','nonempty_ai_without_target'],dropna=False).size().reset_index(name='raw_rows')
    save('drug_dictionary_audit.csv',dictionary)
    results=[]; influence=[]; concentrations=[]; sourcedist=[]
    for definition in ['primary','brand_inclusive']:
      for role in ['PS','ALL']:
        flag=definition+('_PS' if role=='PS' else '')
        for drug in DRUGS:
          exp=exposures[exposures.drug.eq(drug)&exposures[flag]].copy()
          for endpoint,eventcol in ENDPOINTS.items():
            for scenario,excluded in [('full',set()),('exclude_SA_DE',{'SA','DE'})]:
              bsub=bg[~bg.reporter_country.isin(excluded)]; esub=exp[~exp.reporter_country.isin(excluded)]
              rec=dict(mapping=definition,role=role,drug=drug,endpoint=endpoint,scenario=scenario)
              rec.update(statistics(*counts(bsub,esub,eventcol),seed=SEED+len(results)));results.append(rec)
            # 同一固定主映射下，所有贡献 >=2 事件报告的国家都逐个排除，不能挑结果。
            if definition=='primary':
              for country,n in exp.loc[exp[eventcol],'reporter_country'].value_counts().items():
                if n<2 or not country: continue
                rec=dict(mapping=definition,role=role,drug=drug,endpoint=endpoint,excluded_reporter_country=country,excluded_exposed_events=int(n))
                rec.update(statistics(*counts(bg[bg.reporter_country.ne(country)],exp[exp.reporter_country.ne(country)],eventcol)))
                influence.append(rec)
            event=er[er.drug.eq(drug)&er[flag]&er[eventcol]]
            for field in ['reporter_country','year','source_quarter','occp_cod','lit_ref','auth_num']:
              if field not in event: continue
              vals=event[field].fillna('').map(norm); known=vals[vals.ne('')]; ct=known.value_counts(); n=len(event); nk=len(known)
              concentrations.append(dict(mapping=definition,role=role,drug=drug,endpoint=endpoint,field=field,
                  event_reports=n,nonmissing_reports=nk,missing_reports=n-nk,distinct_nonmissing=len(ct),
                  top_nonmissing_count=int(ct.iloc[0]) if nk else 0,top_share_all_reports=float(ct.iloc[0]/n) if nk else np.nan,
                  HHI_nonmissing=float(((ct/nk)**2).sum()) if nk else np.nan,
                  top_value=str(ct.index[0]) if nk else '',scope='endpoint_only_source_description'))
              for value,num in ct.items():
                sourcedist.append(dict(mapping=definition,role=role,drug=drug,endpoint=endpoint,field=field,value=value,n=int(num)))
    result=pd.DataFrame(results);save('ROR_PRR_BCPNN_all_scenarios.csv',result)
    save('country_leave_one_out.csv',pd.DataFrame(influence));save('source_concentration_summary.csv',pd.DataFrame(concentrations));save('source_distribution.csv',pd.DataFrame(sourcedist))
    print(result.query("mapping=='primary' and role=='PS'")[['drug','endpoint','scenario','a','exposed_reports','sample_ROR','exact_OR025','exact_OR975']].to_string(index=False),flush=True)
    # active comparator 的所有 PS/all-role 同时暴露都按主成分 mapping 排除。
    comps=[]; multi=[]
    for role in ['PS','ALL']:
      flag='primary_PS' if role=='PS' else 'primary'
      anymap=exposures[exposures.primary].groupby('caseid').drug.nunique()
      overlapping=set(anymap[anymap>1].index)
      exp=exposures[exposures[flag]&~exposures.caseid.isin(overlapping)].copy()
      multi.append(dict(role=role,any_role_multiple_study_drug_reports=len(overlapping),eligible_single_study_drug_reports=len(exp)))
      for left,right in [('TADALAFIL',['SILDENAFIL']),('TADALAFIL',['VARDENAFIL']),('SILDENAFIL',['VARDENAFIL']),('TADALAFIL',['SILDENAFIL','VARDENAFIL'])]:
        for endpoint,eventcol in ENDPOINTS.items():
          for scenario,excluded in [('full',set()),('exclude_SA_DE',{'SA','DE'})]:
            x=exp[exp.drug.isin([left]+right)&~exp.reporter_country.isin(excluded)].copy(); x['left']=x.drug.eq(left)
            a=int((x.left&x[eventcol]).sum());b=int((x.left&~x[eventcol]).sum());c=int((~x.left&x[eventcol]).sum());d=int((~x.left&~x[eventcol]).sum())
            rec=dict(role=role,drug=left,comparator='+'.join(right),endpoint=endpoint,scenario=scenario)
            rec.update(statistics(a,b,c,d))
            tables=[]
            for key,g in x.groupby(['reporter_country','year'],dropna=False):
              t=np.array([[int((g.left&g[eventcol]).sum()),int((g.left&~g[eventcol]).sum())],
                          [int((~g.left&g[eventcol]).sum()),int((~g.left&~g[eventcol]).sum())]])
              if min(t.sum(axis=0))>0 and min(t.sum(axis=1))>0: tables.append(t)
            rec['informative_country_year_strata']=len(tables)
            if tables:
              mh=StratifiedTable(np.stack(tables,axis=2),shift_zeros=False)
              with np.errstate(divide='ignore',invalid='ignore'):
                lo,hi=mh.oddsratio_pooled_confint(); rec.update(MH_OR=mh.oddsratio_pooled,MH_OR025=lo,MH_OR975=hi,MH_p=mh.test_null_odds().pvalue)
            comps.append(rec)
    save('active_comparator_results.csv',pd.DataFrame(comps));save('multi_drug_co_report_audit.csv',pd.DataFrame(multi))
    meta=dict(status='complete',quarters=52,N=int(bg.N.sum()),endpoint_counts={k:int(bg[v].sum()) for k,v in ENDPOINTS.items()},
      primary_role='PS',mapped_exposure='generic word or exact declared ingredient token; exclude DN',
      source_concentration_scope='All selected CSCR/RD reports, not full-cohort publication linkage',
      full_background_exclusions='Reporter country only; every four-cell count recalculated',
      patient_independence='unresolved; no new cross-CASEID duplicate collapse',
      BCPNN_CrI='1,000,000 draws, empirical-margin Dirichlet prior, conditional on observed reports',
      software={'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'duckdb':duckdb.__version__},
      protocol_sha256=sha(OUT/'protocol_locked.md'),script_sha256=sha(Path(__file__)),input_audit_sha256=sha(OLD/'run_audit.json'))
    (OUT/'analysis_audit.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
    print('COMPLETE',flush=True)

if __name__=='__main__': main()
