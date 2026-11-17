"""构建眼科PT/类别的全背景国家分层分母，并计算来源敏感性。"""
from pathlib import Path
import json, sys
import duckdb
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];R=HERE/'results'
OLD=ROOT/'data/raw/FAERS/tadalafil_local_results/20261006_210246'
CLASS=ROOT/'review_20261013/mural_microenvironment_extension/FAERS_class_compare'
sys.path.insert(0,str(HERE));from analyze_ocular_PT import table_stats,DRUGS,N

def norm(x):return ' '.join(str(x).strip().upper().split())

def main():
 dic=pd.read_csv(R/'ocular_PT_dictionary.csv'); dic=dic[dic.include_ocular][['pt','ocular_group']]
 ocular=set(dic.pt); group=dict(zip(dic.pt,dic.ocular_group))
 audit=json.loads((OLD/'run_audit.json').read_text(encoding='utf-8'))
 reacs=[x for x in audit['files'] if x['table']=='REAC']
 con=duckdb.connect(str(OLD/'working.duckdb'),read_only=True)
 background=[];group_union=[]
 bg_cache=R/'background_country_PT_counts.csv.gz'; gu_cache=R/'background_country_group_counts.csv'
 if bg_cache.exists() and gu_cache.exists():
  bg=pd.read_csv(bg_cache).fillna({'reporter_country':''}); gu=pd.read_csv(gu_cache).fillna({'reporter_country':''})
  print('Reusing completed ocular background cache',flush=True)
 else:
  for i,e in enumerate(reacs,1):
   cohort=con.execute('select primaryid,upper(trim(coalesce(reporter_country,\'\'))) reporter_country from cohort where source_quarter=?',[e['quarter']]).fetchdf()
   cohort.primaryid=cohort.primaryid.astype(str); ids=set(cohort.primaryid)
   parts=[]
   for ch in pd.read_csv(e['path'],sep='$',encoding='latin1',dtype=str,usecols=['primaryid','pt'],chunksize=500_000,low_memory=False,keep_default_na=False):
    ch.primaryid=ch.primaryid.str.strip();ch.pt=ch.pt.map(norm)
    x=ch[ch.primaryid.isin(ids)&ch.pt.isin(ocular)]
    if len(x):parts.append(x)
   if not parts:continue
   q=pd.concat(parts).drop_duplicates(['primaryid','pt']).merge(cohort,on='primaryid',validate='many_to_one')
   q['ocular_group']=q.pt.map(group)
   background.append(q.groupby(['reporter_country','pt']).primaryid.nunique().rename('event_reports').reset_index())
   group_union.append(q.groupby(['reporter_country','ocular_group']).primaryid.nunique().rename('event_reports').reset_index())
   anyq=q.drop_duplicates(['reporter_country','primaryid']).groupby('reporter_country').primaryid.nunique().rename('event_reports').reset_index()
   anyq['ocular_group']='ANY_OCULAR';group_union.append(anyq)
   print(f'Background ocular {i}/52',flush=True)
  bg=pd.concat(background).groupby(['reporter_country','pt'],as_index=False).event_reports.sum()
  gu=pd.concat(group_union).groupby(['reporter_country','ocular_group'],as_index=False).event_reports.sum()
  bg.to_csv(bg_cache,index=False,compression='gzip');gu.to_csv(gu_cache,index=False)
 # 国家总报告分母来自已审计主队列。
 totals=pd.read_csv(CLASS/'full_background_country_year_counts.csv').fillna({'reporter_country':''})
 totals.reporter_country=totals.reporter_country.map(norm)
 totals=totals.groupby('reporter_country',as_index=False,dropna=False).N.sum()
 assert totals.N.sum()==N
 rows=pd.read_csv(R/'three_drug_all_PT_report_rows.csv.gz',dtype={'primaryid':str},low_memory=False).merge(dic,on='pt')
 exposure=pd.read_csv(CLASS/'eligible_three_drug_flags.csv.gz',dtype={'primaryid':str},low_memory=False)
 # 聚合endpoint正式统计。
 agg=[]
 for role,flag in [('PS','primary_PS'),('ALL','primary')]:
  for drug in DRUGS:
   ex=exposure[exposure.drug.eq(drug)&exposure[flag]]; e=ex.primaryid.nunique()
   rr=rows[rows.drug.eq(drug)&rows[flag]]
   for label in ['ANY_OCULAR']+sorted(dic.ocular_group.unique()):
    a=rr.primaryid.nunique() if label=='ANY_OCULAR' else rr.loc[rr.ocular_group.eq(label),'primaryid'].nunique()
    ev=int(gu.loc[gu.ocular_group.eq(label),'event_reports'].sum())
    rec={'role':role,'drug':drug,'endpoint':label};rec.update(table_stats(a,e-a,ev-a,N-e-ev+a));agg.append(rec)
 ag=pd.DataFrame(agg);ag['BH_q_fisher']=np.nan
 for k,g in ag.groupby('role'):ag.loc[g.index,'BH_q_fisher']=multipletests(g.fisher_p,method='fdr_bh')[1]
 ag.to_csv(R/'ocular_group_disproportionality.csv',index=False)
 # 候选信号逐国剔除；同时重算四格。SA+DE固定联合敏感性也包含。
 signal=pd.read_csv(R/'ocular_PT_primary_PS_signals.csv')
 influence=[]
 for x in signal.itertuples(index=False):
  ex=exposure[exposure.drug.eq(x.drug)&exposure.primary_PS]
  er=rows[rows.drug.eq(x.drug)&rows.primary_PS&rows.pt.eq(x.pt)]
  countries=[set(['SA','DE'])]
  countries += [{c} for c,n in er.reporter_country.fillna('').map(norm).value_counts().items() if c and n>=2]
  seen=set()
  for excluded in countries:
   label='+'.join(sorted(excluded))
   if label in seen:continue
   seen.add(label)
   excluded_N=int(totals.loc[totals.reporter_country.isin(excluded),'N'].sum()); n=N-excluded_N
   ee=ex[~ex.reporter_country.fillna('').map(norm).isin(excluded)]; e=ee.primaryid.nunique()
   aa=er[~er.reporter_country.fillna('').map(norm).isin(excluded)].primaryid.nunique()
   excluded_ev=int(bg.loc[bg.pt.eq(x.pt)&bg.reporter_country.isin(excluded),'event_reports'].sum())
   ev=int(bg.loc[bg.pt.eq(x.pt),'event_reports'].sum())-excluded_ev
   rec={'drug':x.drug,'pt':x.pt,'ocular_group':x.ocular_group,'excluded_country':label,
        'original_a':int(x.a),'removed_exposed_events':int(x.a-aa)}
   rec.update(table_stats(aa,e-aa,ev-aa,n-e-ev+aa));influence.append(rec)
 pd.DataFrame(influence).to_csv(R/'signal_country_influence.csv',index=False)
 print(json.dumps({'status':'complete','background_ocular_PT_country_rows':len(bg),
  'background_group_country_rows':len(gu),'group_statistics':len(ag),'country_influence_rows':len(influence)},indent=2))

if __name__=='__main__':main()
