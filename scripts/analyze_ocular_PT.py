"""三种PDE5抑制剂全部眼科PT的报告不成比例及直接药间比较。

输入词典在结果生成前冻结。主表以PS暴露为主，all-role为敏感性。
"""
from pathlib import Path
import hashlib, json, math, sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.multitest import multipletests

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]; R=HERE/'results'
CLASS=ROOT/'review_20261013/mural_microenvironment_extension/FAERS_class_compare'
OLD=ROOT/'data/raw/FAERS/tadalafil_local_results/20261006_210246'
sys.path.insert(0,str(ROOT/'analysis'))
from faers_bayes import bcpnn
N=16_687_482; DRUGS=['TADALAFIL','SILDENAFIL','VARDENAFIL']
MC=200_000; SEED=20261014

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def table_stats(a,b,c,d,mc=False,seed=SEED):
 a,b,c,d=map(int,(a,b,c,d)); rec={'a':a,'b':b,'c':c,'d':d,'N':a+b+c+d,
  'exposed_reports':a+b,'event_reports':a+c,'status':'estimable','zero_correction':'none'}
 if min(a+b,c+d,a+c,b+d)<=0:
  rec['status']='zero_margin_not_estimable'; return rec
 exact=odds_ratio([[a,b],[c,d]],kind='conditional'); ci=exact.confidence_interval(.95)
 rec.update(sample_ROR=(a*d/(b*c) if b*c else (math.inf if a*d else 0.0)),
            conditional_OR=float(exact.statistic),exact_OR025=float(ci.low),exact_OR975=float(ci.high),
            fisher_p=float(fisher_exact([[a,b],[c,d]]).pvalue))
 if a>0 and c>0:
  prr=(a/(a+b))/(c/(c+d)); se=np.sqrt(1/a-1/(a+b)+1/c-1/(c+d))
  rec.update(PRR=float(prr),PRR025=float(np.exp(np.log(prr)-1.95996398454*se)),
             PRR975=float(np.exp(np.log(prr)+1.95996398454*se)))
 else: rec.update(PRR=np.nan,PRR025=np.nan,PRR975=np.nan)
 # 全表用可复算解析BCPNN；满足其他条件者追加Dirichlet MC区间。
 bc=bcpnn(a,b,c,d,mc_samples=1000,seed=seed,batch_size=1000)
 rec.update(IC_legacy=bc['IC_mean_legacy_normal'],IC025_legacy=bc['IC025_legacy_normal'],
            IC975_legacy=bc['IC975_legacy_normal'])
 if mc:
  bc=bcpnn(a,b,c,d,mc_samples=MC,seed=seed,batch_size=20_000)
  rec.update(IC=bc['IC_mean_dirichlet_MC'],IC025=bc['IC025_dirichlet_MC'],
             IC975=bc['IC975_dirichlet_MC'],IC_MCSE025=bc['IC025_dirichlet_MCSE_approx'],
             IC_MC_draws=MC,IC_MC_seed=seed)
 return rec

def main():
 dictionary=pd.read_csv(R/'ocular_PT_dictionary.csv')
 ocular=dictionary.loc[dictionary.include_ocular.astype(str).str.lower().eq('true'),['pt','ocular_group']]
 rows=pd.read_csv(R/'three_drug_all_PT_report_rows.csv.gz',dtype={'primaryid':str,'caseid':str},low_memory=False)
 rows=rows.merge(ocular,on='pt',how='inner',validate='many_to_one')
 exposure=pd.read_csv(CLASS/'eligible_three_drug_flags.csv.gz',dtype={'primaryid':str,'caseid':str},low_memory=False)
 expden={
  ('primary','PS'):exposure[exposure.primary_PS].groupby('drug').primaryid.nunique().to_dict(),
  ('primary','ALL'):exposure[exposure.primary].groupby('drug').primaryid.nunique().to_dict(),
  ('brand_inclusive','PS'):exposure[exposure.brand_inclusive_PS].groupby('drug').primaryid.nunique().to_dict(),
  ('brand_inclusive','ALL'):exposure[exposure.brand_inclusive].groupby('drug').primaryid.nunique().to_dict()}
 # 全库PT分母来自已按报告×PT去重的all_PT_counts；同一pt的值对三药相同。
 counts=pd.read_csv(R/'all_PT_counts_long.csv.gz')
 background=counts.drop_duplicates('pt').set_index('pt').background_event_reports.to_dict()
 results=[]
 for mapping in ['primary','brand_inclusive']:
  for role in ['PS','ALL']:
   flag=mapping+('_PS' if role=='PS' else '')
   for drug in DRUGS:
    e=int(expden[(mapping,role)].get(drug,0)); sub=rows[rows.drug.eq(drug)&rows[flag]]
    ec=sub.groupby('pt').primaryid.nunique()
    for ix,x in ocular.iterrows():
     pt=x.pt; a=int(ec.get(pt,0)); ev=int(background[pt])
     rec={'mapping':mapping,'role':role,'drug':drug,'pt':pt,'ocular_group':x.ocular_group}
     # 先解析计算；MC仅在其他三指标均有正向下限且a>=3时运行。
     s=table_stats(a,e-a,ev-a,N-e-ev+a,mc=False,seed=SEED+len(results))
     pre=(a>=3 and s.get('exact_OR025',0)>1 and s.get('PRR025',0)>1 and s.get('IC025_legacy',-np.inf)>0)
     if pre: s=table_stats(a,e-a,ev-a,N-e-ev+a,mc=True,seed=SEED+len(results))
     rec.update(s); results.append(rec)
 out=pd.DataFrame(results)
 out['BH_q_fisher']=np.nan
 for key,g in out.groupby(['mapping','role','drug'],sort=False):
  out.loc[g.index,'BH_q_fisher']=multipletests(g.fisher_p.fillna(1),method='fdr_bh')[1]
 out['multi_algorithm_signal']=(out.a>=3)&(out.exact_OR025>1)&(out.PRR025>1)&(out.IC025.fillna(-np.inf)>0)&(out.BH_q_fisher<.05)
 out.to_csv(R/'ocular_PT_disproportionality_all.csv.gz',index=False,compression='gzip')
 primary=out[(out.mapping=='primary')&(out.role=='PS')].copy()
 primary.to_csv(R/'ocular_PT_primary_PS.csv',index=False)
 primary[primary.multi_algorithm_signal].sort_values(['drug','BH_q_fisher','a'],ascending=[True,True,False]).to_csv(R/'ocular_PT_primary_PS_signals.csv',index=False)
 # 任一眼科事件和九组事件均按report而不是PT行计数。
 aggregate=[]
 cats=[('ANY_OCULAR',None)]+[(g,g) for g in sorted(ocular.ocular_group.unique())]
 for role,flag in [('PS','primary_PS'),('ALL','primary')]:
  for drug in DRUGS:
   e=int(expden[('primary',role)][drug]); drows=rows[rows.drug.eq(drug)&rows[flag]]
   for label,cat in cats:
    ids=set(drows.primaryid if cat is None else drows.loc[drows.ocular_group.eq(cat),'primaryid'])
    # 全背景聚合终点从全部暴露药报告不足以得到；用全报告PT行重建另脚本前，标记待计算。
    aggregate.append({'role':role,'drug':drug,'endpoint':label,'exposed_reports':e,
                      'exposed_event_reports':len(ids),'background_event_reports':'REQUIRES_FULL_REPORT_UNION'})
 pd.DataFrame(aggregate).to_csv(R/'ocular_group_exposed_counts.csv',index=False)
 # 三药排他报告的直接PT比较；同一报告含多药全部排除。
 primary_exp=exposure[exposure.primary_PS].copy()
 ndrug=primary_exp.groupby('caseid').drug.nunique(); overlap=set(ndrug[ndrug>1].index)
 exclusive=primary_exp[~primary_exp.caseid.isin(overlap)]
 erows=rows[rows.primary_PS&~rows.caseid.isin(overlap)]
 comps=[]
 pairs=[('TADALAFIL','SILDENAFIL'),('TADALAFIL','VARDENAFIL'),('SILDENAFIL','VARDENAFIL')]
 for left,right in pairs:
  el=set(exclusive.loc[exclusive.drug.eq(left),'primaryid']); er=set(exclusive.loc[exclusive.drug.eq(right),'primaryid'])
  for x in ocular.itertuples(index=False):
   z=erows[erows.pt.eq(x.pt)]; a=z[z.drug.eq(left)].primaryid.nunique(); c=z[z.drug.eq(right)].primaryid.nunique()
   rec={'left':left,'right':right,'pt':x.pt,'ocular_group':x.ocular_group}
   rec.update(table_stats(a,len(el)-a,c,len(er)-c,mc=False)); comps.append(rec)
 comp=pd.DataFrame(comps); comp['BH_q_fisher']=np.nan
 for key,g in comp.groupby(['left','right']): comp.loc[g.index,'BH_q_fisher']=multipletests(g.fisher_p.fillna(1),method='fdr_bh')[1]
 comp['direct_difference_evaluable']=(comp.a>=3)&(comp.c>=3)
 comp.to_csv(R/'active_comparator_ocular_PT.csv.gz',index=False,compression='gzip')
 meta={'status':'complete','N':N,'ocular_PT':len(ocular),'primary_PS_signal_counts':
       primary.groupby('drug').multi_algorithm_signal.sum().astype(int).to_dict(),
       'BCPNN':'legacy normal interval for all rows; 200,000-draw Dirichlet interval for prefiltered rows',
       'MC_prefilter_is_not_final_signal':True,'dictionary_sha256':sha(R/'ocular_PT_dictionary.csv'),
       'input_report_PT_sha256':sha(R/'three_drug_all_PT_report_rows.csv.gz')}
 (R/'analysis_audit.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
 print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
