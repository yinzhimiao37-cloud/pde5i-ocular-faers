"""眼科PT的FDA版本衰减、国家影响及文献/authority来源聚集审计。"""
from pathlib import Path
from collections import defaultdict
import concurrent.futures as cf
import hashlib, json, re, subprocess
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]; R=HERE/'results'
OLD=ROOT/'data/raw/FAERS/tadalafil_local_results/20261006_210246'
CLASS=ROOT/'review_20261013/mural_microenvironment_extension/FAERS_class_compare'

def norm(x): return ' '.join(str(x or '').upper().split())
def lit_family(x):
 u=norm(x)
 if not u:return ''
 rules=[
  ('ALSARHANI_CSCR_TADALAFIL',['ALSARHANI','CENTRAL SEROUS']),
  ('GARGALLO_TADALAFIL_RD',['GARGALLO','RETINAL','DETACHMENT']),
  ('ETMINAN_PDE5I_OCULAR',['ETMINAN','RISK OF OCULAR ADVERSE']),
  ('MOHAMMADPOUR_SILDENAFIL_CSCR',['MOHAMMADPOUR','CENTRAL SEROUS']),
 ]
 for key,tokens in rules:
  if all(t in u for t in tokens): return key
 # 未经人工确认的文献不做模糊合并；保留规范化字面值。
 return 'EXACT_'+hashlib.sha256(u.encode()).hexdigest()[:20]

def extract_demo(job):
 ent,ids=job; path=Path(ent['path'])
 with path.open('rb') as f: cols=f.readline().decode('latin1').rstrip('\r\n').lower().split('$')
 pattern=r'^(?:'+'|'.join(sorted(map(re.escape,ids)))+r')\$'
 run=subprocess.run(['rg','--text','--no-heading','--no-line-number','--color','never',pattern,str(path)],capture_output=True)
 if run.returncode not in (0,1):raise RuntimeError(run.stderr.decode(errors='replace'))
 rows=[]
 for line in run.stdout.decode('latin1').splitlines():
  v=line.split('$')
  if len(v)!=len(cols):raise ValueError(f'Malformed DEMO row {path}')
  z=dict(zip(cols,v))
  if z.get('primaryid') in ids:
   z['source_quarter']=ent['quarter'];rows.append(z)
 return rows

def historical_reac(job):
 quarter,path,ids,ocular=job
 chunks=[]
 for ch in pd.read_csv(path,sep='$',encoding='latin1',dtype=str,usecols=['primaryid','pt'],
                       chunksize=500_000,low_memory=False,keep_default_na=False):
  ch.primaryid=ch.primaryid.str.strip(); ch.pt=ch.pt.map(norm)
  x=ch[ch.primaryid.isin(ids)&ch.pt.isin(ocular)]
  if len(x):chunks.append(x)
 if not chunks:return pd.DataFrame(columns=['primaryid','pt','source_quarter'])
 x=pd.concat(chunks).drop_duplicates(['primaryid','pt']);x['source_quarter']=quarter;return x

def main():
 dic=pd.read_csv(R/'ocular_PT_dictionary.csv'); ocular=set(dic.loc[dic.include_ocular,'pt'])
 final=pd.read_csv(R/'three_drug_all_PT_report_rows.csv.gz',dtype={'primaryid':str,'caseid':str},low_memory=False)
 final=final[final.pt.isin(ocular)]
 audit=json.loads((OLD/'run_audit.json').read_text(encoding='utf-8'))
 demos=[x for x in audit['files'] if x['table']=='DEMO']; reacs={x['quarter']:x['path'] for x in audit['files'] if x['table']=='REAC'}
 idsmap=final.groupby('source_quarter').primaryid.agg(lambda x:set(x)).to_dict()
 jobs=[(e,idsmap[e['quarter']]) for e in demos if e['quarter'] in idsmap]
 source=[]
 with cf.ThreadPoolExecutor(max_workers=3) as pool:
  for i,x in enumerate(pool.map(extract_demo,jobs),1): source+=x; print(f'DEMO {i}/{len(jobs)}',flush=True)
 demo=pd.DataFrame(source).fillna('')
 if demo.duplicated(['primaryid','source_quarter']).any():raise ValueError('Duplicate selected DEMO')
 keep=['primaryid','source_quarter']+[c for c in ['lit_ref','auth_num','mfr_num','mfr_sndr','occr_country','rept_cod'] if c in demo]
 ledger=final.merge(demo[keep],on=['primaryid','source_quarter'],how='left',validate='many_to_one')
 for c in keep[2:]:ledger[c]=ledger[c].map(norm)
 ledger['literature_family']=ledger.get('lit_ref','').map(lit_family)
 # 优先已核查/精确文献，再authority；两者均无则保持报告级，不猜患者重复。
 ledger['source_family']=np.where(ledger.literature_family.ne(''),'LIT:'+ledger.literature_family,
                           np.where(ledger.get('auth_num','').ne(''),'AUTH:'+ledger.auth_num,'REPORT:'+ledger.caseid))
 ledger.to_csv(R/'ocular_report_source_ledger.csv.gz',index=False,compression='gzip')
 conc=[]; attenuation=[]
 for role,flag in [('PS','primary_PS'),('ALL','primary')]:
  for (drug,pt),g in ledger[ledger[flag]].groupby(['drug','pt'],observed=True):
   n=g.primaryid.nunique(); fam=g.drop_duplicates(['primaryid','source_family']).source_family.value_counts()
   # 每报告只对应一个family；保险断言防止多PT行以外的重复。
   one=g.drop_duplicates('primaryid'); fam=one.source_family.value_counts(); collapsed=len(fam)
   attenuation.append({'role':role,'drug':drug,'pt':pt,'report_count':n,'source_family_count':collapsed,
      'absolute_attenuation':n-collapsed,'relative_attenuation':(n-collapsed)/n if n else np.nan,
      'top_source_count':int(fam.iloc[0]),'top_source_share':float(fam.iloc[0]/n),
      'source_HHI':float(((fam/n)**2).sum())})
   for field in ['reporter_country','year','source_quarter','occp_cod','literature_family','auth_num','mfr_num']:
    vals=one[field].fillna('').map(norm); known=vals[vals.ne('')]; ct=known.value_counts()
    conc.append({'role':role,'drug':drug,'pt':pt,'field':field,'reports':n,'nonmissing':len(known),
      'distinct_nonmissing':len(ct),'top_value':ct.index[0] if len(ct) else '',
      'top_count':int(ct.iloc[0]) if len(ct) else 0,'top_share_all':float(ct.iloc[0]/n) if len(ct) else np.nan,
      'HHI_nonmissing':float(((ct/len(known))**2).sum()) if len(known) else np.nan})
 pd.DataFrame(attenuation).to_csv(R/'source_family_attenuation.csv',index=False)
 pd.DataFrame(conc).to_csv(R/'source_concentration.csv.gz',index=False,compression='gzip')

 # 历史版本：严格使用相同primary成分映射和PS角色；统计report version×PT，不是病例/患者。
 raw=pd.read_csv(CLASS/'mapped_drug_rows.csv.gz',dtype=str,keep_default_na=False,low_memory=False)
 hist=raw[raw.primary.astype(str).str.lower().eq('true')&raw.role.eq('PS')][['primaryid','source_quarter','drug']].drop_duplicates()
 hids=hist.groupby('source_quarter').primaryid.agg(set).to_dict()
 hjobs=[(q,reacs[q],ids,ocular) for q,ids in hids.items()]
 hparts=[]
 with cf.ThreadPoolExecutor(max_workers=3) as pool:
  for i,x in enumerate(pool.map(historical_reac,hjobs),1):hparts.append(x);print(f'Historical REAC {i}/{len(hjobs)}',flush=True)
 hr=pd.concat(hparts,ignore_index=True).merge(hist,on=['primaryid','source_quarter'],validate='many_to_many')
 hr=hr.drop_duplicates(['primaryid','drug','pt'])
 hc=hr.groupby(['drug','pt']).primaryid.nunique().rename('historical_version_reports')
 fc=final[final.primary_PS].groupby(['drug','pt']).primaryid.nunique().rename('latest_case_reports')
 va=pd.concat([hc,fc],axis=1).fillna(0).astype(int).reset_index()
 va['absolute_attenuation']=va.historical_version_reports-va.latest_case_reports
 va['relative_attenuation']=np.where(va.historical_version_reports>0,va.absolute_attenuation/va.historical_version_reports,np.nan)
 va.to_csv(R/'FDA_version_attenuation.csv',index=False)
 meta={'status':'complete','ocular_source_reports':int(final.primaryid.nunique()),
       'source_family_rule':'literature family if nonempty; otherwise authority; otherwise report',
       'source_family_is_not_patient_deduplication':True,'historical_version_rows':len(hr),
       'FDA_version_attenuation_is_not_patient_deduplication':True}
 (R/'source_version_audit.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
 print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
