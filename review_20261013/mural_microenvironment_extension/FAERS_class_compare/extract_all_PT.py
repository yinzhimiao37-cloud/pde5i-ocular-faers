"""提取三种PDE5抑制剂报告的全部PT，并计算去重主队列的全库PT分母。

运行：py review_20261014/FAERS_PDE5i_ocular_spectrum/extract_all_PT.py
输出是全部PT，不在此脚本中按结果挑眼科词。
"""
from pathlib import Path
from collections import Counter
import json
import hashlib
import re
import duckdb
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'data/raw/FAERS/tadalafil_local_results/20261006_210246'
CLASS=ROOT/'review_20261013/mural_microenvironment_extension/FAERS_class_compare'
OUT=HERE/'results'; OUT.mkdir(parents=True,exist_ok=True)
CHUNK=500_000

# 只用于生成高敏感性候选，不能冒充MedDRA官方Eye disorders SOC。
OCULAR_PATTERN=re.compile(
 r'\b(EYE|EYES|OCULAR|OPHTHALM|OPTIC|RETIN|MACUL|CHOROID|CORNE|CONJUNCTIV|EYELID|LACRIM|'
 r'LENS|IRIS|UVE|VITRE|SCLER|PUPIL|ORBIT|VISUAL|VISION|BLIND|DIPLOPIA|PHOTOPHOBIA|'
 r'PHOTOPSIA|CHROMATOPSIA|XANTHOPSIA|CYANOPSIA|DYSCHROMATOPSIA|AMAUROSIS|SCOTOMA|'
 r'MYOPIA|HYPERMETROPIA|ASTIGMATISM|PRESBYOPIA|CATARACT|GLAUCOMA|NYSTAGMUS|AMBLYOPIA|'
 r'BLEPHAR|KERAT|IRID|CYCLITIS|ENDOPHTHALMITIS|PAPILLOEDEMA|PALINOPSIA|METAMORPHOPSIA|'
 r'ENTOPSIA|ASTHENOPIA|MIOSIS|MYDRIASIS|EPIPHORA|DACRY|HORDEOLUM|CHALAZION|PINGUECULA|'
 r'PTERYGIUM|ACCOMMODATION)\w*\b',re.I)

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()

def norm_pt(x):
 # MedDRA PT按大小写不敏感处理；FAERS偶有全大写与标题式变体，不能拆成两个事件。
 return ' '.join(str(x).strip().upper().split())

def main():
 audit=json.loads((OLD/'run_audit.json').read_text(encoding='utf-8'))
 files={x['quarter']:Path(x['path']) for x in audit['files'] if x['table']=='REAC'}
 if len(files)!=52: raise ValueError(f'Expected 52 REAC files, got {len(files)}')
 exp=pd.read_csv(CLASS/'eligible_three_drug_flags.csv.gz',dtype={'primaryid':str,'caseid':str,'source_quarter':str},low_memory=False)
 exp=exp[['primaryid','caseid','source_quarter','drug','primary','primary_PS','brand_inclusive','brand_inclusive_PS',
          'reporter_country','occp_cod','year','sex','age','age_cod']]
 # 一个drug×report只保留一行；同报告多药保留多行以供逐药分析与共同报告审计。
 if exp.duplicated(['primaryid','source_quarter','drug']).any(): raise ValueError('Exposure flags are not unique')
 con=duckdb.connect(str(OLD/'working.duckdb'),read_only=True)
 cohort=con.execute('select primaryid,caseid,source_quarter from cohort').fetchdf()
 cohort.primaryid=cohort.primaryid.astype(str); cohort.source_quarter=cohort.source_quarter.astype(str)
 total_reports=len(cohort)
 bg=Counter(); exposed=[]; manifests=[]
 for i,(quarter,path) in enumerate(sorted(files.items()),1):
  selected=cohort.loc[cohort.source_quarter.eq(quarter),['primaryid','caseid']].copy()
  selected_ids=set(selected.primaryid)
  em=exp.loc[exp.source_quarter.eq(quarter)].copy()
  chunks=[]; rows_in=0; rows_selected=0
  for chunk in pd.read_csv(path,sep='$',encoding='latin1',dtype=str,usecols=['primaryid','caseid','pt'],
                           chunksize=CHUNK,low_memory=False,keep_default_na=False):
   rows_in+=len(chunk); chunk.primaryid=chunk.primaryid.str.strip(); chunk.pt=chunk.pt.map(norm_pt)
   chunk=chunk.loc[chunk.primaryid.isin(selected_ids)&chunk.pt.ne(''),['primaryid','caseid','pt']]
   rows_selected+=len(chunk); chunks.append(chunk)
  q=pd.concat(chunks,ignore_index=True).drop_duplicates(['primaryid','pt'])
  if not set(q.primaryid).issubset(selected_ids): raise AssertionError('Non-cohort reaction retained')
  bg.update(q.pt.value_counts().to_dict())
  if len(em):
   z=q.merge(em,on=['primaryid','caseid'],how='inner',validate='many_to_many')
   # q已按报告×PT唯一，em按报告×drug唯一，因此z按报告×drug×PT唯一。
   if z.duplicated(['primaryid','drug','pt']).any(): raise AssertionError('Duplicate exposed PT')
   exposed.append(z)
  manifests.append({'quarter':quarter,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path),
                    'raw_REAC_rows':rows_in,'selected_REAC_rows_before_PT_dedup':rows_selected,
                    'selected_report_PT_rows':len(q),'unique_PT':q.pt.nunique()})
  print(f'{i}/52 {quarter}: {len(q):,} selected report-PT rows',flush=True)
 all_exp=pd.concat(exposed,ignore_index=True)
 all_exp.to_csv(OUT/'three_drug_all_PT_report_rows.csv.gz',index=False,compression='gzip')
 pt=pd.DataFrame({'pt':list(bg),'background_event_reports':[bg[x] for x in bg]})
 pt['background_non_event_reports']=total_reports-pt.background_event_reports
 ecounts=all_exp.groupby(['drug','pt'],observed=True).agg(
  exposed_event_reports=('primaryid','nunique'),PS_event_reports=('primary_PS','sum'),
  allrole_event_reports=('primary','sum')).reset_index()
 # PS_event_reports/allrole_event_reports按report×drug×PT唯一，布尔求和即报告数。
 pt.merge(ecounts,on='pt',how='left').to_csv(OUT/'all_PT_counts_long.csv.gz',index=False,compression='gzip')
 observed=sorted(all_exp.pt.unique())
 candidates=pd.DataFrame({'pt':observed})
 candidates['keyword_candidate']=candidates.pt.str.upper().map(lambda x:bool(OCULAR_PATTERN.search(x)))
 candidates['include_ocular']='REVIEW_REQUIRED'
 candidates['ocular_group']='REVIEW_REQUIRED'
 candidates['review_reason']='Generated without viewing counts or disproportionality; semantic review required'
 candidates.loc[~candidates.keyword_candidate,['include_ocular','ocular_group']]='NOT_SCREENED'
 candidates.to_csv(OUT/'ocular_PT_dictionary_review_template.csv',index=False)
 pd.DataFrame(manifests).to_csv(OUT/'REAC_input_manifest.csv',index=False)
 meta={'status':'complete','window':'2013Q1-2025Q4','quarters':52,'cohort_reports':total_reports,
       'exposed_report_PT_rows':len(all_exp),'observed_PT_in_three_drug_reports':len(observed),
       'keyword_candidates':int(candidates.keyword_candidate.sum()),
       'ocular_dictionary_status':'REVIEW_REQUIRED; not official MedDRA SOC',
       'exposure_input_sha256':sha(CLASS/'eligible_three_drug_flags.csv.gz')}
 (OUT/'extraction_audit.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
 print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__':main()
