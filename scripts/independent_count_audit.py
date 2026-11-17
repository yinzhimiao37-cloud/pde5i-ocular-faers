"""Independent count-level reconstruction for the final ocular FAERS tables."""
from pathlib import Path
import json
import pandas as pd

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]; R=HERE/"results"
CLASS=ROOT/"review_20261013/mural_microenvironment_extension/FAERS_class_compare"

def main():
    dictionary=pd.read_csv(R/"ocular_PT_dictionary.csv")
    ocular=dictionary[dictionary.include_ocular.astype(str).str.lower().eq("true")][["pt","ocular_group"]]
    rows=pd.read_csv(R/"three_drug_all_PT_report_rows.csv.gz",dtype={"primaryid":str,"caseid":str},low_memory=False)
    rows=rows.merge(ocular,on="pt",how="inner",validate="many_to_one")
    rows=rows[rows.primary_PS].drop_duplicates(["drug","primaryid","pt"])
    tab=pd.read_csv(R/"ocular_PT_primary_PS.csv")
    observed=rows.groupby(["drug","pt"]).primaryid.nunique().rename("a_recomputed").reset_index()
    chk=tab.merge(observed,on=["drug","pt"],how="left",validate="one_to_one").fillna({"a_recomputed":0})
    pt_counts_match=chk.a.eq(chk.a_recomputed)

    group=pd.read_csv(R/"ocular_group_disproportionality.csv")
    group=group[group.role.eq("PS")][["drug","endpoint","a"]]
    parts=[]
    for drug,g in rows.groupby("drug"):
        parts.append({"drug":drug,"endpoint":"ANY_OCULAR","a_recomputed":g.primaryid.nunique()})
        for name,h in g.groupby("ocular_group"):
            parts.append({"drug":drug,"endpoint":name,"a_recomputed":h.primaryid.nunique()})
    gchk=group.merge(pd.DataFrame(parts),on=["drug","endpoint"],how="left",validate="one_to_one").fillna({"a_recomputed":0})
    group_counts_match=gchk.a.eq(gchk.a_recomputed)

    exposure=pd.read_csv(CLASS/"eligible_three_drug_flags.csv.gz",dtype={"primaryid":str,"caseid":str},low_memory=False)
    pe=exposure[exposure.primary_PS]
    overlap=set(pe.groupby("caseid").drug.nunique().loc[lambda x:x>1].index)
    exclusive=pe[~pe.caseid.isin(overlap)]
    erows=rows[~rows.caseid.isin(overlap)]
    direct=pd.read_csv(R/"active_comparator_ocular_PT.csv.gz")
    mismatches=[]
    for z in direct.itertuples(index=False):
        a=erows[(erows.drug==z.left)&(erows.pt==z.pt)].primaryid.nunique()
        c=erows[(erows.drug==z.right)&(erows.pt==z.pt)].primaryid.nunique()
        nl=exclusive[exclusive.drug==z.left].primaryid.nunique()
        nr=exclusive[exclusive.drug==z.right].primaryid.nunique()
        if (a,nl-a,c,nr-c)!=(z.a,z.b,z.c,z.d): mismatches.append([z.left,z.right,z.pt])

    out={
        "status":"pass" if pt_counts_match.all() and group_counts_match.all() and not mismatches else "fail",
        "PT_cells_checked":int(len(chk)),"PT_count_mismatches":int((~pt_counts_match).sum()),
        "group_cells_checked":int(len(gchk)),"group_count_mismatches":int((~group_counts_match).sum()),
        "direct_cells_checked":int(len(direct)),"direct_table_mismatches":len(mismatches),
        "multi_drug_CASEIDs_excluded_from_direct_comparisons":len(overlap),
    }
    (R/"independent_count_audit.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
    if out["status"]!="pass": raise AssertionError(out)

if __name__=="__main__": main()
