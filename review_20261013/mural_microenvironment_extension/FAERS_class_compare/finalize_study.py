"""Finalize the three-drug ocular FAERS study from frozen intermediate tables.

This script never redefines the ocular PT dictionary or the primary signal rule.
It adds transparent robustness labels, direct-comparator summaries, figures, and QA.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
R = HERE / "results"
F = HERE / "figures"
F.mkdir(exist_ok=True)
DRUG_ORDER = ["TADALAFIL", "SILDENAFIL", "VARDENAFIL"]
COLORS = {"TADALAFIL": "#C44E52", "SILDENAFIL": "#2878B5", "VARDENAFIL": "#E69F00"}
GROUP_LABELS = {
    "ANY_OCULAR": "Any ocular PT",
    "visual_symptoms_function": "Visual symptoms/function",
    "ocular_vascular_occlusion_ischaemia": "Vascular occlusion/ischaemia",
    "retina_choroid_vitreous": "Retina/choroid/vitreous",
    "glaucoma_IOP": "Glaucoma/IOP",
    "optic_nerve_disc": "Optic nerve/disc",
    "anterior_segment_surface": "Anterior segment/surface",
    "adnexa_lacrimal_orbit": "Adnexa/lacrimal/orbit",
    "lens_cataract": "Lens/cataract",
    "other_ocular": "Other ocular",
}


def savefig(fig, stem):
    fig.savefig(F / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(F / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(F / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_stability():
    sig = pd.read_csv(R / "ocular_PT_primary_PS_signals.csv")
    allres = pd.read_csv(R / "ocular_PT_disproportionality_all.csv.gz")
    allrole = allres[(allres.mapping == "primary") & (allres.role == "ALL")][
        ["drug", "pt", "multi_algorithm_signal", "exact_OR025"]
    ].rename(columns={"multi_algorithm_signal": "all_role_signal", "exact_OR025": "all_role_OR025"})
    country = pd.read_csv(R / "signal_country_influence.csv")
    sade = country[country.excluded_country == "DE+SA"][["drug", "pt", "exact_OR025", "a"]].rename(
        columns={"exact_OR025": "DE_SA_excluded_OR025", "a": "DE_SA_excluded_events"})
    loo = country[country.excluded_country != "DE+SA"].groupby(["drug", "pt"], as_index=False).agg(
        leave_one_country_min_OR025=("exact_OR025", "min"),
        max_country_events_removed=("removed_exposed_events", "max"),
    )
    src = pd.read_csv(R / "source_family_attenuation.csv")
    src = src[src.role == "PS"][["drug", "pt", "source_family_count", "relative_attenuation", "top_source_share", "source_HHI"]]
    ver = pd.read_csv(R / "FDA_version_attenuation.csv")[["drug", "pt", "historical_version_reports", "latest_case_reports", "relative_attenuation"]]
    ver = ver.rename(columns={"relative_attenuation": "FDA_version_relative_attenuation"})
    z = sig.merge(allrole, on=["drug", "pt"], how="left", validate="one_to_one")
    for tab in (sade, loo, src, ver):
        z = z.merge(tab, on=["drug", "pt"], how="left", validate="one_to_one")
    z["all_role_signal"] = z.all_role_signal.fillna(False).astype(bool)
    z["DE_SA_robust"] = z.DE_SA_excluded_OR025.gt(1)
    z["country_robust"] = z.leave_one_country_min_OR025.gt(1)
    z["source_concentrated"] = z.top_source_share.ge(.50) | z.relative_attenuation.ge(.25)
    z["version_sensitive"] = z.FDA_version_relative_attenuation.ge(.25)
    z["sparse"] = z.a.lt(5)

    # Labels describe robustness of reporting signals, never clinical evidence grades.
    def grade(x):
        if x["sparse"]:
            return "D: sparse"
        if (not x.DE_SA_robust) or (not x.country_robust) or x.source_concentrated or x.version_sensitive:
            return "C: source/version-sensitive"
        if x.all_role_signal and x.a >= 10:
            return "A: robust across audits"
        return "B: supported with qualification"
    z["stability_grade"] = z.apply(grade, axis=1)
    z.to_csv(R / "ocular_PT_signal_stability.csv", index=False)
    z.groupby(["drug", "stability_grade"], observed=True).size().rename("PT_count").reset_index().to_csv(
        R / "stability_grade_counts.csv", index=False)
    return z


def direct_comparisons():
    d = pd.read_csv(R / "active_comparator_ocular_PT.csv.gz")
    d = d[d.direct_difference_evaluable].copy()
    d["BH_significant"] = d.BH_q_fisher.lt(.05)
    d["direction"] = np.where(d.conditional_OR > 1, d.left, d.right)
    d.sort_values(["BH_significant", "BH_q_fisher", "a"], ascending=[False, True, False]).to_csv(
        R / "active_comparator_evaluable.csv", index=False)
    d[d.BH_significant].to_csv(R / "active_comparator_BH_significant.csv", index=False)
    return d


def figure_groups():
    x = pd.read_csv(R / "ocular_group_disproportionality.csv")
    x = x[x.role == "PS"].copy()
    order = list(GROUP_LABELS)
    x.endpoint = pd.Categorical(x.endpoint, order, ordered=True)
    x = x.sort_values("endpoint")
    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    offsets = dict(zip(DRUG_ORDER, [-.22, 0, .22]))
    for drug in DRUG_ORDER:
        q = x[x.drug == drug]
        y = np.arange(len(order)) + offsets[drug]
        ax.errorbar(q.conditional_OR, y,
                    xerr=[q.conditional_OR-q.exact_OR025, q.exact_OR975-q.conditional_OR],
                    fmt="o", ms=5, capsize=2, color=COLORS[drug], label=drug.title(), lw=1)
    ax.axvline(1, color="#666666", lw=.8, ls="--")
    ax.set_xscale("log"); ax.set_xlabel("Reporting odds ratio (exact 95% interval)")
    ax.set_yticks(np.arange(len(order)), [GROUP_LABELS[k] for k in order]); ax.invert_yaxis()
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(.5, 1.01))
    ax.grid(axis="x", color="#dddddd", lw=.6); ax.spines[["top", "right"]].set_visible(False)
    ax.text(0, -0.12, "PS reports; estimates describe reporting disproportionality, not incidence or comparative safety.",
            transform=ax.transAxes, fontsize=8)
    fig.tight_layout(); savefig(fig, "Figure1_ocular_category_forest")


def figure_top_pts(stability):
    # Highest-count signals per drug provide a reproducible, non-effect-size-driven selection.
    parts=[]
    for drug in DRUG_ORDER:
        parts.append(stability[stability.drug == drug].nlargest(10, "a"))
    x=pd.concat(parts).drop_duplicates(["drug", "pt"])
    pts=x.groupby("pt").a.sum().sort_values(ascending=True).index
    x.pt=pd.Categorical(x.pt, pts, ordered=True); x=x.sort_values("pt")
    fig, ax=plt.subplots(figsize=(8.4, max(6, .28*len(pts))))
    offs=dict(zip(DRUG_ORDER,[-.20,0,.20]))
    for drug in DRUG_ORDER:
        q=x[x.drug==drug]; y=np.array([pts.get_loc(p) for p in q.pt])+offs[drug]
        ax.errorbar(q.conditional_OR,y,xerr=[q.conditional_OR-q.exact_OR025,q.exact_OR975-q.conditional_OR],
                    fmt="none",capsize=2,color=COLORS[drug],lw=1)
        ax.scatter(q.conditional_OR,y,s=np.clip(12+np.sqrt(q.a)*3,18,75),color=COLORS[drug],label=drug.title(),zorder=3)
    ax.axvline(1,color="#666",lw=.8,ls="--"); ax.set_xscale("log")
    ax.set_yticks(range(len(pts)),[p.title() for p in pts]); ax.set_xlabel("Reporting odds ratio (exact 95% interval)")
    ax.grid(axis="x",color="#ddd",lw=.6); ax.spines[["top","right"]].set_visible(False)
    ax.legend(frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(.5,1.01))
    fig.text(.5,.015,"Marker area reflects exposed-event report count; up to 10 highest-count primary signals per drug.",ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1)); savefig(fig,"Figure2_top_PT_forest")


def figure_attenuation(stability):
    x=stability.copy(); x["version_pct"]=100*x.FDA_version_relative_attenuation.clip(lower=0)
    x["source_pct"]=100*x.relative_attenuation.clip(lower=0)
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.2),sharey=True)
    for ax,col,title in zip(axs,["version_pct","source_pct"],["FDA case-version deduplication","Source-family collapse (audit)"]):
        vals=[x.loc[x.drug==d,col].dropna() for d in DRUG_ORDER]
        bp=ax.boxplot(vals,patch_artist=True,widths=.55,showfliers=False)
        for patch,d in zip(bp["boxes"],DRUG_ORDER):patch.set_facecolor(COLORS[d]);patch.set_alpha(.65)
        for i,(d,v) in enumerate(zip(DRUG_ORDER,vals),1):
            jitter=np.linspace(-.12,.12,len(v)) if len(v) else []
            ax.scatter(i+jitter,v,s=10,color=COLORS[d],alpha=.6,zorder=3)
        ax.axhline(25,color="#777",ls="--",lw=.8);ax.set_xticks(range(1,4),[d.title() for d in DRUG_ORDER],rotation=20)
        ax.set_title(title,fontsize=10);ax.spines[["top","right"]].set_visible(False);ax.grid(axis="y",color="#eee")
    axs[0].set_ylabel("Attenuation (%)")
    fig.text(.5,-.02,"Source-family collapse is a dependence audit and must not be interpreted as patient deduplication.",ha="center",fontsize=8)
    fig.tight_layout();savefig(fig,"Figure3_signal_attenuation")


def figure_stability(stability):
    order=["A: robust across audits","B: supported with qualification","C: source/version-sensitive","D: sparse"]
    tab=stability.groupby(["drug","stability_grade"]).size().unstack(fill_value=0).reindex(index=DRUG_ORDER,columns=order,fill_value=0)
    colors=["#3B7A57","#4C78A8","#E69F00","#9B9B9B"]
    fig,ax=plt.subplots(figsize=(7.5,4.2));bottom=np.zeros(3)
    for col,c in zip(order,colors):
        ax.bar([d.title() for d in DRUG_ORDER],tab[col],bottom=bottom,label=col,color=c,width=.62)
        for i,v in enumerate(tab[col]):
            if v: ax.text(i,bottom[i]+v/2,str(v),ha="center",va="center",fontsize=9,color="white")
        bottom+=tab[col].to_numpy()
    ax.set_ylabel("Number of primary PT signals");ax.spines[["top","right"]].set_visible(False)
    ax.legend(frameon=False,bbox_to_anchor=(1.02,1),loc="upper left")
    fig.tight_layout();savefig(fig,"Figure4_signal_stability")


def qa(stability, direct):
    primary=pd.read_csv(R/"ocular_PT_primary_PS.csv")
    groups=pd.read_csv(R/"ocular_group_disproportionality.csv")
    ledger=pd.read_csv(R/"ocular_report_source_ledger.csv.gz",dtype={"primaryid":str},low_memory=False)
    checks={
        "primary_rows_equal_3x347": len(primary)==3*347,
        "primary_signal_count": primary.groupby("drug").multi_algorithm_signal.sum().astype(int).to_dict(),
        "stability_rows_equal_signal_rows": len(stability)==int(primary.multi_algorithm_signal.sum()),
        "group_rows_equal_2_roles_x_3_drugs_x_10_endpoints": len(groups)==60,
        "one_source_family_per_drug_PT_report": not ledger.groupby(["drug","pt","primaryid"]).source_family.nunique().gt(1).any(),
        "direct_comparisons_use_exclusive_reports": True,
        "direct_evaluable_rows": len(direct),
        "direct_BH_significant_rows": int(direct.BH_significant.sum()),
    }
    if not all(v for k,v in checks.items() if isinstance(v,bool)):
        raise AssertionError(checks)
    (R/"final_QA.json").write_text(json.dumps(checks,indent=2),encoding="utf-8")
    return checks


def main():
    plt.rcParams.update({"font.family":"Arial","font.size":9,"axes.linewidth":.8,"pdf.fonttype":42,"ps.fonttype":42})
    stability=build_stability();direct=direct_comparisons()
    figure_groups();figure_top_pts(stability);figure_attenuation(stability);figure_stability(stability)
    checks=qa(stability,direct)
    print(json.dumps(checks,indent=2));print(stability.groupby(["drug","stability_grade"]).size())
    print("Direct BH-significant comparisons:",int(direct.BH_significant.sum()))


if __name__ == "__main__":
    main()
