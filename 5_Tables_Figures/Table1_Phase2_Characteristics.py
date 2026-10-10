"""Table 1: Phase 2 participant characteristics by multidimensional schema group (quantitative CT variant),
in RNA-seq participants with complete Phase 2 classification and model covariates (n = 3,751)."""
import pandas as pd, numpy as np, csv
from scipy import stats
from copdgene_phenotypes import load
from multidimensional_schema import schema, G

#### Edit paths to match your file system ####
RNA_MASTER = 'path/to/master.file.freeze5.phase2.csv'
OUT_FILE = 'path/to/output/Table1_PHASE2classified_quantCT.csv'

df=load()
master=pd.read_csv(RNA_MASTER)
rna=set(master.loc[master['final.analyzed']==1,'aliasID'].astype(str).str.strip())
p2=df[df.Phase_study==2].drop_duplicates('sid').set_index('sid')
diffs=['neutrophl_pct','lymphcyt_pct','monocyt_pct','eosinphl_pct']

d=p2.loc[[s for s in p2.index if s in rna]].copy()
d['group']=schema(d,'quant')
keep=d['group'].notna() & d['smoking_status'].notna() & d[diffs].notna().all(axis=1) & d['age_visit'].notna()
d=d[keep].copy()
d['group']=pd.Categorical(d['group'],categories=G,ordered=True)

d['female']=(d['gender'].astype(str).str.strip().isin(['2','2.0'])).astype(float)
d['aa']=(d['race'].astype(str).str.strip().isin(['2','2.0'])).astype(float)
d['current_smoker']=(d['smoking_status']==2).astype(float)
d['mmrc2']=(d['MMRCDyspneaScor']>=2).astype(float)
d['cb']=(d['Chronic_Bronchitis']==1).astype(float)
d['distwalked_m']=pd.to_numeric(d['distwalked'],errors='coerce')*0.3048
gold=pd.to_numeric(d['finalgold_visit'],errors='coerce')
d['gold_cat']=np.select([gold<=0,gold==1,gold==2,gold>=3],['GOLD 0','GOLD 1','GOLD 2','GOLD 3-4'],default=None)

SPEC=[('SECTION','CLINICAL'),('CONT','age_visit','Age, years',1),('CAT','female','Female, n (%)'),
 ('CAT','aa','African American, n (%)'),('CONT','BMI','BMI, kg/m2',1),
 ('CAT','current_smoker','Current smoker, n (%)'),('CONT','ATS_PackYears','Pack-years',1),
 ('CONT','SGRQ_scoreTotal','SGRQ total score',1),('CAT','mmrc2','mMRC dyspnea >=2, n (%)'),
 ('CAT','cb','Chronic bronchitis, n (%)'),('CONT','distwalked_m','6MWD, m',0),
 ('SECTION','SPIROMETRY'),('CONT','FEV1pp_post','FEV1 % predicted',1),('CONT','FEV1_FVC_post','FEV1/FVC',2),
 ('CONT','FVCpp_post','FVC % predicted',1),
 ('MULTI','gold_cat','GOLD stage at Phase 2, n (%)',['GOLD 0','GOLD 1','GOLD 2','GOLD 3-4']),
 ('SECTION','QUANTITATIVE CT'),('CONT','pctEmph_Thirona','Emphysema % (LAA<-950 HU)',1),
 ('CONT','WallAreaPct_seg_Thirona','Segmental wall area %',1),('CONT','Pi10_Thirona','Pi10 airway wall, mm',2),
 ('SECTION','BLOOD CELL COUNTS'),('CONT','wbc','White blood cells, x10^9/L',1),
 ('CONT','neutrophl_pct','Neutrophils %',1),('CONT','lymphcyt_pct','Lymphocytes %',1),
 ('CONT','monocyt_pct','Monocytes %',1),('CONT','eosinphl_pct','Eosinophils %',1)]

fmt_p=lambda p,s=False:("" if not np.isfinite(p) else ("<0.001" if p<0.001 else f"{p:.3f}")+("*" if s else ""))
def iqr(v,dp):
    f=lambda x:f"{x:,.{dp}f}"
    return f"{f(v.median())} ({f(v.quantile(.25))}-{f(v.quantile(.75))})"

gr={g:d[d.group==g] for g in G}
hdr=["Characteristic"]+[f"{g} (n={len(gr[g]):,})" for g in G]+["P value"]
rows=[]
for sp in SPEC:
    k=sp[0]
    if k=='SECTION': rows.append([sp[1]]+[""]*(len(G)+1))
    elif k=='CONT':
        _,col,lab,dp=sp; cells=[];sm=[]
        for g in G:
            v=pd.to_numeric(gr[g][col],errors='coerce').dropna(); sm.append(v.values)
            cells.append(iqr(v,dp) if len(v) else "--")
        sm=[x for x in sm if len(x)>1]
        rows.append([lab]+cells+[fmt_p(stats.kruskal(*sm).pvalue if len(sm)>=2 else np.nan)])
    elif k=='CAT':
        _,col,lab=sp; cells=[];tab=[]
        for g in G:
            v=pd.to_numeric(gr[g][col],errors='coerce').dropna(); n1=int(v.sum())
            cells.append(f"{n1:,} ({100*n1/len(v):.1f})" if len(v) else "--"); tab.append([n1,len(v)-n1])
        tab=np.array(tab,float)
        rows.append([lab]+cells+[fmt_p(stats.chi2_contingency(tab[tab.sum(1)>0])[1])])
    elif k=='MULTI':
        _,col,lab,lv=sp
        tab=np.array([[int((gr[g][col]==l).sum()) for l in lv] for g in G],float)
        rows.append([lab]+[""]*len(G)+[fmt_p(stats.chi2_contingency(tab[:,tab.sum(0)>0])[1],True)])
        for j,l in enumerate(lv):
            tot=[len(gr[g][gr[g][col].notna()]) for g in G]
            rows.append([l]+[f"{int(tab[i,j]):,} ({100*tab[i,j]/tot[i]:.1f})" for i in range(len(G))]+[""])

out=OUT_FILE
with open(out,'w',newline='') as f:
    w=csv.writer(f); w.writerow(hdr); w.writerows(rows)
print('\n=== PHASE 2 CLASSIFIED TABLE 1 (quantitative CT schema) ===')
print(' | '.join(hdr))
for r in rows: print(' | '.join(r))
print('\nwrote', out)
