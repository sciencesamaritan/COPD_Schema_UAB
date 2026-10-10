"""Load the COPDGene longitudinal phenotype file (SAS) and keep the variables used for the schema and Table 1."""
import pandas as pd, os
KEEP = ['sid','Phase_study','MultiDim_COPD_Subtype_P1','MultiDim_COPD_score_P1',
 'age_visit','gender','race','BMI','smoking_status','ATS_PackYears',
 'FEV1pp_post','FEV1_FVC_post','FVCpp_post','SGRQ_scoreTotal','distwalked',
 'MMRCDyspneaScor','pctEmph_Thirona','Pi10_Thirona','WallAreaPct_seg_Thirona',
 'Chronic_Bronchitis','finalgold_visit','wbc','neutrophl_pct','lymphcyt_pct',
 'monocyt_pct','eosinphl_pct','basophl_pct',
 'DL_CLE','DL_PSE','DL_BWT','CT_Visual_Emph_Severity','CT_Visual_Wall_Thickening',
 'CT_Visual_Emph_Paraseptal','CT_Visual_Adjudicated','FEV1_FVC_pre','Height_CM','Weight_KG']
#### Edit paths to match your file system ####
PHENO_SAS = 'path/to/COPDGn_P1P2P3P4_SM_NS_Long_Nov25.sas7bdat'
CACHE = 'copdgene_pheno_slim.pkl'
def load():
    if os.path.exists(CACHE): return pd.read_pickle(CACHE)
    rd = pd.read_sas(PHENO_SAS, format='sas7bdat', chunksize=20000, encoding='latin-1')
    df = pd.concat([ch[[c for c in KEEP if c in ch.columns]] for ch in rd], ignore_index=True)
    df['sid']=df['sid'].astype(str).str.strip()
    df.to_pickle(CACHE)
    return df
