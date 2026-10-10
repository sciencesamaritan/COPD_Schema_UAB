"""Multidimensional COPD diagnostic schema (Bhatt et al., JAMA 2025) computed from phenotype variables.

schema(ph, mode): 'visual' = CT visual emphysema >= mild and definite bronchial wall thickening;
                  'quant'  = emphysema >= 5% (LAA < -950 HU) and Pi10 >= 2.24 mm.
Run as a script to validate against the published Phase 1 counts, evaluate the schema at Phase 2 in the
RNA-seq sample, and tabulate Phase 1 -> Phase 2 migration (Supplementary Table 7).
"""
import pandas as pd, numpy as np
from copdgene_phenotypes import load
#### Edit paths to match your file system ####
RNA_MASTER = 'path/to/master.file.freeze5.phase2.csv'


def schema(ph, mode):
    """mode 'visual' = Fleischner CLE>=mild & definite BWT; 'quant' = %emph>=5 & Pi10>=2.24"""
    if mode=='visual':
        emph = ph['CT_Visual_Emph_Severity']>=2; bwt = ph['CT_Visual_Wall_Thickening']>=2
        img  = ['CT_Visual_Emph_Severity','CT_Visual_Wall_Thickening']
    else:
        emph = ph['pctEmph_Thirona']>=5.0;       bwt = ph['Pi10_Thirona']>=2.24
        img  = ['pctEmph_Thirona','Pi10_Thirona']
    sd=ph['MMRCDyspneaScor']>=2; sq=ph['SGRQ_scoreTotal']>=25; sb=ph['Chronic_Bronchitis']==1
    nm = emph.astype(int)+bwt.astype(int)+sd.astype(int)+sq.astype(int)+sb.astype(int)
    AO = ph['FEV1_FVC_post']<0.70
    comp = ph[img+['MMRCDyspneaScor','SGRQ_scoreTotal','Chronic_Bronchitis','FEV1_FVC_post']].notna().all(axis=1)
    g = pd.Series(np.where(AO&(nm>=1),'Established COPD',
                  np.where((~AO)&(nm>=3),'Newly diagnosed',
                  np.where(AO,'Potentially excluded','No COPD'))), index=ph.index)
    return g.where(comp)

G=['No COPD','Newly diagnosed','Potentially excluded','Established COPD']
def show(g,label):
    v=g.value_counts().reindex(G).fillna(0).astype(int)
    print(f'{label:52s} n={int(g.notna().sum()):5d}  '+' / '.join(f'{v[k]:,}' for k in G))

if __name__ == '__main__':
    df = load()
    p1=df[df.Phase_study==1].drop_duplicates('sid').set_index('sid')
    print('=== VALIDATION at Phase 1, full COPDGene ===')
    print('eTable E9 target (quant): None 4,779 | Major 3,859 | Minor 776 ; potentially excluded 305')
    gq=schema(p1,'quant');  show(gq, 'Phase 1 quantitative (Pi10>=2.24, emph>=5%)')
    print(f'   -> None = NoCOPD+PotExcl = {int((gq=="No COPD").sum()+(gq=="Potentially excluded").sum()):,}')
    gv=schema(p1,'visual'); show(gv, 'Phase 1 visual (CLE>=mild, definite BWT)')
    print('   -> reference: 4,721 None (4,439+282) / 811 minor / 3,884 major')

    print()
    print('=== PHASE 2, RNA-seq analytic sample ===')
    master=pd.read_csv(RNA_MASTER)
    rna=set(master.loc[master['final.analyzed']==1,'aliasID'].astype(str).str.strip())
    p2=df[df.Phase_study==2].drop_duplicates('sid').set_index('sid')
    diffs=['neutrophl_pct','lymphcyt_pct','monocyt_pct','eosinphl_pct']

    ids=[s for s in p2.index if s in rna]
    b=p2.loc[ids]
    print(f'visit-2 RNA-seq records linked to phenotype: {len(b):,}')

    g2 = schema(b,'quant')
    covok = b['smoking_status'].notna() & b[diffs].notna().all(axis=1) & b['age_visit'].notna()
    keep = g2.notna() & covok
    print(f'complete on Phase 2 schema + covariates:     {int(keep.sum()):,}')
    show(g2[keep], 'PHASE 2 quantitative schema')

    g1q = schema(p1,'quant').reindex(b.index)
    g1v = schema(p1,'visual').reindex(b.index)
    show(g1q[keep], 'Phase 1 quantitative (same participants)')
    show(g1v[keep & g1v.notna()], 'Phase 1 visual (same participants)')

    print()
    print('Migration Phase 1 -> Phase 2 (quantitative schema, same participants):')
    m=pd.crosstab(g1q[keep].rename('Phase 1'), g2[keep].rename('Phase 2')).reindex(index=G,columns=G).fillna(0).astype(int)
    print(m.to_string())
    print(f'\nstable: {np.trace(m.values):,} of {m.values.sum():,} ({100*np.trace(m.values)/m.values.sum():.1f}%)')
