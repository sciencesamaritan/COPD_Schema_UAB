import re
THEMES = [
 ("Cytokine / NF-κB", ["INTERLEUKIN", "_IL1", "_IL6", "_IL2", "IL10", "IL12", "IL17", "TNF", "NFKB", "CYTOKINE", "CHEMOKINE", "INFLAMM", "JAK_STAT"]),
 ("T / B cell receptor", ["TCR", "CD28", "T_CELL", "BCR", "B_CELL", "ZAP70", "LCK", "CTLA4", "PD1", "ANTIGEN", "MHC", "IMMUNOREGULATORY", "COSTIMULATION"]),
 ("Interferon", ["INTERFERON", "IFN", "ANTIVIRAL", "ISG15"]),
 ("Innate / neutrophil", ["NEUTROPHIL", "TOLL", "TLR", "FCER", "FCGR", "LECTIN", "COMPLEMENT", "DEGRANULATION", "INNATE", "MYD88", "INFLAMMASOME", "NOD", "DECTIN"]),
 ("GPCR / Ca2+ / NFAT", ["GPCR", "G_PROTEIN", "CA2", "CALCIUM", "PAR1", "EDG", "LPAR", "GNB", "GNA", "NFAT", "PLCG", "ITPR", "ETS_PATHWAY", "PYK2"]),
 ("MAPK / Ras / ERK", ["MAPK", "_ERK", "RAS", "RAF", "MEK", "MAP_KINASE", "P38", "JNK"]),
 ("Growth-factor RTK / PI3K", ["FGFR", "EGFR", "PDGF", "_KIT", "ERBB", "INSULIN", "PI3K", "_AKT", "MTOR", "IGF", "VEGF", "PIP"]),
 ("Proteasome / ubiquitin", ["PROTEASOME", "UBIQUITIN", "DEGRADATION", "_SCF", "APC_C", "NEDDYLATION", "DEUBIQUITINATION"]),
 ("DNA repair / TP53", ["TP53", "DNA_REPAIR", "HDR", "HOMOLOGOUS", "MISMATCH", "EXCISION", "FANCONI", "DOUBLE_STRAND", "BASE_EXCISION", "TRANSLESION", "GAP_FILLING"]),
 ("Cell cycle / mitosis", ["MITOTIC", "CELL_CYCLE", "_G1", "_G2", "KINETOCHORE", "PROMETAPHASE", "CHECKPOINT", "MITOSIS", "CENTROSOME", "SEPARATION", "DNA_REPLICATION", "S_PHASE"]),
 ("Chromatin / methylation", ["METHYLATION", "HISTONE", "CHROMATIN", "SUMOYLATION", "ACETYLATION", "HDAC", "PRC2", "EPIGENETIC", "SIRT"]),
 ("Translation / ribosome", ["TRANSLATION", "RIBOSOM", "ELONGATION", "_EIF", "NONSENSE", "PEPTIDE_CHAIN", "SRP_DEPENDENT", "RRNA", "SELENO"]),
 ("Glycosylation", ["GLYCOSYL", "GLYCAN", "MUCIN", "HYPUSIN", "CARBOXYLATION", "GPI"]),
 ("Mitochondria / OXPHOS", ["RESPIRATORY", "ELECTRON", "COMPLEX_I", "OXIDATIVE_PHOSPHORYLATION", "MITOCHONDRIAL", "TCA", "CITRIC"]),
 ("Metabolism", ["METABOLISM", "GLYCOLYSIS", "LIPID", "FATTY", "CHOLESTEROL", "AMINO_ACID", "GLUCOKINASE", "BIOSYNTHESIS", "PURINE", "PYRIMIDINE"]),
 ("Trafficking / vesicles", ["TRANSPORT", "VESICLE", "GOLGI", "ENDOSOM", "CLATHRIN", "ENDOCYTOSIS", "LYSOSOM", "ER_TO_GOLGI", "COPI", "COPII"]),
 ("RNA processing", ["SPLICING", "MRNA", "TRNA", "RNA_POLYMERASE", "SNRNP", "CAPPING"]),
 ("Development / neuronal", ["NEURON", "AXON", "NOTCH", "WNT", "HEDGEHOG", "NEUROTRANSMITTER", "SYNAP", "PROGENITOR"]),
 ("TGF-β / ECM", ["TGF", "TGFB", "SMAD", "COLLAGEN", "EXTRACELLULAR_MATRIX", "INTEGRIN", "TGFBR3"]),
 ("Apoptosis", ["APOPTO", "CASPASE", "PROGRAMMED_CELL_DEATH", "NECROPTOSIS"]),
]
def theme_of(name):
    u = "_" + name.upper(); return [t for t, ks in THEMES if any(k in u for k in ks)]
def label_program(names):
    from collections import Counter
    c = Counter(t for nm in names for t in theme_of(nm)); n = len(names)
    if not c: return "Mixed"
    (a, ca), *rest = c.most_common(3)
    if rest and rest[0][1] >= 0.6 * ca: return f"{a} + {rest[0][0]}"
    return a

def label_program_enriched(names, universe):
    """Theme label by enrichment over the tested universe (log fold of member share vs background share),
    damped when a theme covers < 15% of members. Adds a second theme when it scores >= 75% of the top."""
    import numpy as np
    from collections import Counter
    bg = Counter(t for p in universe for t in theme_of(p)); N = len(universe); n = len(names)
    c = Counter(t for p in names for t in theme_of(p))
    sc = sorted(((np.log(((c[t] + .5) / n) / (bg[t] / N)) * min(1, c[t] / (0.15 * n)), t) for t in c), reverse=True)
    if not sc or sc[0][0] <= 0: return "Mixed"
    if len(sc) > 1 and sc[1][0] >= 0.75 * sc[0][0]: return f"{sc[0][1]} + {sc[1][1]}"
    return sc[0][1]
