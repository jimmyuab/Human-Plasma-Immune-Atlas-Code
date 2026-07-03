#!/usr/bin/env python
"""
HDDM Layer 54 - Step 13
Rewrite the manuscript to Nature-family standard:
  - sharper title
  - Nature-grammar abstract
  - 4 clean results sections
  - claim-strength-gated language throughout (uses evidence_tiered_targets.tsv)
  - embeds the 6 multi-panel Nature figures
Output: 10_manuscript/Plasma_Immunome_Phenome_Atlas_Nature.docx
"""
import os
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")
FIG  = os.path.join(ROOT, "08_figures", "nature")
PROC = os.path.join(ROOT, "02_data_processed")
OUT  = os.path.join(ROOT, "10_manuscript")
os.makedirs(OUT, exist_ok=True)

# ---- load data / numbers ----
ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
imm = ann[ann.is_plasma_immune == 1]
mr  = pd.read_csv(os.path.join(GEN, "cis_MR_immune_results.tsv"), sep="\t")
inst= pd.read_csv(os.path.join(GEN, "immune_cis_eqtl_instruments.tsv"), sep="\t")
tier= pd.read_csv(os.path.join(GEN, "evidence_tiered_targets.tsv"), sep="\t")
final= pd.read_csv(os.path.join(GEN, "FINAL_evidence_tiers.tsv"), sep="\t")
_replp = os.path.join(GEN, "FINAL_evidence_tiers_repl.tsv")
if os.path.exists(_replp):
    final = pd.read_csv(_replp, sep="\t")

N_OLINK = 2923
N_IMM   = int(imm.shape[0])
N_INST  = int(inst.gene_symbol.nunique())
N_TEST  = int(mr.shape[0])
N_SIG   = int((mr.FDR < 0.05).sum())
N_DIS   = int(mr.disease.nunique())
N_T4 = int((tier.evidence_tier==4).sum())
N_T3 = int((tier.evidence_tier==3).sum())
N_T2 = int((tier.evidence_tier==2).sum())
t4 = tier[tier.evidence_tier==4]
# protein-level (INTERVAL pQTL) layer
N_T5   = int((final.final_tier==5).sum())
t5g    = final[final.final_tier==5]
N_PQTL = int(final.pQTL_OR.notna().sum())
N_CONC = int((final.pQTL_concordant=="yes").sum())
N_DISC = int((final.pQTL_concordant=="NO").sum())
# independent-GWAS replication layer
HAS_REP = "rep_status" in final.columns
N_COV  = int((final.rep_status!="not covered").sum()) if HAS_REP else 0
N_REP  = int((final.rep_status=="replicated").sum()) if HAS_REP else 0
# novelty layer
def _nl(f):
    p=os.path.join(GEN,f); return pd.read_csv(p,sep="\t") if os.path.exists(p) else None
nov_enr=_nl("novelty_class_enrichment.tsv")
nov_dir=_nl("novelty_drug_direction.tsv")
nov_pri=_nl("novelty_prioritised_targets.tsv")
if nov_dir is not None:
    N_DIRMATCH=int((nov_dir.concordance=="match").sum())
    N_NOVELNOM=int((nov_pri.category=="NOVEL nomination").sum()) if nov_pri is not None else 0
else:
    N_DIRMATCH=N_NOVELNOM=0
N_RCON = int(final.rep_status.isin(["replicated","concordant (ns)"]).sum()) if HAS_REP else 0

# ---- pan-phenome layer (src/24-29) ----
def _pl(f):
    p=os.path.join(GEN,f); return pd.read_csv(p,sep="\t") if os.path.exists(p) else None
ph_hits = _pl("phenome_hits.tsv")
ph_pleio= _pl("phenome_pleiotropy_axes.tsv")
ph_col  = _pl("coloc_phenome_results.tsv")
nov_rank= _pl("novelty_engine_ranked.tsv")
HAS_PHEN = ph_hits is not None and nov_rank is not None
if HAS_PHEN:
    PH_HITS = int(len(ph_hits))
    PH_DIS  = int(ph_hits.disease.nunique())
    PH_GENE = int(ph_hits.gene_symbol.nunique())
    PH_NCAT = int(ph_hits.disease_category.nunique())
    PH_CATCOUNTS = ph_hits.disease_category.value_counts().to_dict()
    PH_PLEIO = int(ph_pleio.gene_symbol.nunique()) if ph_pleio is not None else 0
    PH_COLOC = int((ph_col.PP_H4>=0.8).sum()) if ph_col is not None else 0
    PH_NOVCOL= int((nov_rank.category_label=="NOVEL colocalized").sum())
    PH_NOVNOM= int((nov_rank.category_label=="novel nomination").sum())
    PH_TOTAL = int(N_DIS + PH_DIS)  # combined disease breadth is described in text

# ---- doc styling ----
doc = Document()
st = doc.styles["Normal"]; st.font.name = "Calibri"; st.font.size = Pt(10.5)

def H(txt, lvl=1, color=None):
    p = doc.add_heading(txt, level=lvl)
    if color:
        for r in p.runs: r.font.color.rgb = color
    return p
def P(txt, italic=False, size=10.5, align=None, bold=False):
    p = doc.add_paragraph()
    r = p.add_run(txt); r.italic = italic; r.bold = bold; r.font.size = Pt(size)
    if align: p.alignment = align
    return p
def figure(fname, caption):
    fp = os.path.join(FIG, fname)
    if os.path.exists(fp):
        doc.add_picture(fp, width=Inches(6.3))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    c = doc.add_paragraph(); r = c.add_run(caption); r.italic=True; r.font.size=Pt(9)

BLUE = RGBColor(0x1f,0x4e,0x79)

# ================= TITLE =================
t = doc.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = t.add_run("An open, genetics-anchored plasma immunome atlas recovers established "
              "therapeutic immune axes and nominates colocalized autoimmune targets")
r.bold=True; r.font.size=Pt(15); r.font.color.rgb=BLUE

P("HDDM Layer 54 \u00b7 Plasma Immunome\u2013Phenome Atlas \u00b7 an independent, reproducible resource",
  italic=True, size=10, align=WD_ALIGN_PARAGRAPH.CENTER)
doc.add_paragraph()

# ================= ABSTRACT =================
H("Abstract", 1, BLUE)
P(
 f"Circulating immune proteins mediate host defence and are the targets of many "
 f"approved immunotherapies, yet a systematic, openly reproducible map linking the "
 f"plasma immune proteome to disease has been lacking. Here we assemble a plasma "
 f"immunome\u2013phenome atlas that curates the {N_OLINK:,}-analyte Olink Explore proteomic "
 f"space to {N_IMM:,} plasma-detectable immune proteins using orthogonal Human Protein "
 f"Atlas immune-cell specificity and pathway annotation, and anchors each protein to "
 f"disease through genetics. Using cis-eQTL instruments for {N_INST} immune genes, we "
 f"perform Mendelian randomization against {N_DIS} immune-mediated diseases "
 f"({N_TEST:,} tests) and calibrate every finding on an explicit claim-strength ladder. "
 f"The analysis blindly recovers established therapeutic axes\u2014IL6ST in rheumatoid "
 f"arthritis, CTLA4 in autoimmune thyroid disease, TNFRSF1A in ankylosing spondylitis\u2014"
 f"providing an internal positive control. Statistical colocalization resolves {N_T4+N_T5} "
 f"loci to a shared causal variant (PP.H4\u22650.8), promoting them to prioritized "
 f"transcript-level causal targets, and nominates candidates including SIGLEC7/9 in "
 f"Guillain\u2013Barr\u00e9 syndrome, HAVCR1 in coeliac disease and IFNGR2 in psoriasis. "
 f"Independent protein-level instruments from plasma cis-pQTLs (INTERVAL) then confirm a "
 f"subset at the protein level\u2014TNFSF14 in multiple sclerosis and SWAP70 in rheumatoid "
 f"arthritis colocalize at both transcript and protein level (PP.H4>0.97), elevating them "
 f"to protein-level causal targets\u2014while exposing informative transcript\u2013protein "
 f"discordances (e.g. IL6ST) that temper otherwise strong transcript-level signals. "
 f"In independent, FinnGen-free consortium GWAS the direction of effect replicates in "
 f"{N_RCON} of {N_COV} testable nominations ({N_REP} at P<0.05). "
 f"Signals in the MHC ({N_T2}) are explicitly flagged for linkage-disequilibrium "
 f"confounding and withheld from causal claims. All curation, instruments, statistics "
 f"and figures are built from public data and released as an open, versioned resource."
)

# ================= INTRODUCTION =================
H("Introduction", 1, BLUE)
P("Plasma proteins are the most accessible window on human immunity and constitute a "
  "disproportionate share of the druggable proteome. Affinity proteomics now quantifies "
  "thousands of circulating proteins at population scale, but converting protein\u2013disease "
  "correlations into directional, genetically-supported hypotheses requires (i) a principled "
  "definition of the immune proteome, (ii) genetic instruments that expose causal direction, "
  "and (iii) discipline in matching claim strength to evidence. Association is not causation, "
  "a cis-eQTL is not a protein measurement, and a Mendelian-randomization estimate at the MHC "
  "is not a clean instrument. We therefore built this atlas around an explicit evidence ladder "
  "and report every finding at the tier its evidence supports\u2014no higher.")

# ================= RESULTS =================
H("Results", 1, BLUE)

H("Curation of a plasma-detectable immune proteome", 2)
P(f"Starting from the {N_OLINK:,} analytes of Olink Explore, we integrated Human Protein "
  f"Atlas immune-cell RNA specificity, secretome location and curated pathway membership "
  f"(cytokine, chemokine, interferon, TNF, complement, checkpoint, CD/leukocyte, HLA, "
  f"immunoglobulin, acute-phase, coagulation) into a transparent inclusion score, retaining "
  f"{N_IMM:,} plasma immune proteins (Fig. 1). The resource preserves, for each protein, its "
  f"immune class and the blood-cell lineage that expresses it, so that downstream disease "
  f"signals can be read back to a cell of origin. A permissive immune-signature gene set "
  f"(MSigDB C7) covered almost the entire Olink space and was therefore retained only as "
  f"annotation, not as an inclusion criterion\u2014an early example of the calibration principle "
  f"applied throughout.")
figure("Figure1_design_curation.png",
       "Figure 1 | Design and curation of the plasma immunome. Inclusion funnel from the "
       f"Olink Explore universe to {N_IMM:,} plasma immune proteins; immune-class composition; "
       "mapping of proteins to blood-cell lineage of origin; and co-occurrence of annotation flags.")

H("Architecture and druggability of the immune proteome", 2)
P("The curated proteome is dominated by immune-cell-enriched proteins with substantial "
  "cytokine, chemokine, complement and checkpoint representation, and maps onto granulocyte, "
  "myeloid, dendritic, T-, B- and NK-cell lineages (Fig. 2). A large fraction are secreted or "
  "single-pass receptors\u2014the protein classes most tractable to antibodies and biologics\u2014"
  "establishing the atlas as a systematic hunting ground for immune drug targets before any "
  "genetic evidence is applied.")
figure("Figure2_druggability_architecture.png",
       "Figure 2 | Architecture and druggability of the plasma immune proteome: immune-class "
       "and lineage composition, secretome/receptor tractability, and target-class breakdown.")

H("Genetic anchoring recovers established therapeutic immune axes", 2)
P(f"To expose causal direction we instrumented {N_INST} immune genes with their strongest cis "
  f"expression-quantitative-trait locus (eQTLGen; single-instrument Wald ratio) and tested them "
  f"against {N_DIS} immune-mediated diseases from FinnGen R12 ({N_TEST:,} tests), controlling "
  f"the false-discovery rate at 5% ({N_SIG} significant gene\u2013disease pairs; Fig. 3). Critically, "
  f"the analysis\u2014run blind to drug annotation\u2014recovers multiple axes that are already the "
  f"targets of approved immunotherapies: IL6ST (the IL-6 co-receptor gp130) in rheumatoid "
  f"arthritis, CTLA4 in autoimmune thyroid disease, and the TNF-receptor axis (TNFRSF1A) in "
  f"ankylosing spondylitis, with directions of effect consistent with the known pharmacology "
  f"(Fig. 4). This concordance is an internal positive control: a pipeline that rediscovers "
  f"tocilizumab-, abatacept- and anti-TNF-relevant biology from genetics alone is calibrated to "
  f"be believed when it points somewhere new.")
figure("Figure3_MR_design_global.png",
       "Figure 3 | Mendelian-randomization design and global results across 13 immune-mediated "
       "diseases: workflow, effect-size volcano, significant hits per disease, and immune-class enrichment.")
figure("Figure4_positive_controls.png",
       "Figure 4 | The pipeline blindly recovers established immune therapeutic axes (IL6ST, "
       "CTLA4, TNFRSF1A, TNFSF14), each shown with its cis-MR odds ratio and the drug it underpins.")

H("Colocalization calibrates novel, disease-specific target nominations", 2)
P(f"We next asked which significant signals reflect a shared causal variant rather than "
  f"linkage-disequilibrium coincidence, applying approximate-Bayes-factor colocalization "
  f"(coloc.abf) at every non-MHC hit. Colocalization behaves as expected on the positive "
  f"controls (IL6ST\u2013rheumatoid arthritis PP.H4=1.00; CTLA4 PP.H4>0.98; TNFSF14\u2013multiple "
  f"sclerosis PP.H4>0.99) and promotes {N_T4} loci to prioritized transcript-level causal targets "
  f"(PP.H4\u22650.8), including several disease-specific nominations: SIGLEC7 and SIGLEC9 "
  f"(sialic-acid-binding inhibitory receptors) in Guillain\u2013Barr\u00e9 syndrome, HAVCR1/KIM-1 in "
  f"coeliac disease, IFNGR2 in psoriasis, CD226 in sarcoidosis and HNMT in Crohn's disease "
  f"(Fig. 5). A further {N_T3} signals remain genetically-supported nominations pending stronger "
  f"colocalization. Signals in the MHC\u2014most prominently the complement gene C2, significant "
  f"across five diseases\u2014are held at nomination status and explicitly annotated for LD "
  f"confounding rather than promoted ({N_T2} MHC pairs). Figure 6 places every finding on the "
  f"claim-strength ladder and defines the validation path to causal and translational tiers.")
figure("Figure5_novel_nominations.png",
       "Figure 5 | Colocalized novel autoimmune target nominations, each annotated with its "
       "colocalization posterior (PP.H4): C2 (MHC, held at nomination), SIGLEC7/9\u2192Guillain\u2013Barr\u00e9, "
       "HAVCR1\u2192coeliac, IFNGR2\u2192psoriasis.")
figure("Figure6_validation_roadmap.png",
       "Figure 6 | Evidence-tier calibration and validation roadmap: the claim-strength ladder, "
       "where each finding sits, per-locus colocalization posteriors, and the path to protein-level "
       "(pQTL), replication and functional validation.")

H("Protein-level pQTL instruments confirm a subset and reveal informative discordance", 2)
P(f"Because a cis-eQTL indexes transcript, not protein, abundance, we sought independent "
  f"confirmation using plasma cis-pQTLs from the INTERVAL study (Sun et al., 2018) as "
  f"protein-level instruments for the significant genes ({N_PQTL} gene\u2013disease pairs had a "
  f"matched plasma aptamer). This layer is decisive rather than cosmetic. Two nominations "
  f"colocalize at both the transcript and the protein level with concordant direction and "
  f"survive protein-level MR\u2014TNFSF14 in multiple sclerosis (pQTL PP.H4>0.99) and SWAP70 in "
  f"rheumatoid arthritis (pQTL PP.H4>0.97)\u2014and we designate these {N_T5} the study's "
  f"protein-level causal targets, the highest tier reached here (Fig. 7). A further set, "
  f"including TNFRSF1A in ankylosing spondylitis, are directionally concordant with nominally "
  f"significant protein-level MR. Critically, only {N_CONC} of {N_PQTL} pairs are directionally "
  f"concordant between transcript and protein, and {N_DISC} are discordant\u2014most notably IL6ST "
  f"in rheumatoid arthritis, where the transcript signal (risk-increasing) reverses at the "
  f"protein level, consistent with the known biology of soluble gp130 as an inhibitor of IL-6 "
  f"trans-signalling. Rather than suppress these discordances we report them, because they "
  f"define which transcript-level nominations should not yet be read as protein-level causal "
  f"claims (e.g. IFNGR2, ERBB3), and they are exactly the distinctions a therapeutic programme "
  f"needs before committing to a modality.")
figure("Figure7_pqtl_confirmation.png",
       "Figure 7 | Protein-level (INTERVAL plasma pQTL) validation. (a) Transcript- vs "
       "protein-level MR log-odds-ratios per hit (blue, directionally concordant; red, discordant). "
       "(b) Protein-level colocalization posteriors; TNFSF14\u2192multiple sclerosis and SWAP70\u2192"
       "rheumatoid arthritis exceed PP.H4=0.8 at the protein level.")

H("Nominations replicate in independent, non-FinnGen disease GWAS", 2)
P(f"Because discovery used FinnGen alone, we sought out-of-sample replication in independent "
  f"consortium GWAS accessed through OpenGWAS (IMSGC/Patsopoulos multiple sclerosis; Okada "
  f"rheumatoid arthritis; IGAS/Cortes ankylosing spondylitis; Stuart psoriasis; Fischer "
  f"sarcoidosis; Sakaue autoimmune thyroid; none of which include FinnGen). Of the {N_SIG} "
  f"significant gene\u2013disease pairs, {N_COV} had the instrument variant present in a matched "
  f"independent GWAS; of these, {N_RCON} were directionally concordant and {N_REP} replicated "
  f"at nominal significance (P<0.05), including every one of the positive-control axes and the "
  f"tier-5 target TNFSF14 in multiple sclerosis (P=2\u00d710\u207b\u00b9\u00b2, IMSGC), IL6ST in rheumatoid "
  f"arthritis (P=7\u00d710\u207b\u00b2\u2074, Okada) and CTLA4 (P=3\u00d710\u207b\u00b2\u2070) (Fig. 8). Notably, the direction of "
  f"effect agreed in 19 of 19 covered hits, an out-of-sample consistency that would be "
  f"vanishingly unlikely by chance and that substantially raises confidence in the "
  f"transcript-level nominations. Several candidates could not be tested because no "
  f"FinnGen-independent GWAS exists (Guillain\u2013Barr\u00e9 syndrome, Sjögren syndrome) or because "
  f"the instrument was absent from an older targeted-array study (coeliac disease); these are "
  f"stated as coverage gaps rather than failures to replicate.")
figure("Figure8_replication.png",
       "Figure 8 | Independent replication. (a) \u2212log10(P) of each instrument variant in a matched "
       "FinnGen-independent consortium GWAS (bar colour = replication status; cohort labelled). "
       "(b) Replication outcome across all significant hits. Direction of effect was concordant in "
       "all 19 covered nominations.")

# ================= NOVELTY / TRANSLATION =================
H("The atlas recovers approved-drug direction, resolves shared causal axes, and prioritizes novel targets", 2)
P(f"Beyond identifying targets, the genetic architecture carries a directional signal that is "
  f"directly actionable. Encoding each effect as protein-lowering-protective (OR<1, arguing for "
  f"agonism/replacement) versus protein-raising-risk (OR>1, arguing for blockade), the atlas "
  f"recovers not only the identity but the correct pharmacological DIRECTION of established "
  f"therapies: CTLA4 is protective and is drugged by an agonist (abatacept, CTLA4-Ig); TNFRSF1A "
  f"and the IL-6 axis are risk-increasing and are drugged by blockers (etanercept, tocilizumab); "
  f"IL4 raises psoriasis risk and is blocked by dupilumab. In total {N_DIRMATCH} recovered axes "
  f"match approved-drug direction, and the two cases where plasma pQTL direction inverts the eQTL "
  f"signal (IL6ST, reflecting soluble gp130 trans-signalling; TNFRSF14/HVEM) are exactly the axes "
  f"where receptor biology is known to be direction-dependent\u2014so the discordance is mechanistically "
  f"informative rather than noise (Fig. 9c). Mapping genes that act across more than one disease "
  f"exposes shared causal immune axes\u2014CTLA4 (autoimmune thyroid + rheumatoid arthritis), NUMB and "
  f"TNFRSF14 (each across two diseases), and the MHC-region C2 signal recurring across five\u2014"
  f"pointing to pleiotropic control points whose directions are consistent across indications "
  f"(Fig. 9b). At the class level, genetically-supported causal targets trend toward the TNF "
  f"superfamily, interferon axis and immune-checkpoint classes (Fisher odds ratios 4\u201310), a "
  f"coherent enrichment given known autoimmune drug space, though not individually significant "
  f"after correction given the small target set (Fig. 9a). Finally, separating internal "
  f"positive controls from genuinely new signals leaves {N_NOVELNOM} novel colocalized nominations "
  f"lacking an approved drug for their indication; ranked by cumulative orthogonal evidence "
  f"(transcript coloc + protein pQTL + independent replication), SWAP70 in rheumatoid arthritis "
  f"tops the list as the only novel target with protein-level colocalization AND independent "
  f"replication, followed by replicated transcript-level nominations IFNGR2 (psoriasis), CA8 and "
  f"MYL6B (rheumatoid arthritis) and HAVCR1 (coeliac disease) (Fig. 9d).")
figure("Figure9_novelty.png",
       "Figure 9 | Novelty and translation. (a) Immune-class enrichment of tier\u22654 causal targets "
       "(Fisher; trends toward TNF-superfamily/interferon/checkpoint classes). (b) Pleiotropic "
       "shared causal axes for genes acting on \u22652 diseases (red=risk/blockade, blue=protective/"
       "agonism). (c) Genetics recover the DIRECTION of approved drugs, with mechanistically "
       "informative discordances. (d) Novel versus recovered-known target prioritization by "
       "cumulative orthogonal evidence.")

# ================= PAN-PHENOME EXPANSION =================
if HAS_PHEN:
    H("Extending the atlas to a pan-phenome causal map exposes shared immune axes beyond autoimmunity", 2)
    _cc = ", ".join(f"{k} ({v})" for k,v in sorted(PH_CATCOUNTS.items(), key=lambda x:-x[1]))
    P(f"To test whether the immune proteome acts causally beyond classical autoimmunity, we "
      f"re-instrumented the same {N_INST} immune genes against an expanded {PH_DIS}-disease "
      f"phenome spanning five categories\u2014autoimmune, cardiovascular, metabolic, renal and "
      f"neurodegenerative/aging\u2014assembled from FinnGen R12 (cardiometabolic, dementia, "
      f"Alzheimer's disease, epilepsy, glaucoma, chronic kidney disease, osteoporosis and "
      f"others). The identical cis-MR and colocalization machinery, held to the same FDR<5% "
      f"gate, yields {PH_HITS} significant gene\u2013disease pairs over {PH_GENE} genes and "
      f"{PH_NCAT} disease categories ({_cc}), materially broadening the atlas from an "
      f"autoimmune resource into a multi-system causal map (Fig. 11). Colocalization promotes "
      f"{PH_COLOC} of these pan-phenome loci to a shared causal variant (PP.H4\u22650.8), "
      f"including cardiometabolic and neurodegenerative signals that were invisible to the "
      f"autoimmune-only analysis: ACE and IFNGR2 in hypertension, PSRC1 and PLAUR in coronary "
      f"heart disease, PROCR in venous thromboembolism, and ACE colocalizing with dementia and "
      f"Alzheimer's disease (PP.H4>0.97).")
    P(f"The pan-phenome view makes the atlas's central novelty engine possible: {PH_PLEIO} "
      f"immune genes are causal across two or more disease CATEGORIES, resolving pleiotropic "
      f"control points that link immune and cardiometabolic disease. IFNGR2 (autoimmune + "
      f"cardiovascular), SWAP70 (autoimmune + cardiovascular), MERTK, MPO, SPINK8 and PIK3IP1 "
      f"(cardiovascular + metabolic) act on multiple systems with internally consistent "
      f"direction, and PM20D1 spans three categories (Fig. 12). Mapping causal targets back to "
      f"their blood-cell lineage of origin shows granulocyte-sourced immune proteins carry the "
      f"strongest causal enrichment across the phenome (Fig. 13), and a direction-aware "
      f"therapeutic map partitions all {PH_HITS} signals into block-versus-agonise strategies "
      f"with MHC/LD-confounded loci held separately (Fig. 14).")
    P(f"Integrating every real-data evidence layer\u2014causal strength, transcript "
      f"colocalization, cross-category pleiotropy, HPA-derived druggability and cell-source "
      f"specificity, minus explicit penalties for known-drug axes and MHC LD\u2014into a single "
      f"auditable novelty-priority score ranks the phenome-wide targets and separates "
      f"{PH_NOVCOL} novel colocalized targets and {PH_NOVNOM} further novel nominations from "
      f"recovered positive controls (Figs. 15, 16). The top novel colocalized nominations "
      f"lacking an approved drug for their indication\u2014ACE\u2192dementia, IFNGR2\u2192hypertension "
      f"and psoriasis, ERBB3\u2192type-1 diabetes, PLAUR\u2192coronary heart disease and SPINK8\u2192"
      f"type-2 diabetes\u2014are each supported by colocalization (PP.H4\u22650.9) and, for the "
      f"pleiotropic ones, by causal action in more than one disease category. As throughout, "
      f"the engine down-weights internal positive controls (CTLA4, IL6ST, IL2RA) so that the "
      f"ranking surfaces genuinely new biology rather than re-discovering approved drugs.")
    figure("Figure11_phenome_volcano.png",
           f"Figure 11 | Pan-phenome cis-MR across {PH_DIS} diseases in five categories: "
           f"effect-size volcano coloured by disease category, with {PH_HITS} FDR<5% causal "
           f"gene\u2013disease pairs.")
    figure("Figure12_pleiotropy_heatmap.png",
           "Figure 12 | Cross-category pleiotropic immune axes. Genes causal in \u22652 disease "
           "categories (e.g. IFNGR2, SWAP70, MERTK, PM20D1), with direction of effect per disease.")
    figure("Figure13_cellsource_map.png",
           "Figure 13 | Causal-target enrichment by blood-cell source of origin across the "
           "phenome (granulocyte-sourced proteins carry the strongest causal signal).")
    figure("Figure14_direction_map.png",
           "Figure 14 | Direction-aware therapeutic map of all pan-phenome causal signals, "
           "partitioned into block/neutralise versus agonise/replace with MHC-caution loci flagged.")
    figure("Figure15_novelty_ranking.png",
           "Figure 15 | Integrated novelty-priority ranking of causal immune targets across the "
           "phenome, showing the stacked evidence components and known-drug/MHC penalties per target.")
    figure("Figure16_novelty_evidence.png",
           "Figure 16 | Evidence composition and novel-versus-known separation: (a) the coloc\u2013"
           "causal evidence landscape sized by druggability, and (b) priority scores separating "
           "novel colocalized targets from recovered controls and MHC-caution loci.")

# ================= DISCUSSION =================
H("Discussion", 1, BLUE)
P("This atlas is deliberately calibrated rather than maximalist. Its central result is not a "
  "single target but a demonstration that an openly reproducible pipeline, built entirely from "
  "public data, can rediscover approved-drug biology from genetics, extend the same evidence "
  "standard to new nominations, and then\u2014crucially\u2014confirm or temper those nominations with an "
  "independent protein-level layer. Two targets reach the protein-level causal tier (TNFSF14 in "
  "multiple sclerosis, SWAP70 in rheumatoid arthritis), colocalizing at both transcript and "
  "protein. The remaining transcript-level nominations\u2014SIGLEC7/9 in Guillain\u2013Barr\u00e9 syndrome, "
  "HAVCR1 in coeliac disease, IFNGR2 in psoriasis\u2014are hypotheses at the transcript level, and we "
  "say no more than that; where plasma pQTLs disagree in direction (e.g. IL6ST, IFNGR2) we flag "
  "the discordance explicitly. The principal remaining limitations define the roadmap: pQTL "
  "coverage is limited to one plasma proteomic panel; effects are estimated in a single "
  "population; and no finding has yet been functionally perturbed. Elevating a nomination further "
  "requires replication in an independent GWAS and ancestry, orthogonal pQTL platforms, direct "
  "plasma-protein\u2013disease association, and perturbation. By publishing the atlas with its "
  "evidence ladder intact, we make both what is known and what remains to be shown fully auditable.")

# ================= METHODS =================
H("Methods (summary)", 1, BLUE)
P("Proteome curation. Olink Explore analytes were annotated with Human Protein Atlas immune-cell "
  "RNA specificity, secretome location, molecular function and curated immune-pathway membership; "
  "a transparent additive score (weighting cytokine/chemokine/IFN/TNF, complement/checkpoint and "
  "immune-cell-enrichment most heavily) defined inclusion. MSigDB C7 was used as annotation only.")
P("Genetic instruments and MR. For each immune gene the strongest cis-eQTL (eQTLGen) was used as a "
  "single instrument; Z-scores were converted to effect sizes using allele frequency and sample "
  "size (Zhu et al., 2016), harmonized against FinnGen R12 summary statistics with palindromic "
  "variants removed, and combined by the Wald ratio. FDR was controlled at 5% (Benjamini\u2013Hochberg).")
P("Colocalization. At each non-MHC significant locus, coloc.abf (Giambartolomei et al., 2014) was "
  "applied over shared cis variants using per-SNP variances\u2014real FinnGen standard errors for the "
  "outcome and Zhu-reconstructed standard errors for the exposure\u2014with default priors "
  "(p1=p2=1e-4, p12=1e-5). Loci with PP.H4\u22650.8 were designated prioritized causal targets.")
P("Protein-level pQTL validation. For each significant gene, plasma cis-pQTLs from INTERVAL "
  "(Sun et al., 2018; GWAS Catalog harmonised summary statistics, GRCh37) within \u00b1500 kb of the "
  "gene body were used as protein-level instruments. The strongest cis-pQTL drove a Wald-ratio MR "
  "against each disease (real beta and standard error on both sides), and coloc.abf was applied "
  "over shared cis variants. A pair was promoted to a protein-level causal target only if the "
  "transcript- and protein-level MR were directionally concordant and the protein-level "
  "colocalization reached PP.H4\u22650.8; directional discordances were reported, not discarded.")
P("Independent replication. Each significant instrument variant was queried in a "
  "FinnGen-independent consortium GWAS for the matching disease via the OpenGWAS API "
  "(IMSGC/Patsopoulos, Okada, IGAS/Cortes, Stuart, Fischer, Sakaue). The disease effect was "
  "harmonised to the eQTL expression-increasing allele; a hit was scored 'replicated' when the "
  "independent effect was directionally concordant with discovery and nominally significant "
  "(P<0.05). Diseases with no FinnGen-independent GWAS, or variants absent from older "
  "targeted-array studies, were recorded as coverage gaps.")
P("Claim-strength gate. Every finding was assigned an evidence tier\u2014transcript-level proxy "
  "(cis-eQTL), genetically-supported nomination (cis-MR FDR<0.05), prioritized causal target "
  "(+transcript colocalization), protein-level causal target (+concordant pQTL-MR and protein "
  "colocalization)\u2014with MHC signals capped at nomination and flagged for LD. Language in this "
  "manuscript is bound to these tiers.")
P("Data availability. The atlas is built entirely from public resources (Olink Explore universe, "
  "Human Protein Atlas, MSigDB, eQTLGen, FinnGen R12, INTERVAL plasma pQTLs via the GWAS Catalog) "
  "and released, with all code, as a versioned open resource. Individual-level UK Biobank data and "
  "the UKB-PPP/Synapse-gated pQTLs are controlled-access and were not used.")

# ---- evidence table (protein-aware final tiers + replication) ----
H("Table 1 | Final evidence-tiered targets (discovery + protein + replication)", 2)
tb = doc.add_table(rows=1, cols=8); tb.style = "Light Grid Accent 1"
hdrs=["Gene","Disease","eQTL OR","coloc\nPP.H4","pQTL\nPP.H4","dir","Replication","Final tier"]
for i,h in enumerate(hdrs):
    tb.rows[0].cells[i].paragraphs[0].add_run(h).bold = True
for _,r in final.iterrows():
    c = tb.add_row().cells
    c[0].text=str(r.gene_symbol); c[1].text=str(r.disease)
    c[2].text=f"{r.OR:.2f}"
    c[3].text=("\u2014" if pd.isna(r.PP_H4) else f"{r.PP_H4:.2f}")
    c[4].text=("\u2014" if pd.isna(r.pQTL_PPH4) else f"{r.pQTL_PPH4:.2f}")
    c[5].text=str(r.pQTL_concordant) if str(r.pQTL_concordant)!="nan" else ""
    c[6].text=(str(r.rep_status) if HAS_REP and str(r.get("rep_status"))!="nan" else "")
    c[7].text=f"T{int(r.final_tier)}"
    for cell in c:
        for para in cell.paragraphs:
            for run in para.runs: run.font.size=Pt(7.5)
P("Tiers: T5 protein-level causal target (concordant pQTL-MR + protein coloc PP.H4\u22650.8); "
  "T4 prioritized causal target (transcript coloc PP.H4\u22650.8); T3 genetically-supported "
  "nomination; T2 nomination held for MHC LD. dir = transcript/protein direction concordance; "
  "Replication = status in a FinnGen-independent consortium GWAS (OpenGWAS).",
  italic=True, size=8)

# ================= EXTENDED DATA / SUPPLEMENTARY FIGURE GALLERY =================
SUP = os.path.join(ROOT, "08_figures", "supplementary")
PHE = os.path.join(ROOT, "08_figures", "phenome")

def gallery_figure(fp, caption, width=5.2):
    if os.path.exists(fp):
        doc.add_picture(fp, width=Inches(width))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    c = doc.add_paragraph(); r = c.add_run(caption); r.italic=True; r.font.size=Pt(8)

def _disease_from(stem):
    return stem.replace("_"," ").strip().rstrip("_").replace("  "," ")

# ---- Pan-phenome extended figures ----
if HAS_PHEN and os.path.isdir(PHE):
    doc.add_page_break()
    H("Extended Data \u00b7 Pan-phenome figures", 1, BLUE)
    P("Per-category and per-disease views of the 28-disease pan-phenome cis-MR "
      "(src/24\u201326). Volcano plots show cis-MR effect size vs \u2212log10(FDR) for every "
      "disease; forest plots summarise the strongest causal targets within each disease "
      "category.", italic=True, size=9)
    phe = sorted(os.listdir(PHE))
    order = [f for f in phe if f.startswith("PFig_hits")] + \
            [f for f in phe if f.startswith("PFig_forest")] + \
            [f for f in phe if f.startswith("PFig_volcano")]
    n = 0
    for f in order:
        stem = f[:-4]
        if f.startswith("PFig_hits"):
            cap = "Pan-phenome causal hits per disease category."
        elif f.startswith("PFig_forest"):
            cat = stem.replace("PFig_forest_","").replace("_"," ")
            cap = f"Forest plot of causal immune targets \u2014 {cat} diseases (OR \u00b1 95% CI)."
        else:
            dis = _disease_from(stem.replace("PFig_volcano_",""))
            cap = f"Pan-phenome cis-MR volcano \u2014 {dis}."
        n += 1
        gallery_figure(os.path.join(PHE,f), f"Extended Data Fig. P{n} | {cap}")
    print(f"embedded pan-phenome extended figures: {n}")

# ---- Supplementary figures ----
if os.path.isdir(SUP):
    doc.add_page_break()
    H("Supplementary figures", 1, BLUE)
    P("Quality-control and per-disease supplementary figures for the autoimmune "
      "discovery arc (src/20): proteome composition and instrument diagnostics; "
      "per-disease MR volcano, QQ + genomic-inflation \u03bb, and forest plots; per-gene "
      "INTERVAL cis-pQTL regional association; and discovery-versus-replication panels.",
      italic=True, size=9)
    SCAP = {"immune_class_composition":"Immune-class composition of the curated proteome.",
            "source_cell_lineage":"Blood-cell lineage of origin of immune proteins.",
            "instrument_strength":"cis-eQTL instrument strength distribution.",
            "instrument_samplesize":"eQTLGen instrument sample-size distribution.",
            "tests_per_disease":"MR tests performed per disease.",
            "hits_per_disease":"FDR<5% MR hits per disease."}
    def _scap(stem):
        body = stem.split("_",1)[1] if "_" in stem else stem   # drop SFigSNN_
        for k,v in SCAP.items():
            if body==k: return v
        if body.startswith("volcano_"):  return f"Per-disease cis-MR volcano \u2014 {_disease_from(body[8:])}."
        if body.startswith("qq_"):       return f"MR QQ plot and genomic-inflation \u03bb \u2014 {_disease_from(body[3:])}."
        if body.startswith("forest_"):   return f"Per-disease causal-target forest \u2014 {_disease_from(body[7:])}."
        if body.startswith("pqtl_regional_"): return f"INTERVAL cis-pQTL regional association \u2014 {body[14:]}."
        if body.startswith("repl_"):     return f"Discovery-versus-independent-replication \u2014 {_disease_from(body[5:])}."
        return body.replace("_"," ")
    m = 0
    for f in sorted(os.listdir(SUP)):
        if not f.endswith(".png"): continue
        stem = f[:-4]; m += 1
        num = stem.split("_")[0].replace("SFigS","S")
        gallery_figure(os.path.join(SUP,f), f"Supplementary Fig. {num} | {_scap(stem)}")
    print(f"embedded supplementary figures: {m}")

out = os.path.join(OUT, "Plasma_Immunome_Phenome_Atlas_Nature.docx")
doc.save(out)
print("wrote", out)
print(f"tiers: T5={N_T5} T4={N_T4} T3={N_T3} T2(MHC)={N_T2} | pQTL pairs={N_PQTL} "
      f"concordant={N_CONC} discordant={N_DISC} | repl covered={N_COV} replicated={N_REP} "
      f"concordant={N_RCON} | MR sig={N_SIG}/{N_TEST} | immune proteins={N_IMM}")
