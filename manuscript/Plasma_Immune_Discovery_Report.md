# Disease-Trained PIRS — Plasma Immune Discovery Report

## Intelligence layer over the genetics-anchored plasma immunome atlas

### Mode

**Causal-atlas-only mode.** No trained PIRS was found in `05_machine_learning/`, so PIRS-dependent columns (coefficient, CV stability, sensitivity/specificity/AUROC contribution) are reported as `NA (train PIRS)`. Every causal, novelty, direction and tier column is populated from real genetic data. Train a PIRS (`src/22_train_pirs.py`) on an authorised cohort to activate the predictive columns — no values are fabricated.

### Result summary

- Gene–disease pairs assessed: **176** across **25 diseases** (5 categories).

- High-novelty targets (Tier 4/5): **45** (Tier 5 protein-level novel: **0**; Tier 4 prioritized: **45**).

- Known-drug positive controls (Tier 1): **12** — recovered from genetics alone, validating the pipeline.

- Independently replicated pairs: **17**.

- Protein-level causal pairs (transcript coloc + concordant plasma pQTL): **6**.

### Why Tier 5 is currently empty (an honest gate, not a gap in the run)

Tier 5 requires **all five** of: prediction-ready + colocalization + concordant plasma pQTL + specificity + druggability. Two hits reach protein-level causality (transcript coloc + concordant INTERVAL pQTL): **TNFSF14 → multiple sclerosis** and **SWAP70 → rheumatoid arthritis**. TNFSF14 is flagged as a **known-drug axis** → it is reported as a **Tier 1 positive control** (correct pipeline recovery, not a novel target). SWAP70 is **novel and protein-level causal** but has **druggability = 0** (intracellular, not secreted/membrane) → it lands at **Tier 4**, not Tier 5. Separately, the **pan-phenome diseases (cardiovascular / metabolic / renal / neuro) have no plasma pQTL layer computed yet** (INTERVAL pQTL was run only against the autoimmune arc), so none of those hits can reach Tier 5 until a plasma pQTL panel is colocalized against them. This is a data-coverage boundary, not a fabricated ceiling — running cis-pQTL MR+coloc across the phenome is the single step that can promote Tier-4 targets to Tier 5.

### Claim discipline

Claims are bound to the evidence level actually reached:

- PIRS weight alone → **biomarker** (predictive, not causal).

- + cis-MR → **causal nomination**.

- + colocalization → **prioritized causal target**.

- + direction-concordant plasma pQTL → **protein-level causal target**.

- + perturbation → **proof of mechanism** (not asserted here; proposed as next experiment).

All signals restricted to **plasma immune proteins**; no single-cell, tissue or intracellular claim is made unless such data are separately supplied.

### Tier 4 — prioritized causal targets (obtain plasma pQTL to upgrade)

- **ACE → Dementia**: OR=0.44, PP.H4=0.995, Agonize / recombinant-protein replacement (protective, secreted) | High (FDA-precedented class, tractable).

- **IFNGR2 → Hypertension**: OR=1.066, PP.H4=0.985, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **ERBB3 → Type 1 diabetes**: OR=0.424, PP.H4=0.901, Replace / decoy-receptor or agonist (protective, soluble) | Moderate.

- **IL2RA → Type 1 diabetes**: OR=0.302, PP.H4=0.994, Replace / decoy-receptor or agonist (protective, soluble) | High (FDA-precedented class, tractable).

- **PLAUR → Coronary heart disease**: OR=1.587, PP.H4=0.899, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **IFNGR2 → Psoriasis**: OR=0.876, PP.H4=0.974, Replace / decoy-receptor or agonist (protective, soluble) | High (FDA-precedented class, tractable).

- **PSRC1 → Coronary heart disease**: OR=0.85, PP.H4=0.988, Agonize (protective; intracellular route caution) | Uncertain.

- **FES → Hypertension**: OR=0.847, PP.H4=0.978, Agonize (protective; intracellular route caution) | Uncertain.

- **ACE → Atrial fibrillation**: OR=1.724, PP.H4=0.878, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **SWAP70 → Rheumatoid arthritis**: OR=0.829, PP.H4=0.966, Agonize (protective; intracellular route caution) | Uncertain.

- **ACE → Alzheimer's disease**: OR=0.43, PP.H4=0.978, Agonize / recombinant-protein replacement (protective, secreted) | High (FDA-precedented class, tractable).

- **SPINK8 → Type 2 diabetes**: OR=0.834, PP.H4=0.974, Agonize / recombinant-protein replacement (protective, secreted) | Low-moderate.

- **MERTK → Hypertension**: OR=1.062, PP.H4=0.872, Block (antagonist / neutralizing antibody / small molecule) | Low-moderate.

- **PM20D1 → Type 1 diabetes**: OR=0.758, PP.H4=0.918, Agonize / recombinant-protein replacement (protective, secreted) | Low-moderate.

- **ANXA2 → Hypertension**: OR=1.11, PP.H4=0.997, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

### Figure panels

- Panel A — plasma immune discovery workflow.

- Panel B — disease-specific PIRS performance (causal signal strength if untrained).

- Panel C — plasma immune signature (high-novelty targets, direction-coloured).

- Panel D — causal-evidence concordance ladder (MR→coloc→pQTL→replication).

- Panel E — novelty-priority map (PP.H4 vs PINS, tier-coloured).

- Panel F — validation plan for the top novel targets.

### Final Required Output Table

See `intelligence_layer_final_table.tsv` (ranked; high-novelty Tier 4/5 first). Columns: Rank, Disease, Protein, Gene, Protein class, Plasma detectability, PIRS coefficient, PIRS direction, CV stability, Sensitivity/Specificity/AUROC contribution, MR OR/FDR/support, Coloc PP.H4/support, pQTL support, Replication, Known-drug status, Druggability, Novelty score (PINS), Novelty tier (+label), Causal-predictive concordance, Therapeutic direction, Biomarker-or-target, Best figure panel, Best validation experiment, Final recommendation.
