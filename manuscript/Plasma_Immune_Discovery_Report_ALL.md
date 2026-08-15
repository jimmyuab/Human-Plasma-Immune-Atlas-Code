# Disease-Trained PIRS — Plasma Immune Discovery Report

## Intelligence layer over the genetics-anchored plasma immunome atlas

### Mode

**Causal-atlas-only mode.** No trained PIRS was found in `05_machine_learning/`, so PIRS-dependent columns (coefficient, CV stability, sensitivity/specificity/AUROC contribution) are reported as `NA (train PIRS)`. Every causal, novelty, direction and tier column is populated from real genetic data. Train a PIRS (`src/22_train_pirs.py`) on an authorised cohort to activate the predictive columns — no values are fabricated.

### Result summary

- Gene–disease pairs assessed: **1016** across **380 diseases** (5 categories).

- High-novelty targets (Tier 4/5): **364** (Tier 5 protein-level novel: **5**; Tier 4 prioritized: **359**).

- Known-drug positive controls (Tier 1): **90** — recovered from genetics alone, validating the pipeline.

- Independently replicated pairs: **4**.

- Protein-level causal pairs (transcript coloc + concordant plasma pQTL): **39**.

### Claim discipline

Claims are bound to the evidence level actually reached:

- PIRS weight alone → **biomarker** (predictive, not causal).

- + cis-MR → **causal nomination**.

- + colocalization → **prioritized causal target**.

- + direction-concordant plasma pQTL → **protein-level causal target**.

- + perturbation → **proof of mechanism** (not asserted here; proposed as next experiment).

All signals restricted to **plasma immune proteins**; no single-cell, tissue or intracellular claim is made unless such data are separately supplied.

### Tier 5 — novel protein-level plasma-immune targets

**IL2RA → Type 1 diabetes, definitions combined** (PINS 6.154). MR OR=0.302 (FDR 1.58e-06), coloc PP.H4=0.994, plasma pQTL concordant, replication: not tested. Direction: Replace / decoy-receptor or agonist (protective, soluble). Next: cis-pQTL colocalization already met; confirm with CRISPRi/CRISPRa over-expression/supplementation in primary immune cells + plasma NPX dose-response.

**PPP3R1 → Varicose veins** (PINS 5.596). MR OR=1.415 (FDR 3.27e-36), coloc PP.H4=0.929, plasma pQTL concordant, replication: replicated (FinnGen + UK Biobank). Direction: Block (antagonist / neutralizing antibody / small molecule). Next: cis-pQTL colocalization already met; confirm with CRISPRi/CRISPRa neutralization/knockdown in primary immune cells + plasma NPX dose-response.

**PPP3R1 → Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified** (PINS 5.563). MR OR=1.246 (FDR 1.68e-24), coloc PP.H4=0.896, plasma pQTL concordant, replication: replicated (FinnGen + UK Biobank). Direction: Block (antagonist / neutralizing antibody / small molecule). Next: cis-pQTL colocalization already met; confirm with CRISPRi/CRISPRa neutralization/knockdown in primary immune cells + plasma NPX dose-response.

**ANXA2 → Hypertension** (PINS 4.02). MR OR=1.11 (FDR 2.42e-03), coloc PP.H4=0.997, plasma pQTL concordant, replication: not tested. Direction: Block (antagonist / neutralizing antibody / small molecule). Next: cis-pQTL colocalization already met; confirm with CRISPRi/CRISPRa neutralization/knockdown in primary immune cells + plasma NPX dose-response.

**ANXA2 → Hypertension, essential** (PINS 3.721). MR OR=1.098 (FDR 4.02e-02), coloc PP.H4=0.942, plasma pQTL concordant, replication: not tested. Direction: Block (antagonist / neutralizing antibody / small molecule). Next: cis-pQTL colocalization already met; confirm with CRISPRi/CRISPRa neutralization/knockdown in primary immune cells + plasma NPX dose-response.

### Tier 4 — prioritized causal targets (obtain plasma pQTL to upgrade)

- **ACE → Antihypertensive medication - note that there are other indications**: OR=2.121, PP.H4=0.998, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **ERBB3 → Disorders of the thyroid gland**: OR=0.78, PP.H4=0.917, Replace / decoy-receptor or agonist (protective, soluble) | Moderate.

- **ACE → Cardiovascular diseases (excluding rheumatic etc)**: OR=1.584, PP.H4=0.997, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **ERBB2 → Asthma/COPD (KELA code 203)**: OR=2.535, PP.H4=0.868, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **ERBB3 → Type 1 diabetes, definitions combined**: OR=0.424, PP.H4=0.901, Replace / decoy-receptor or agonist (protective, soluble) | Moderate.

- **LTBP3 → Otosclerosis**: OR=0.336, PP.H4=1.0, Agonize / recombinant-protein replacement (protective, secreted) | Low-moderate.

- **ERBB3 → Asthma/COPD (KELA code 203)**: OR=0.799, PP.H4=0.983, Replace / decoy-receptor or agonist (protective, soluble) | Moderate.

- **IL2RA → Dermatitis and eczema**: OR=1.487, PP.H4=1.0, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **IL2RA → Dermatitis**: OR=1.647, PP.H4=1.0, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **ACE → Dementia**: OR=0.44, PP.H4=0.995, Agonize / recombinant-protein replacement (protective, secreted) | High (FDA-precedented class, tractable).

- **ACE → Any dementia (more control exclusions)**: OR=0.414, PP.H4=0.994, Agonize / recombinant-protein replacement (protective, secreted) | High (FDA-precedented class, tractable).

- **ACE → Hypertension, essential**: OR=1.534, PP.H4=0.917, Block (antagonist / neutralizing antibody / small molecule) | High (FDA-precedented class, tractable).

- **CFH → Age-related macular degeneration (whether dry or wet)**: OR=0.004, PP.H4=0.981, Agonize / recombinant-protein replacement (protective, secreted) | Low-moderate.

- **CFH → Other retinal disorders**: OR=0.16, PP.H4=0.978, Agonize / recombinant-protein replacement (protective, secreted) | Low-moderate.

- **CFH → Moderate visual impairment, binocular**: OR=0.018, PP.H4=0.977, Agonize / recombinant-protein replacement (protective, secreted) | Low-moderate.

### Figure panels

- Panel A — plasma immune discovery workflow.

- Panel B — disease-specific PIRS performance (causal signal strength if untrained).

- Panel C — plasma immune signature (high-novelty targets, direction-coloured).

- Panel D — causal-evidence concordance ladder (MR→coloc→pQTL→replication).

- Panel E — novelty-priority map (PP.H4 vs PINS, tier-coloured).

- Panel F — validation plan for the top novel targets.

### Final Required Output Table

See `intelligence_layer_final_table.tsv` (ranked; high-novelty Tier 4/5 first). Columns: Rank, Disease, Protein, Gene, Protein class, Plasma detectability, PIRS coefficient, PIRS direction, CV stability, Sensitivity/Specificity/AUROC contribution, MR OR/FDR/support, Coloc PP.H4/support, pQTL support, Replication, Known-drug status, Druggability, Novelty score (PINS), Novelty tier (+label), Causal-predictive concordance, Therapeutic direction, Biomarker-or-target, Best figure panel, Best validation experiment, Final recommendation.
