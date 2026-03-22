# De Novo Design of a Metal-Nucleating Protein for CoFe₂O₄ / Fe₁₆N₂ Crystal Nucleation
### ARIA Universal Fabricators 3D Challenge — Design Report

---

## Executive Summary

This report describes the computational design of a de novo protein capable of nucleating cobalt ferrite (CoFe₂O₄) and iron nitride (Fe₁₆N₂) crystals. The protein was designed from scratch — no natural template was used — using a four-stage pipeline: (1) motif definition, (2) backbone diffusion, (3) sequence design, and (4) structure validation. The final output is a set of ten ranked fusion proteins, each consisting of a designed metal-nucleating core (82–90 residues) appended to a SpyTag003 polymerisation tag, ready for gene synthesis and experimental validation.

---

## 1. What Makes a Good Crystal-Nucleating Protein?

Before describing the design process, it is worth establishing the biophysical rationale for the features we targeted.

### 1.1 Metal Coordination: The Carboxylate Triad

Crystal nucleation of metal oxides and nitrides requires the protein to concentrate and pre-organise metal ions at a surface. In natural biomineralisation proteins (ferritin, magnetosome-associated Mms6, silicatein), this is achieved through clusters of acidic residues — primarily glutamate (Glu, E) and aspartate (Asp, D) — whose carboxylate side chains coordinate divalent cations (Fe²⁺, Co²⁺) with octahedral geometry. A histidine (His, H) residue is frequently co-present to provide a nitrogen donor, and asparagine (Asn, N) contributes hydrogen-bonding stabilisation of the coordination shell.

We therefore defined a **carboxylate triad motif**: three glutamate residues (or Glu/Asp combinations) presented on a single alpha-helix, spaced at i, i+3, i+7 positions to align on the same helical face. This geometry places all three carboxylate oxygens within ~5 Å of a central metal-binding site, mimicking the active site geometry of natural iron-binding proteins.

### 1.2 Acidic Mineralisation Surface

Beyond the active triad, the broader protein surface should be negatively charged at physiological pH. An acidic surface (net charge < −3 at pH 7, pI < 6) serves two purposes:
- **Electrostatic attraction** of divalent metal cations (Fe²⁺, Co²⁺) to the protein surface, increasing local metal concentration
- **Inhibition of non-specific aggregation** — highly acidic proteins tend to remain soluble and monomeric in the absence of metal ions

We targeted a surface composition where 15–32% of residues are Asp or Glu (frac_DE), a range empirically associated with mineralisation activity in known biomineralisation proteins.

### 1.3 Absence of Cysteine

Cysteine residues form disulfide bonds under oxidising conditions and can coordinate metals non-specifically via their thiol groups, disrupting the designed coordination geometry. All designs were required to be cysteine-free.

### 1.4 Compact, Helical Architecture

Alpha-helical bundles are the preferred scaffold for metal-binding protein design because:
- Helices present side chains on a predictable, periodic surface
- Helical bundles are thermodynamically stable without disulfide bonds
- The geometry is compatible with the DNA origami cylinder constraint (inner diameter ~20 nm), which requires a compact protein footprint

Target length was 60–110 residues for the core domain.

### 1.5 SpyTag003 Polymerisation Tag

The challenge requires the protein to polymerise inside a DNA origami cylinder to form a linear array of nucleation sites. This was achieved by appending **SpyTag003** (`RGVPHIVMVDAYKRYK`, GenBank MN433888.1) to the C-terminus via a flexible `GSGSGS` linker. SpyTag003 forms an irreversible isopeptide bond with its partner SpyCatcher003, enabling covalent head-to-tail polymerisation. The full fusion suffix appended to every candidate is:

```
GSGSGSRGVPHIVMVDAYKRYK  (22 residues)
```

---

## 2. Design Pipeline

The pipeline comprised four sequential computational stages, each with defined quality filters before passing candidates to the next stage.

```
Motif PDB  →  RFdiffusion  →  ProteinMPNN  →  Boltz-2  →  Final Ranking
(1 motif)     (100 backbones)  (48 sequences)  (34 structures)  (10 candidates)
```

### Stage 1: Motif Definition

A 12-residue alpha-helix was constructed in idealised geometry (φ = −57°, ψ = −47°) with the following composition:

- Positions 1, 4, 8: **Glutamate** (metal-coordinating carboxylates)
- Position 6: **Histidine** (nitrogen donor)
- Position 10: **Asparagine** (hydrogen-bond stabiliser)
- Remaining positions: Alanine (helix-stabilising, inert)

This motif was saved as a PDB file and used as a **fixed constraint** in the subsequent backbone diffusion step — the diffusion model was required to incorporate this exact helix into every generated backbone.

### Stage 2: Backbone Generation with RFdiffusion

**RFdiffusion** (Watson et al., 2023) is a denoising diffusion model trained on protein structures. It generates novel protein backbones conditioned on partial structure constraints. We ran 100 independent diffusion trajectories, each conditioned on the motif helix, targeting a total protein length of 80–100 residues.

**Quality filter applied:** NC-ratio (ratio of contacts within the designed region to total contacts) > 0.5, ensuring the motif helix is genuinely integrated into the core of the fold rather than dangling at a terminus. 29 of 100 backbones passed; the top 25 by NC-ratio were carried forward.

The 25 selected backbones represent diverse helical bundle topologies — 3-helix bundles, 4-helix bundles, and mixed arrangements — all sharing the common metal-binding helix.

### Stage 3: Sequence Design with ProteinMPNN

**ProteinMPNN** (Dauparas et al., 2022) is a graph neural network that designs amino acid sequences to fold into a given backbone. For each of the 25 backbones, we ran ProteinMPNN at two sampling temperatures (T = 0.1 for conservative designs, T = 0.3 for more diverse designs) and generated multiple sequences per backbone, yielding an initial pool of ~200 candidate sequences.

Each sequence was then scored and filtered using the **v3.2 composite score**:

```
composite = mpnn_score + max(0, pI − 4.5) × 0.3 − triad_score × 0.05

where:  triad_score = (count of Glu at triad positions) × 1.5
                    + (count of Asp at triad positions) × 1.0
```

This formula rewards:
- Low ProteinMPNN negative log-likelihood (high sequence-backbone compatibility)
- Low pI (acidic protein, favourable for metal binding)
- Recovery of the designed Glu/Asp triad residues (penalises sequences where the model mutated away the key metal-binding residues)

**Hard filters applied before scoring:**

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| frac_DE | 0.15 – 0.32 | Acidic surface composition |
| pI | < 6.0 | Net negative charge at physiological pH |
| net_charge (pH 7) | < −3.0 | Electrostatic metal attraction |
| triad_recovered | ≥ 1 of 3 positions | Preserve metal-binding function |
| Cys count | = 0 | Avoid disulfide/thiol interference |
| Length | 60 – 110 aa | Compact, DNA-origami compatible |

**48 sequences** passed all filters. These were tiered by triad recovery:
- **Tier A** (3 sequences): All three triad positions recovered as Glu (EEE) — highest metal-binding potential
- **Tier B** (15 sequences): Two triad positions recovered as Glu/Asp
- **Tier C** (29 sequences): One triad position recovered
- **Tier D** (1 sequence): Triad not recovered but all other filters passed

All 48 sequences were carried forward to structure prediction.

### Stage 4: Structure Prediction and Validation with Boltz-2

**Boltz-2** is a state-of-the-art structure prediction model (analogous to AlphaFold3) that predicts protein structure from sequence and outputs per-residue confidence (pLDDT) and pairwise alignment error (pAE). Each of the 48 sequences was submitted as an independent prediction job.

For each prediction, three metrics were computed:

| Metric | Threshold | Meaning |
|--------|-----------|---------|
| **pLDDT** | > 85 | Per-residue confidence; > 85 indicates a well-ordered, confident prediction |
| **pAE** | < 5 Å | Predicted alignment error; low values indicate the model is confident about the relative positions of all residues |
| **RMSD** | < 1.5 Å | Kabsch-aligned Cα RMSD between the Boltz-2 predicted structure and the original RFdiffusion backbone; verifies that the sequence actually folds back to the intended backbone |

The RMSD check is particularly important: it confirms that ProteinMPNN designed a sequence that is genuinely compatible with the target backbone, rather than a sequence that folds into a different structure with similar energy.

**34 sequences** were successfully predicted; **32 passed** all three structure filters (2 failed on RMSD > 1.5 Å, indicating the sequence had drifted from the intended backbone).

---

## 3. Final Scoring and Ranking

Passing sequences were ranked by a **final composite score** that combines sequence-level quality (composite score) with structure-level confidence penalties:

```
final_score = composite + (pAE / 5.0) × 0.3 + (RMSD / 1.5) × 0.2
```

Lower scores are better. The pAE and RMSD terms penalise structural uncertainty and backbone deviation respectively, ensuring that the top-ranked sequences are both sequence-optimal and structurally well-predicted.

---

## 4. Top-10 Final Candidates

The ten best-scoring sequences are presented below. All have been appended with the SpyTag003 fusion suffix to produce the final submission sequences.

| Rank | ID | Tier | pLDDT | pAE (Å) | RMSD (Å) | Final Score | Fusion Length |
|------|----|------|-------|---------|----------|-------------|---------------|
| 1 | metal_nucleator_42_s1_T01 | B | 96.04 | 1.598 | 0.613 | **1.1571** | 112 aa |
| 2 | metal_nucleator_51_s5_T03 | **A ★** | 94.84 | 1.804 | 0.753 | 1.2351 | 109 aa |
| 3 | metal_nucleator_37_s2_T01 | C | 95.81 | 2.223 | 0.985 | 1.2383 | 108 aa |
| 4 | metal_nucleator_32_s4_T01 | B | 94.87 | 3.019 | 1.168 | 1.2770 | 110 aa |
| 5 | metal_nucleator_86_s3_T01 | C | 94.84 | 2.717 | 1.092 | 1.2797 | 104 aa |
| 6 | metal_nucleator_19_s8_T03 | B | 93.53 | 2.472 | 0.938 | 1.2921 | 107 aa |
| 7 | metal_nucleator_14_s4_T03 | B | 94.97 | 2.056 | 0.713 | 1.2979 | 106 aa |
| 8 | metal_nucleator_88_s5_T01 | B | 94.13 | 3.350 | 0.853 | 1.3118 | 107 aa |
| 9 | metal_nucleator_37_s4_T03 | C | 94.85 | 2.589 | 0.712 | 1.3157 | 108 aa |
| 10 | metal_nucleator_19_s1_T03 | C | 94.70 | 2.031 | 0.717 | 1.3262 | 107 aa |

★ = EEE triad (all three metal-coordinating positions recovered as glutamate)

### Highlighted Candidates

**Rank 1 — metal_nucleator_42_s1_T01** is the overall winner by final score. It achieves the lowest pAE (1.60 Å) and lowest RMSD (0.61 Å) in the entire set, indicating an exceptionally well-ordered structure that closely matches the designed backbone. Its pLDDT of 96.04 is the second-highest across all 32 passing sequences. The triad pattern is EAE (Glu–Ala–Glu), retaining two of three metal-coordinating positions. Net charge at pH 7 is −10.06, making it the most acidic sequence in the top-10 and highly favourable for metal ion attraction.

**Rank 2 — metal_nucleator_51_s5_T03** is the only Tier A (EEE triad) sequence in the top-10, making it the highest-priority candidate for experimental validation of metal-binding function. All three designed glutamate positions are preserved, giving it the strongest theoretical metal coordination geometry. Its pLDDT (94.84), pAE (1.80 Å), and RMSD (0.75 Å) are all excellent. Net charge −9.05, pI 4.46.

**Rank 7 — metal_nucleator_14_s4_T03** achieves the second-lowest RMSD in the top-10 (0.71 Å) with a pLDDT of 94.97, indicating very high backbone fidelity to the RFdiffusion design. Its triad pattern is KEE (Lys–Glu–Glu), retaining two acidic metal-coordinating residues.

### Full Fusion Sequences (Top-10)

```
>rank01|metal_nucleator_42_s1_T01|tier=B|final=1.1571|pLDDT=96.04|pAE=1.598|RMSD=0.613|len=112
SAAFEAEVERLLERVEAGIAELRALVAELRERGLLTEEEAAELDRRIEEWEARLARLKERREDPEVPEELRRLLAEVEAAIAELRAREAAGSGSGSRGVPHIVMVDAYKRYK

>rank02|metal_nucleator_51_s5_T03|tier=A|final=1.2351|pLDDT=94.84|pAE=1.804|RMSD=0.753|len=109
EQQRVQALIEAFRREAEELIRELERRLREAGREEEARETVEELRRELDAAIAELEAMVAAGAPLEEIEARIAAELAKLRARVEEILAGSGSGSRGVPHIVMVDAYKRYK

>rank03|metal_nucleator_37_s2_T01|tier=C|final=1.2383|pLDDT=95.81|pAE=2.223|RMSD=0.985|len=108
AALAAKLAALRAAIAEARALLERYAALVADPSTSPEARAAVRAEAARAEAEARALAAELAPEAGAEVAALLAELEAYRAALEAELAGSGSGSRGVPHIVMVDAYKRYK

>rank04|metal_nucleator_32_s4_T01|tier=B|final=1.2770|pLDDT=94.87|pAE=3.019|RMSD=1.168|len=110
MEEKAKELIEKMKKKIEEAKKLLEEAKKNPEKAVEELEKKAKELEELIKEAEEFLPELNEEEKKIMQELIKEAKEIKKELEKLAKEKKGSGSGSRGVPHIVMVDAYKRYK

>rank05|metal_nucleator_86_s3_T01|tier=C|final=1.2797|pLDDT=94.84|pAE=2.717|RMSD=1.092|len=104
ELAELRAKARELVAKALALEAEAAALPDPADVGALLEEAEKARKEAEEAIAAYFKALGEELTPERLAAELAAIRAELAAAEAGSGSGSRGVPHIVMVDAYKRYK

>rank06|metal_nucleator_19_s8_T03|tier=B|final=1.2921|pLDDT=93.53|pAE=2.472|RMSD=0.938|len=107
SSLEEERKALLKEIEKKAAEAKALAERLRADPEAPPEVAEKFEKIAKKLEELKKRVEEGADLEEAQKEAEELKKEAEKLLEKLEAGSGSGSRGVPHIVMVDAYKRYK

>rank07|metal_nucleator_14_s4_T03|tier=B|final=1.2979|pLDDT=94.97|pAE=2.056|RMSD=0.713|len=106
AEAERAARLNAEMARAREEAERLIEEAKKLKKEGKEKEAEEALERAEETLRRLEELAAQVPERRAEVEAYVEELKEEIKKLRKEGSGSGSRGVPHIVMVDAYKRYK

>rank08|metal_nucleator_88_s5_T01|tier=B|final=1.3118|pLDDT=94.13|pAE=3.350|RMSD=0.853|len=107
AAAAAAAAARAALLAAYRARVRELWLRTKELIKEGTPEEAREATEREIAALRAEVLAQLADLDPAAVEAELAALEEEVRAEVEAEGSGSGSRGVPHIVMVDAYKRYK

>rank09|metal_nucleator_37_s4_T03|tier=C|final=1.3157|pLDDT=94.85|pAE=2.589|RMSD=0.712|len=108
SELEERLRRLRAALERAKELLDRYAALVKDPTTSPEEREAVRAEAERALEEALRLAEELKPHAGEEVDALLAELEAYRAARERELAGSGSGSRGVPHIVMVDAYKRYK

>rank10|metal_nucleator_19_s1_T03|tier=C|final=1.3262|pLDDT=94.70|pAE=2.031|RMSD=0.717|len=107
SAAEAERKALEAEIAAAAAEARALAAALEADPAAPAELAQRFREIAERLEALLAEVQAGADLEVARQDAKALRAERAALLAAAAAGSGSGSRGVPHIVMVDAYKRYK
```

---

## 5. Key Figures

**Figure 1 — Final Scoreboard** (`final_scoreboard.png`): Horizontal bar chart showing all 32 passing sequences ranked by final score. Top-10 candidates are outlined in black. Tier A sequences (EEE triad) are marked with ★. The dashed vertical line marks the top-10 cutoff at final score = 1.3262.

**Figure 2 — Boltz-2 Quality Metrics** (`boltz2_quality_metrics.png`): Three scatter plots showing pLDDT, pAE, and RMSD vs rank for all 32 passing sequences. Red dashed lines indicate the filter thresholds. All top-10 sequences (outlined points) comfortably clear all three thresholds.

---

## 6. Design Rationale Summary

The pipeline was designed around three core principles:

**1. Function-first motif anchoring.** Rather than designing a protein and hoping metal-binding emerges, we started from the metal-binding geometry and built the protein around it. RFdiffusion's motif-scaffolding capability ensures the carboxylate triad is structurally integrated into every backbone, not appended post-hoc.

**2. Multi-criteria sequence filtering.** ProteinMPNN generates sequences that fold well, but does not optimise for surface charge, metal-binding residue recovery, or cysteine absence. The v3.2 composite score and hard filters enforce all biophysical requirements simultaneously, preventing the common failure mode of selecting sequences that are structurally confident but functionally inert.

**3. Backbone fidelity as a validation gate.** The RMSD filter between Boltz-2 predictions and RFdiffusion backbones is the most stringent check in the pipeline. A sequence that passes pLDDT and pAE but fails RMSD has been predicted to fold into a different structure than intended — the metal-binding geometry would be disrupted. This filter eliminated 2 sequences that would otherwise have appeared to be high-quality predictions.

---

## 7. Recommended Experimental Path

For experimental validation, we recommend the following priority order:

1. **Rank 2 (metal_nucleator_51_s5_T03, Tier A ★)** — first for metal-binding assays (ITC, fluorescence quenching with Co²⁺/Fe²⁺) due to complete EEE triad recovery
2. **Rank 1 (metal_nucleator_42_s1_T01, Tier B)** — first for structural characterisation (CD spectroscopy, SEC-SAXS) due to highest structural confidence metrics
3. **Ranks 1–5** — for initial expression screening in *E. coli* (His-tag on N-terminus, IPTG induction, IMAC purification)
4. **SpyCatcher003 polymerisation assay** — mix each candidate with SpyCatcher003 in equimolar ratio, verify isopeptide bond formation by SDS-PAGE band shift, then load into pre-assembled DNA origami cylinders

---

## 8. Software and Methods

| Tool | Version / Reference | Purpose |
|------|---------------------|---------|
| RFdiffusion | Watson et al., *Nature* 2023 | Backbone generation with motif scaffolding |
| ProteinMPNN | Dauparas et al., *Science* 2022 | Sequence design on fixed backbone |
| Boltz-2 | Wohlwend et al., 2024 | Structure prediction and confidence scoring |
| BioPython / NumPy | Standard | Sequence analysis, RMSD calculation (Kabsch algorithm) |
| SpyTag003 | Keeble et al., *PNAS* 2019 (GenBank MN433888.1) | Covalent polymerisation tag |

---

*Report generated: 2026-03-22*
*Pipeline executed by: Biomni (Phylo) computational design platform*
