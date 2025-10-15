
# Cross-Modal TF→Gene Regulatory Inference & Cross-Species Comparison

*A practical summary for integrating RNA-based GRN inference (GENIE3/GRNBoost2) with scATAC motif evidence, peak–gene linking, and extending to human–mouse comparisons.*

**Author:** (auto-generated for Jussi Kupari)  
**Date:** 2025-10-15

---

## Table of Contents
1. [Problem Framing](#problem-framing)
2. [Linking Peaks to Genes: TSS Windows vs Co‑accessibility](#linking-peaks-to-genes)
3. [Core Integration Concept (RNA GRN + ATAC Motifs)](#core-integration-concept)
4. [Scoring: Motif Evidence Score (MES) and Combined Regulatory Confidence (CRC)](#scoring)
5. [Python-Oriented Implementation Plan](#implementation-plan)
6. [Data Alignment Across Datasets](#alignment)
7. [Cross-Species (Human–Mouse) Comparison](#cross-species)
8. [Motif-Level Comparisons](#motif-level)
9. [Hypothesizing Mechanisms for Species-Biased Expression](#mechanism)
10. [Diagnostics, Visualizations, and Deliverables](#viz)
11. [Practical Tips & Gotchas](#tips)
12. [Next Steps & Options](#next-steps)

---

## 1) Problem Framing <a name="problem-framing"></a>
We want to identify **Transcription Factors (TFs)** that regulate **specific genes**, by integrating:
- **RNA-based GRN inference** (e.g., **GENIE3/GRNBoost2**) → TF→gene edges with weights (correlation/importance).
- **scATAC-seq evidence** → TF **motif presence** and **accessibility** in peaks **linked to genes**.

Goal: For a predicted TF→gene edge, ask *“Is the TF’s motif found in accessible regulatory peaks linked to that gene (promoter or distal) in the relevant cell type?”* Then combine the evidence into a ranked regulatory network.

---

## 2) Linking Peaks to Genes: TSS Windows vs Co‑accessibility <a name="linking-peaks-to-genes"></a>
Two common strategies:

**A. Fixed windows around TSS (e.g., ±10 kb; promoter-biased)**
- *Pros:* Simple, fast, good for promoter-driven regulation.
- *Cons:* Misses **distal enhancers**; equates linear proximity with regulation; blind to 3D chromatin.

**B. Co-accessibility / correlation-based linking (Cicero, Signac `LinkPeaks`, ArchR)**
- *Pros:* Captures **distal enhancer–promoter** interactions via cell-to-cell co-variation; more realistic for complex regulation.
- *Cons:* Heavier compute; requires sufficient cell numbers; can reflect compositional confounders.

**Recommendation for TF–gene validation:** Prefer **co‑accessibility-based links** (B) for distal coverage; optionally include promoter windows as baseline.

---

## 3) Core Integration Concept (RNA GRN + ATAC Motifs) <a name="core-integration-concept"></a>
1. From RNA (GENIE3/GRNBoost2): ranked **TF→gene** edges with weights.
2. From ATAC: **peaks linked to genes** + **motif annotations** (TF motif present?), with **accessibility** and optional **co‑accessibility** weights.
3. For each edge (TF, gene, cell type): aggregate motif evidence across the gene’s linked peaks → produce a **Motif Evidence Score (MES)**.
4. Combine **GENIE3/GRN** weight and **MES** into a **Combined Regulatory Confidence (CRC)**.

---

## 4) Scoring: Motif Evidence Score (MES) and Combined Regulatory Confidence (CRC) <a name="scoring"></a>
Let \(P(g)\) be the set of peaks linked to gene \(g\). For TF \(t\) and cell type \(c\):

**Motif Evidence Score (MES)**

\[
\mathrm{MES}_{t,g,c} 
= \sum_{p \in P(g)} \underbrace{\mathbf{1}[\text{motif}(t) \in p]}_{\text{motif presence}} 
\times \underbrace{A_{p,c}}_{\text{accessibility weight}} 
\times \underbrace{L_{p,g}}_{\text{link/co-accessibility score}} 
\times \underbrace{M_{t,p}}_{\text{(optional) motif strength}}
\]

**Combined Regulatory Confidence (CRC)**

Normalize within cell type (z-scores) and combine:

\[
\mathrm{CRC}_{t,g,c} = z\big(\mathrm{GENIE3\ weight}_{t,g}\big) \times z\big(1 + \mathrm{MES}_{t,g,c}\big)
\]

- Product emphasizes edges supported by **both** RNA and ATAC layers.
- Variants: weighted sums, inclusion of **TF activity** (chromVAR deviations) as a multiplicative factor.

---

## 5) Python-Oriented Implementation Plan <a name="implementation-plan"></a>
**Inputs**
- **scRNA:** expression matrix; TF list; cell-type labels; (optional) pseudobulk per cell type.
- **scATAC:** peak×cell matrix; peak metadata; **peak→gene links** (Cicero/Signac/ArchR or distance-based backup); **motif hits** per peak (pycistarget/gimmemotifs/MOODS); (optional) chromVAR deviations.

**Workflow**
1. **GRN inference (RNA):** run `arboreto` (**GRNBoost2**/**GENIE3**) to get TF→gene edges with importance.
2. **Peak→gene linking (ATAC):** import Cicero/Signac links; if absent, use TSS windows with distance-decay.
3. **Motif annotation:** annotate peaks with TF motifs (prefer **pycistarget** with cisTarget DBs; or scan sequences).
4. **Pseudobulk accessibility:** sum per cell type and normalize (CPM) for stability across datasets.
5. **Compute MES per edge per cell type:** aggregate motif+accessibility+link evidence across linked peaks.
6. **Compute CRC:** combine z-scored GENIE3 weight and MES (and optionally TF activity).
7. **Filter & rank:** retain edges with consistent TF/gene expression & accessible peaks in the cell type.
8. **Export:** ranked edges, plus per-edge diagnostics.

**Notes for your stack**
- You prefer **Python tools** (arboreto/GRNBoost2 + **pycistarget**). Rcistopic topic regions can seed motif analyses.
- You can skip re-running arboreto if you already have **GENIE3** edges from a prior SCENIC run—just load and proceed to MES/CRC.

---

## 6) Data Alignment Across Datasets <a name="alignment"></a>
Because RNA and ATAC come from **different experiments**, align at the **cell-type level**:
- Harmonize labels (ontology mapping, marker-based matching, or label transfer).
- Create **pseudobulk RNA** and **pseudobulk ATAC** per cell type.
- Restrict inference to cell types where **TF** and **target gene** are expressed and the **linked peaks are accessible**.

---

## 7) Cross-Species (Human–Mouse) Comparison <a name="cross-species"></a>
To compare TF regulation for orthologous genes between human and mouse:
1. **Orthology mapping:** map genes and TFs to 1:1 orthologs (e.g., Ensembl BioMart / HomoloGene).
2. **Cell-type alignment:** ensure comparable cell-type labels; if not exact, use coarser ontology.
3. **Run pipeline per species:** compute **CRC** for each TF→gene per cell type.
4. **Compare per gene-of-interest:**
   - Top TFs (ranked by CRC) across species.
   - Overlap metrics (e.g., Jaccard of top-N TFs), rank correlation (Spearman).
   - Presence/absence of **motif evidence** in linked peaks in each species.
5. **Optional positional conservation:** liftOver peaks to check syntenic motif conservation (nice-to-have, not required for functional comparison).

---

## 8) Motif-Level Comparisons <a name="motif-level"></a>
- Motifs are represented as **PWMs** from databases like **JASPAR/CIS‑BP/HOCOMOCO**. They capture a TF’s sequence preference.
- **Direct sequence match ≠ binding**: chromatin context and cofactors matter; treat motif hits as necessary but not sufficient evidence.
- For cross-species work, operate at the level of **motif families/TF families** (e.g., MEF2, ETS) to avoid over-penalizing small PWM differences.
- Practical motif conservation metric for a gene:

\[
\mathrm{Motif\ Conservation}(g) = \frac{|\mathcal{F}_\text{human}(g) \cap \mathcal{F}_\text{mouse}(g)|}{|\mathcal{F}_\text{human}(g) \cup \mathcal{F}_\text{mouse}(g)|}
\]

where \(\mathcal{F}(g)\) is the set of TF families with motif hits in peaks linked to gene \(g\).

---

## 9) Hypothesizing Mechanisms for Species‑Biased Expression <a name="mechanism"></a>
When a gene is **broader** in one species than the other:
1. **CRC deltas:** identify TFs with high CRC in the broad species and low in the restricted species.
2. **Motif presence deltas:** check whether those TF motifs exist in gene-linked peaks in the broad species but are absent/weak in the other → **enhancer gain/loss** hypothesis.
3. **Accessibility deltas:** inspect whether the relevant peaks are accessible across more cell types in the broad species → **cis** changes.
4. **TF activity deltas:** if motifs/accessibility are similar but TF activity (chromVAR) differs → **trans** changes.
5. **Synthesize:** propose candidate mechanisms (e.g., “gain of distal MEF2 sites drives broader expression in species A”).

**A simple hypothesis score per TF:**

\[
\mathrm{HypothesisScore}_{t,g} = \Delta \mathrm{CRC} + w_1\, \Delta \mathrm{MotifPresence} + w_2\, \Delta \mathrm{Accessibility}
\]

with weights \(w_1, w_2\) chosen for balance.

---

## 10) Diagnostics, Visualizations, and Deliverables <a name="viz"></a>
**Per gene:**
- Side-by-side **bar plots** of top TF CRCs (human vs mouse, per aligned cell type).
- **Genome-style track**: linked peaks, motif hits, accessibility (per cell type), and co-accessibility arcs.
- **Heatmaps**: peak accessibility across cell types; motif family presence across species.

**Per TF:**
- Top targets by CRC per species/cell type; GO/pathway enrichment of shared vs species-specific targets.

**Deliverables:**
- `crc_ranked_edges.tsv`: `TF, gene, cell_type, CRC, GENIE3_weight, MES`.
- Edge annotations: counts of motif+linked peaks; promoter vs distal breakdown; optional TF activity deviations.
- Cross-species summary tables (rank correlations, Jaccard overlaps).
- Notebook(s) with plotting utilities.

---

## 11) Practical Tips & Gotchas <a name="tips"></a>
- **Cell-type alignment** is crucial; consider coarsening labels to ensure comparability.
- **Motif redundancy**: many TFs share motifs; map to **TF families** and retain multiple candidate TFs unless expression/activity disambiguates.
- **Peak→gene links**: Prefer Cicero/Signac/ArchR; if using TSS windows, apply **distance-decay** and require promoter support for “direct” edges.
- **Pseudobulk** stabilizes signals across experiments/datasets.
- **Multiple testing**: if thresholding motif hits/enrichments, control FDR at motif or edge level.
- **Replicates**: if available, meta-analyze CRC/MES across replicates to avoid overfitting.

---

## 12) Next Steps & Options <a name="next-steps"></a>
- Build a **notebook scaffold** that:
  1) ingests your GENIE3 edges and Rcistopic/Signac outputs,  
  2) computes MES & CRC per cell type,  
  3) performs cross-species ortholog mapping,  
  4) outputs ranked edges + diagnostic plots.
- Add **chromVAR TF activity** to boost specificity and help disambiguate motif families.
- Implement **liftOver**-assisted motif conservation as an optional confidence layer.

---

### Minimal Pseudocode Snippets

**MES aggregation**
```python
MES[t,g,c] = sum_over_linked_peaks(
    motif_present(t, p) * accessibility[p,c] * link_score[p,g] * motif_strength(t,p)
)
```

**CRC combination**
```python
CRC[t,g,c] = z(GENIE3_weight[t,g]) * z(1 + MES[t,g,c])  # optionally * z(TF_activity[t,c])
```

**Cross-species overlap**
```python
overlap = jaccard(top_TFs_speciesA(g,c), top_TFs_speciesB(g,c))
rank_r  = spearmanr(rank_speciesA(g,c), rank_speciesB(g,c))
```

---

### Personalization Notes
- Tailored to **Python-first** tooling (arboreto/GRNBoost2 + **pycistarget**), consistent with your preference.
- Compatible with **Rcistopic** outputs (use topics/regions for motif enrichment and to guide linked-peak sets).

---

*End of document.*
