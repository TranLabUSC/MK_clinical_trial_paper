# MK Analysis Scripts - Single-Cell RNA-seq Analysis

## Overview
This repository contains R scripts for analyzing single-cell RNA-sequencing data from glioblastoma (GBM) patients, focusing on T cell populations, clonal expansion, and immune dynamics.

---

## Documented Scripts

### ✅ Annotation/annotation.R
**Purpose**: Cell type annotation using marker gene expression  
**Manuscript Figures**: Figure 5a (All cells UMAP), Figure 6a (T cells UMAP)  
**Key Features**:
- Tests clustering resolutions (0.1, 0.3, 1, 3) to find optimal granularity
- Plots marker genes for manual cell type annotation
- Isolates and reclusters T cells for detailed subpopulation analysis
- Removes low-quality/contaminating clusters
## Clonal Expansion Analysis Pipeline (Figure 6d)

### ✅ clonal_expansion_analysis_1.R - Step 1: Data Preparation
**Purpose**: Prepare clonotype count and proportion tables for all T cell subpopulations  
**Pipeline Position**: Step 1 of 4  
**Key Features**:
- Processes TCR VDJ data to extract CDR3β sequences
- Maps T cell clusters to cell type annotations  
- Generates absolute clone count tables (clones × samples)
- Generates proportional clone frequency tables (normalized per sample)
- Filters out low-quality samples (<1 cell)

**Inputs**: T cell Seurat object, TCR contig annotations  
**Outputs**: `clonotype_df_absolute_*.csv`, `clonotype_df_proportion_*.csv` (per cell type)  
**Next Step**: diversity_calculation_2.R
### ✅ diversity_calculation_2.R - Step 2: Diversity Indices
**Purpose**: Calculate Shannon and Simpson diversity indices to quantify clonal expansion  
**Pipeline Position**: Step 2 of 4  
**Key Features**:
- Calculates Shannon diversity index (richness + evenness)
- Calculates Simpson diversity index (emphasizes dominant clones)
- Processes all T cell subpopulations from Step 1

**Biological Interpretation**:
- **LOW diversity** = HIGH expansion = STRONG antigen specificity
- **HIGH diversity** = LOW expansion = WEAK/BROAD response

**Inputs**: `clonotype_df_proportion_*.csv` (from Step 1)  
**Outputs**: `shannon_clonal_diversity_*.txt`, `simpson_clonal_diversity_*.txt`  
**Next Step**: clonal_expansion_calculation_3.R
### ✅ clonal_expansion_calculation_3.R - Step 3: Expansion Metrics
**Purpose**: Calculate clonal expansion by comparing diversity indices between timepoints  
**Pipeline Position**: Step 3 of 4  
**Key Features**:
- Compares diversity indices between two timepoints (C1 vs C2)
- Calculates expansion using ratio method (diversity_C2 / diversity_C1)
- Processes all T cell subpopulations with both Shannon and Simpson indices

**Biological Interpretation**:
- **Ratio < 1**: Clonal expansion occurred (diversity decreased)
- **Ratio > 1**: Clonal diversification occurred (diversity increased)
- Lower ratio = Stronger antigen-specific response

**Inputs**: `shannon_clonal_diversity_*.txt`, `simpson_clonal_diversity_*.txt` (from Step 2)  
**Outputs**: `clonal_expansion_{shannon|simpson}_{celltype}_T_cells_division/C1_vs_C2.txt`  
### ✅ clonal_expansion_plots_4.R - Step 4: Visualization (Figure 6d)
**Purpose**: Generate clonal expansion visualizations across T cell subpopulations, timepoints, and patient cohorts  
**Pipeline Position**: Step 4 of 4 - FINAL STEP  
**Manuscript Figure**: **Figure 6d**
**Key Features**:
- Creates expansion ratio plots (diversity ratio vs 1.0)
- Generates paired diversity boxplots with connecting lines
- Compares multiple timepoint pairs (Pre vs C1, Pre vs C2, C1 vs C2)
- Stratifies by patient cohort (control, short-term, long-term survivors)
- Performs statistical testing (t-test, paired t-test)

**Visualization Types**:
1. Expansion ratio plots: Violin/boxplots comparing ratios to 1.0
2. Paired diversity boxplots: Actual diversity values with patient trajectories

**Inputs**: Expansion ratios (Step 3), diversity indices (Step 2), clinical metadata  
**Outputs**: PDF plots organized by timepoint comparison and cell type

**Next Step**: clonal_expansion_plots_4.R (generates Figure 6d)



- Exports data for Loupe Browser visualization

**Inputs**: Seurat objects, markers.csv  
**Outputs**: Feature plots, UMAP coordinates, cleaned Seurat objects

---

### UMAP_comparisons.R
**Purpose:** Compare chemotaxis gene expression in Non-Classical Monocytes between Pre and C1 timepoints in experimental cohort patients.

**Biological Context:**
- Non-Classical Monocytes: CD14+CD16++ monocyte subset with patrolling behavior and immune surveillance function
- Chemotaxis: Directed cell migration in response to chemical gradients, critical for immune cell recruitment to tumor sites
- Analysis compares baseline (Pre-treatment) vs Cycle 1 (C1) to assess treatment impact on immune cell migration potential

**Workflow:**
1. Load Non-Classical Monocyte Seurat object and filter patient-specific clusters
2. Separate into control (patients 1, 4, 8, 9) vs experimental cohorts (12 patients)
3. Downsample to equal cell numbers between Pre and C1 for balanced comparisons
4. Generate UMAP visualizations for both cohorts showing timepoint distribution
5. Filter experimental cohort cells within UMAP boundaries to remove outliers
6. Define comprehensive chemotaxis gene list (400+ genes: chemokines, receptors, adhesion molecules, migration regulators)
7. Create FeaturePlots showing chemotaxis gene expression patterns (Pre vs C1)
8. Generate interleukin (IL) gene expression FeaturePlots as additional analysis

**Key Parameters:**
- Retained clusters: 0, 1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16 (excludes patient-specific/artifact clusters)
- UMAP limits: X [-7, 7], Y [-6, 8] for standardization
- Expression quantile cutoffs: q10 (min) to q90 (max) for color scaling
- Downsampling ensures balanced Pre/C1 comparison without timepoint bias

**Outputs:**
- `NC_Mono_Control_Pre_C1.pdf`: Control cohort UMAP (Pre=blue, C1=red)
- `NC_Mono_Experiment_Pre_C1.pdf`: Experimental cohort UMAP (Pre=blue, C1=red)
- `chemotaxis_genes.pdf`: 400+ chemotaxis genes split by timepoint (**Manuscript Figure 5d**)
- `IL_genes_test.pdf`: Interleukin genes split by timepoint

**Manuscript Figure:**
- **Figure 5d**: Chemotaxis gene expression comparison in experimental cohort NC monocytes (Pre vs C1)

---

### NC_Monocyte_Optimal_Transport.ipynb
**Purpose:** Employ optimal transport theory to predict cell fate transitions in Non-Classical Monocytes from Pre-treatment to C1, identifying which Pre cells are predecessors to which C1 cells.

**Biological Context:**
- **Cell Fate Prediction**: Inferring developmental/transitional relationships between cells at different timepoints using computational methods
- **Optimal Transport**: Mathematical framework that finds the most efficient way to transform one cell distribution (Pre) into another (C1)
- **Non-Classical Monocytes**: CD14+CD16++ monocyte subset - analysis reveals how treatment alters population dynamics and spatial organization
- **Clinical Relevance**: Understanding how immune cell populations reorganize in response to MK-3475 + MLA treatment

**Mathematical Concepts:**
- **Earth Mover's Distance (EMD)**: Also known as Wasserstein-1 distance, quantifies minimum "cost" to transform source distribution to target distribution
- **Transport Plan**: Matrix T where T[i,j] represents probability/mass moved from source cell i to target cell j
- **Cost Matrix**: Pairwise Euclidean distances between all Pre and C1 cells in UMAP space
- **Predecessor Mapping**: For each C1 cell, identifies the Pre cell with maximum transport probability as its most likely predecessor

**Workflow:**

**PART 1: Experimental Cohort Analysis (Figure 5b)**
1. Load experimental cohort AnnData object (Pre and C1 timepoints)
2. Map timepoints to numeric values (Pre=0, C1=1) and separate cells
3. Extract UMAP coordinates for Pre (source) and C1 (target) cells
4. Compute cost matrix using Euclidean distances in UMAP space
5. Solve optimal transport problem to obtain transport plan matrix
6. Calculate Earth Mover's Distance (EMD) to quantify overall population shift
7. Sample 2400 cells and generate arrow plot showing individual fate transitions
8. Create kernel density estimation (KDE) plot overlaying Pre (blue) and C1 (red) distributions
9. Extract cell barcodes and create predecessor mapping DataFrame
10. Save predecessor mapping to CSV for downstream cluster tracking analysis

**PART 2: Control Cohort Analysis (Comparison)**
11. Load control cohort data (patients 1 and 4 only)
12. Repeat entire optimal transport workflow for control cohort
13. Sample 823 cells for arrow visualization
14. Generate density plot with different colors (Pre=green, C1=red) to distinguish from experimental

**Key Parameters:**
- Distance metric: Euclidean distance in 2D UMAP space
- Transport solver: Earth Mover's Distance (exact OT solver from POT library)
- Arrow plot samples: 2400 cells (experimental), 823 cells (control)
- Visualization: Kernel density estimation with transparency for overlapping distributions

**Outputs:**
- **Arrow Plots**: Visualize individual cell transition trajectories from Pre to C1
- **Density Shift Plots**:
  - `NC_Mono_Exp_Pre_C1_UMAP_Density_Shift.pdf`: Experimental cohort (Pre=blue, C1=red) - **Manuscript Figure 5b**
  - `NC_Mono_Control_Pre_C1_UMAP_Density_Shift.pdf`: Control cohort (Pre=green, C1=red)
- **Predecessor Mapping**: `predecessor_mapping.csv` (maps each C1 cell barcode to its predicted Pre predecessor barcode)
- **Hexbin Density Plots**: Alternative visualization showing spatial density distributions

**Manuscript Figure:**
- **Figure 5b**: UMAP density shift plot showing Non-Classical Monocyte population movement from Pre to C1 in experimental cohort

**Biological Interpretation:**
- Larger EMD value = Greater population shift between timepoints
- Arrow directions reveal dominant cell fate transition patterns
- Density changes show regions of cell accumulation (red increases) or depletion (blue decreases)
- Predecessor mapping enables tracking how specific Pre clusters transform into C1 clusters

**Downstream Applications:**
The predecessor mapping CSV is used in subsequent analyses to:
- Track which Pre clusters give rise to specific C1 clusters
- Identify genes associated with particular fate transitions
- Compare transition patterns between experimental and control cohorts
- Correlate cell fate changes with clinical outcomes

**Dependencies:**
- Data: `MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.h5ad` (experimental), control cohort h5ad
- Libraries: scanpy, anndata, numpy, pandas, POT (Python Optimal Transport), matplotlib, seaborn

---

### nc_monocytes_analysis.R

**Purpose:** Comprehensive differential expression and survival analysis on Non-Classical Monocytes using predecessor relationships from optimal transport, generating three key manuscript figures.

**Manuscript Figures Generated:**
- **Figure 5c**: Heatmap with hierarchical clustering of differentially expressed genes
- **Figure 5e**: GSEA (Gene Set Enrichment Analysis) results
- **Figure 5f**: Cox proportional hazards survival analysis

**Workflow Overview:**

**SECTION 1-2: Data Preparation and Single-Cell DE (Lines 1-287)**
1. Define `plot_predecessor_distribution()` function to visualize optimal transport predecessor patterns
2. Load predecessor mapping CSV from optimal transport notebook
3. Define `perform_de_between_cluster_and_predecessors()` for differential expression
   - Uses Wilcoxon rank-sum test (non-parametric, robust to outliers)
   - Compares C1 cells vs their specific Pre predecessors identified by optimal transport
4. Execute DE analysis for cluster 1 and generate initial visualizations

**SECTION 3-4: Visualization (Lines 310-385)**
5. Generate heatmaps of top differentially expressed genes (upregulated vs downregulated)
6. Create volcano plots showing fold change vs significance
7. Set significance thresholds (adj p-value < 0.05, |log2FC| > 0.25)

**SECTION 5: Pseudo-bulk RNA-seq Aggregation (Lines 388-516)**
8. Load NC monocyte Seurat object (experimental cohort, Pre and C1 timepoints only)
9. Aggregate single-cell counts by Patient × TimePoint (pseudo-bulk approach)
   - Improves statistical power by pooling biological replicates within patients
   - Reduces technical noise inherent in single-cell data
10. Create DGEList object for edgeR analysis
11. Apply TMM normalization (Trimmed Mean of M-values) to account for library size differences
12. Estimate dispersion parameters (common, trended, tagwise) using edgeR framework
13. Fit generalized linear model (GLM) with design matrix
14. Perform likelihood ratio test (LRT) comparing C1 vs Pre timepoints
15. Filter significant genes (FDR < 0.05, |logFC| > 1)

**SECTION 6: Log2 Fold Change Normalization (Lines 520-649)**
16. Calculate patient-specific log2FC values (C1 vs Pre baseline)
    - Pre timepoint serves as within-patient baseline reference
    - Paired design controls for patient-specific baseline differences
17. Add pseudocount (+1) to prevent log(0) errors in low-expression genes
18. Standardize using z-score transformation
    - Formula: z = (x - mean) / SD
    - Enables cross-patient comparison despite different baseline expression levels

**SECTION 7-8: Heatmap Generation - Figure 5c (Lines 653-1047)**
19. Perform hierarchical clustering on genes (Euclidean distance, complete linkage)
    - Identifies co-regulated gene modules
    - Groups genes with similar expression patterns across patients
20. Cut dendrogram into k=20 clusters for biological interpretation
21. Order patients by survival outcomes (short-term → long-term survivors)
22. Generate color-coded heatmap with z-scored log2FC values
    - Red/yellow: Upregulated in C1 (positive z-score)
    - Blue: Downregulated in C1 (negative z-score)
23. Add annotation bars showing patient survival groups
24. **Save as Figure 5c** (manuscript-ready heatmap)

**SECTION 9: GSEA Rank File Preparation - Figure 5e (Lines 1055-1265)**
25. Calculate z-scores for each gene comparing long-term vs short-term survivors
    - Formula: z = (mean_long - mean_short) / SE_diff
    - Positive z-score: Gene expression beneficial for survival (higher in long-term survivors)
26. Compute t-values using Welch's t-test (accounts for unequal variances between groups)
27. Rank genes by z-score in descending order
28. Export ranked gene list in .rnk format for GSEA software
29. GSEA tool uses this ranked list to identify enriched biological pathways in long-term survivors
30. **Results generate Figure 5e** (pathway enrichment visualization)

**SECTION 10: Cox Proportional Hazards Survival Analysis - Figure 5f (Lines 1238-1445)**
31. Calculate pathway activity scores for each patient
    - Mean z-scored log2FC across all genes within specific biological pathways
    - Represents overall pathway upregulation/downregulation in C1 vs Pre
32. Stratify patients by median pathway activity (binary split: high vs low activity)
33. Fit Cox proportional hazards regression model for each pathway
    - Hazard Ratio (HR) > 1: Higher pathway activity → worse survival prognosis
    - Hazard Ratio (HR) < 1: Higher pathway activity → better survival prognosis
34. Calculate p-values using likelihood ratio test
35. Generate visualization showing most significant pathways and their hazard ratios
36. **Save as Figure 5f** (survival curves and forest plots)

**Statistical Methods:**
- **Wilcoxon rank-sum test**: Non-parametric differential expression for single-cell data
- **edgeR GLM framework**: Robust pseudo-bulk RNA-seq analysis with TMM normalization and dispersion estimation
- **Z-score standardization**: Enables cross-patient comparison by removing baseline differences
- **Hierarchical clustering**: Identifies co-regulated gene modules using complete linkage
- **Welch's t-test**: Accounts for unequal variances between survivor groups in GSEA preparation
- **Cox regression**: Estimates hazard ratios linking pathway activities to survival outcomes

**Biological Interpretation:**
- Genes upregulated in C1 (vs Pre predecessors) suggest treatment-induced transcriptional changes
- Gene modules from clustering reveal coordinated biological processes activated by treatment
- Pathway enrichment in long-term survivors identifies protective molecular mechanisms
- Cox analysis quantifies how specific pathway activities predict patient survival

**Key Outputs:**
- Differential expression tables (single-cell and pseudo-bulk approaches)
- Gene cluster assignments (k=20 co-regulated modules)
- GSEA rank file (.rnk format) for pathway enrichment analysis
- **Figure 5c**: Heatmap PDF showing hierarchical clustering of treatment-induced gene changes
- **Figure 5e**: GSEA pathway enrichment visualization
- **Figure 5f**: Cox survival analysis plots (forest plots and Kaplan-Meier curves)

**Dependencies:**
- Seurat object: `MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS`
- Predecessor mapping: `predecessor_mapping.csv` (from NC_Monocyte_Optimal_Transport.ipynb)
- Libraries: Seurat, edgeR, dplyr, tidyr, ggplot2, pheatmap, survival

---

### pathway_activity_comparison_new.R

**Purpose:** Calculate pathway activity scores and compare them between treatment cohorts (MK-3475 Alone vs MK-3475 + MLA) across T cell subpopulations and timepoints.

**Manuscript Figures Generated:**
- **Figure 6b**: Pathway signal change line plots showing patient trajectories between timepoints
- **Figure 6c**: Boxplots comparing pathway signal changes between treatment arms

**Core Concept - Pathway Activity:**
Pathway activity represents the overall functional state of a biological pathway by aggregating gene expression across all pathway members. Instead of examining individual genes, pathway activity captures coordinated biological processes:
- **Calculation**: Mean expression of all genes in the pathway for each cell
- **Biological Rationale**: Reduces noise from individual gene variability, captures true biological signal
- **Interpretation**: Higher pathway activity = Enhanced functional state (e.g., more T cell activation)

**Workflow Overview:**

**PART 1: Function Definitions**

**1. `create_survival_data()` Function (Lines 10-114)**
- Reads GMT file to extract genes for specified pathway
- Calculates pathway activity score for each cell (mean expression of pathway genes)
- Aggregates pathway activity by patient and timepoint
- Updates survival dataframe with pathway activity scores for all timepoints (Pre, C1, C2, C4, C6, C9, C18, C36)

**2. `get_cluster_proportion()` Function (Lines 119-157)**
- Calculates the proportion of a specific T cell subset at a given timepoint
- Formula: Cluster Proportion = (Cells in cluster at timepoint) / (Total cells at timepoint)
- Accounts for inter-patient variability in T cell subset frequencies
- Essential for weighted pathway signal calculation

**PART 2: Main Analysis Pipeline (Two Implementations)**

**Section 1: Original Implementation (Lines 161-426)**

**Setup and Parameters:**
3. Load T cell Seurat object with all subpopulations
4. Load survival data (UF patients only)
5. Define 10 immune pathways from GMT file:
   - GO_ADAPTIVE_IMMUNE_RESPONSE
   - GO_T_CELL_ACTIVATION
   - GO_REGULATION_OF_IMMUNE_RESPONSE
   - GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE
   - GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS
   - REACTOME_IMMUNE_SYSTEM
   - GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY
   - GO_ACTIVATION_OF_IMMUNE_RESPONSE
   - GO_REGULATION_OF_T_CELL_ACTIVATION
   - GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE

6. Map 13 T cell subtypes to Seurat clusters:
   - Activated_CD4 (cluster 0)
   - Effector_CD8 (cluster 1)
   - Effector_Memory_Precursor_CD8 (cluster 2)
   - Exhausted_T (cluster 3)
   - Gamma_Delta_T (cluster 4)
   - Active_CD4 (cluster 5)
   - Naive_CD4 (clusters 6, 9, 18)
   - Memory_CD4 (cluster 7)
   - Stem_Like_CD8 (cluster 8)
   - Effector_Memory_CD8 (cluster 10)
   - Central_Memory_CD8 (cluster 12)
   - GZMK_Effector_Memory_CD8 (cluster 13)
   - Proliferating_Effector (clusters 14, 16, 17)

**Analysis Loop (For Each Pathway × Cell Type × Timepoint Pair):**

7. Subset Seurat object for specific T cell subtype
8. Normalize counts using Relative Counts (RC) method with scaling factor 1e6
9. Calculate pathway activity scores using `create_survival_data()`
10. Filter out IDH-positive patients

**Pathway Signal Calculation with Cluster Proportion Weighting:**

11. Calculate cluster proportions for both timepoints (e.g., Pre and C1)
12. Calculate weighted pathway signal:
    - **Pathway Signal = Mean Expression × Cluster Proportion**
    - Accounts for how abundant this cell type is at each timepoint
    - Higher proportion = Greater contribution to overall immune response

13. Calculate signal change between timepoints:
    - **Ratio method**: Signal_C1 / Signal_Pre (fold-change)
    - **Difference method**: Signal_C1 - Signal_Pre (absolute change)

**Statistical Testing and Visualization:**

14. Perform Wilcoxon rank-sum test comparing MK-3475 Alone vs MK-3475 + MLA cohorts
15. Calculate fold change and log2 fold change between cohort medians
16. Generate line plots showing individual patient trajectories (**Figure 6b** component)
    - Each line represents one patient's pathway signal change
    - Color-coded by treatment arm (red = MK-3475 + MLA, blue = MK-3475 Alone)
17. Generate boxplots comparing signal change distributions (**Figure 6c** component)
    - Shows overall cohort differences
    - Includes p-value annotations
    - Jittered points show individual patient values
18. Save combined grid plots (line plots + boxplots) as PDF files
19. Export comparison statistics to CSV (pathway, celltype, comparison, FC, logFC, p-value)

**Section 2: Refactored Implementation (Lines 432-862)**

**Enhanced Features:**
20. `get_significance()` function converts p-values to publication symbols:
    - p < 0.001: "***"
    - p < 0.01: "**"
    - p < 0.05: "*"
    - p ≥ 0.05: "ns" (not significant)

21. Improved code organization using dplyr pipes for data manipulation
22. Enhanced plotting with removed x-axis labels and legends for cleaner visualizations
23. Better error handling and significance tracking
24. Expanded output CSV includes significance column

**Biological Interpretation:**

**Pathway Activity Changes:**
- **Increased pathway activity** (C1 vs Pre) suggests treatment-induced immune activation
- **Higher activity in MK-3475 + MLA** indicates combination therapy enhances immune pathways
- **Cell type specificity** reveals which T cell subsets drive the immune response

**Cluster Proportion Weighting Importance:**
- Accounts for inter-patient variability in immune cell composition
- Patient A might have 20% Effector CD8, Patient B only 5%
- Weighting ensures fair comparison despite compositional differences
- More accurate representation of pathway's overall contribution to immune response

**Treatment Cohort Comparison:**
- **MK-3475 Alone** (control): Anti-PD-1 monotherapy
- **MK-3475 + MLA** (experimental): Anti-PD-1 + MLA combination
- Significant differences indicate combination therapy's enhanced immunomodulatory effects

**Statistical Framework:**
- **Wilcoxon rank-sum test**: Non-parametric test comparing two independent groups
  - Robust to outliers and non-normal distributions
  - Tests whether cohort distributions differ significantly
- **Fold change analysis**: Quantifies magnitude of difference between cohorts
  - FC > 1: Experimental arm has higher pathway activity
  - FC < 1: Control arm has higher pathway activity

**Key Outputs:**
- **Figures 6b and 6c**: Combined grid plots for each pathway × cell type combination
  - Top row: Line plots showing patient trajectories (Figure 6b style)
  - Bottom row: Boxplots comparing cohort distributions (Figure 6c style)
- **CSV results table**: Comprehensive statistics (pathway, cell type, timepoint comparison, FC, logFC, p-value, significance)
- **Output directory structure**: Organized by method (ratio/difference) and cluster proportion inclusion

**Parameters User Can Modify:**
- `include_cluster_proportion`: Toggle cluster proportion weighting (TRUE/FALSE)
- `method`: Choose signal change calculation ("ratio" or "difference")
- `pathways`: List of pathways to analyze
- `timepoints`: Which timepoint pairs to compare (default: Pre vs C1, Pre vs C2, C1 vs C2)

**Dependencies:**
- T cell Seurat object: `MK_T_Cells_seurat_obj.RDS`
- Survival data: `patient_survival_data.csv`
- GMT pathway file: `subsubset_output.gmt` (Gene Matrix Transposed format with pathway definitions)
- Libraries: Seurat, dplyr, tidyr, ggplot2, cowplot

---

### cluster_proportion_comparison_new.R

**Purpose:** Compare cluster proportions of major cell types and T cell subpopulations between treatment cohorts (MK-3475 Alone vs MK-3475 + MLA) across timepoints.

**Supplemental Tables Generated:**
- **Supplemental Table S3**: Cluster proportion comparison results for 5 major cell types
- **Supplemental Table S4**: Cluster proportion comparison results for 13 T cell subpopulations

**Core Concept - Cluster Proportion:**
Cluster proportion represents the relative abundance of a specific cell type within the total cell population:
- **Calculation**: (Number of cells in cluster) / (Total cells at timepoint)
- **Biological Meaning**: Reflects immune cell composition and dynamic changes in response to treatment
- **Interpretation**:
  - Increased proportion = Cell type expansion (proliferation or recruitment)
  - Decreased proportion = Cell type contraction (depletion or differentiation)

**Script Structure:**

**SECTION 1: Major Cell Types Analysis (Lines 1-303) → Supplemental Table S3**

**Cell Types Analyzed (5 major immune populations):**
1. **Classical Monocytes** (Clusters 1, 3, 16, 12): CD14++ monocytes, inflammatory response
2. **Non-Classical Monocytes** (Clusters 9, 33): CD14+CD16++ patrolling monocytes
3. **cDC** (Cluster 31): Conventional dendritic cells, antigen presentation
4. **pDC** (Cluster 36): Plasmacytoid dendritic cells, interferon production
5. **NK** (Cluster 0): Natural killer cells, innate cytotoxicity

**Workflow:**

1. **Function Definitions:**
   - `get_cluster_proportion()`: Calculates proportion for each patient at specified timepoint
     - Filters cells by cluster and timepoint
     - Computes proportion per patient: cells_in_cluster / total_cells
     - Returns patient-level proportion dataframe
   - `get_significance()`: Converts p-values to publication symbols (*, **, ***, ns)

2. **Data Loading and Setup:**
   - Load all-cells Seurat object (`MK_Cells_seurat_obj.RDS`)
   - Load survival data with treatment arm assignments (UF patients only)
   - Define timepoint pairs to compare (Pre vs C1, Pre vs C2, C1 vs C2)

3. **Analysis Loop (For Each Cell Type × Timepoint Pair):**
   
   **Proportion Calculation:**
   - Calculate cluster proportions for both timepoints (e.g., Pre and C1)
   - Merge with survival data to get treatment arm information
   - Compute proportion change:
     - **Ratio method**: Proportion_B / Proportion_A (fold-change)
     - **Difference method**: Proportion_B - Proportion_A (absolute change)
   
   **Statistical Testing:**
   - Calculate log2 fold change (median across all patients):
     - Formula: median(log2((Proportion_B + ε) / (Proportion_A + ε)))
     - Positive value = Cell type expansion at timepoint B
     - Negative value = Cell type contraction at timepoint B
   - Perform Wilcoxon rank-sum test comparing treatment arms:
     - Tests if proportion changes differ between MK-3475 Alone vs MK-3475 + MLA
     - Requires ≥2 patients per arm for valid statistical comparison
   - Determine significance level (p < 0.001: ***, p < 0.01: **, p < 0.05: *, p ≥ 0.05: ns)
   
   **Visualization Generation:**
   - **Line Plots**: Individual patient trajectories showing proportion change over time
     - Color-coded by treatment arm (red = MK-3475 + MLA, blue = MK-3475 Alone)
     - Each line connects one patient's proportions at two timepoints
   - **Boxplots**: Distribution of proportion changes by treatment arm
     - Jittered points show individual patient values
     - P-value annotation indicates statistical significance
   - Save combined grid plots (line + box) as PDF files

4. **Export Supplemental Table S3:**
   - Columns: celltype, timepoint_A, timepoint_B, comparison, p_value, significance, fold_change
   - CSV file: `comparison_results.csv`
   - Contains statistical results for all cell type × timepoint pair combinations

**SECTION 2: T Cell Subpopulations Analysis (Lines 311-659) → Supplemental Table S4**

**T Cell Subtypes Analyzed (13 functional subpopulations):**

**CD4+ T Helper Cells:**
1. **Activated_CD4** (Cluster 0): Recently activated, early effector functions
2. **Active_CD4** (Cluster 5): Ongoing activation, cytokine production
3. **Naive_CD4** (Clusters 6, 9, 18): Antigen-inexperienced, homeostatic maintenance
4. **Memory_CD4** (Cluster 7): Antigen-experienced, rapid recall response

**CD8+ Cytotoxic T Cells:**
5. **Effector_CD8** (Cluster 1): Terminally differentiated, high cytotoxicity
6. **Effector_Memory_CD8** (Cluster 10): Cytotoxic capacity with memory potential
7. **Effector_Memory_Precursor_CD8** (Cluster 2): Transitional effector-to-memory state
8. **Central_Memory_CD8** (Cluster 12): Long-lived, self-renewal capacity
9. **GZMK_Effector_Memory_CD8** (Cluster 13): GZMK+ cytotoxic memory subset
10. **Stem_Like_CD8** (Cluster 8): Self-renewing, multipotent progenitors

**Specialized T Cell Subsets:**
11. **Exhausted_T** (Cluster 3): High checkpoint expression, impaired function
12. **Gamma_Delta_T** (Cluster 4): Innate-like, non-MHC restricted recognition
13. **Proliferating_Effector** (Clusters 14, 16, 17): Actively dividing effector cells

**Workflow:**

1. Load T cell Seurat object (`MK_T_Cells_seurat_obj.RDS`)
2. Apply identical methodology as Section 1:
   - Calculate proportions for each T cell subtype
   - Compare treatment arms using Wilcoxon test
   - Compute log2 fold changes
   - Generate visualizations (line plots + boxplots)
3. **Export Supplemental Table S4:**
   - Same column structure as Table S3
   - CSV file: `comparison_results_T_Cells.csv`
   - Contains results for all T cell subtype × timepoint pair combinations

**Biological Interpretation:**

**Major Cell Types (Table S3):**
- **Monocyte expansion**: Indicates enhanced antigen presentation and inflammation
- **DC proportion changes**: Reflects alterations in T cell priming capacity
- **NK cell dynamics**: Suggests modulation of innate anti-tumor immunity
- **Cohort differences**: Combination therapy effects on innate immune composition

**T Cell Subpopulations (Table S4):**
- **Effector expansion**: Enhanced anti-tumor cytotoxicity
- **Memory formation**: Long-term immune surveillance potential
- **Exhaustion dynamics**: Checkpoint blockade reversing T cell dysfunction
- **Naive/memory balance**: Immune repertoire diversity and sustainability
- **Proliferation**: Clonal expansion of tumor-reactive T cells
- **Treatment arm comparison**: Identifies combination therapy's impact on T cell differentiation

**Statistical Framework:**
- **Wilcoxon rank-sum test**: Non-parametric comparison of independent groups
  - Robust to outliers and non-normal distributions
  - Tests whether proportion changes differ between treatment arms
- **Log2 fold change**: Quantifies magnitude and direction of change
  - Symmetric scale (doubling = +1, halving = -1)
  - Facilitates comparison across different cell types
- **Significance levels**: Publication-standard thresholds with symbol notation
- **Multiple testing consideration**: Individual comparisons reported; users should apply FDR correction if analyzing all results together

**Change Calculation Methods:**
- **Ratio method** (default for Section 2): Fold-change, interpretable as "X times higher/lower"
- **Difference method** (default for Section 1): Absolute change, sum to zero across all cell types

**Clinical Relevance:**
- Identifies which immune cell populations respond to treatment
- Reveals differential effects of monotherapy vs combination therapy
- Guides understanding of immunotherapy mechanisms of action
- Informs biomarker development for patient stratification

**Key Outputs:**

**Supplemental Table S3 (Major Cell Types):**
- 5 cell types × 3 timepoint pairs = 15 comparisons
- Columns: celltype, timepoint_A, timepoint_B, comparison, p_value, significance, fold_change
- Interpretation: Significant p-values indicate differential treatment effects on innate immunity

**Supplemental Table S4 (T Cell Subpopulations):**
- 13 T cell subtypes × 3 timepoint pairs = 39 comparisons
- Same column structure as Table S3
- Interpretation: Reveals adaptive immune response dynamics and checkpoint blockade effects

**Visualization Outputs:**
- PDF grid plots for each cell type/subtype showing:
  - Top row: Patient trajectory line plots (temporal dynamics)
  - Bottom row: Distribution boxplots (cohort comparisons)
- Organized by method (ratio/difference) subdirectories

**Dependencies:**
- All-cells Seurat object: `MK_Cells_seurat_obj.RDS` (for Section 1)
- T cell Seurat object: `MK_T_Cells_seurat_obj.RDS` (for Section 2)
- Survival data: `patient_survival_data.csv` (contains treatment arm assignments)
- Libraries: Seurat, dplyr, tidyr, ggplot2, cowplot

---

### track_t_cell_subpopulation_clones.R

**Purpose:** Track T cell receptor (TCR) clones across timepoints to visualize spatial distribution and temporal dynamics of antigen-specific T cell responses.

**Manuscript Figure Generated:**
- **Figure 7a**: CD8+ T cell clone distribution across 7 timepoints with density contours colored by patient survival cohort

**Core Concept - TCR Clone Tracking:**
T cell clones are groups of cells sharing identical TCR CDR3β amino acid sequences, indicating they arose from the same antigen-specific precursor cell. Tracking these clones across timepoints reveals:
- **Clonal expansion**: Increase in clone abundance (antigen-specific proliferation)
- **Clonal contraction**: Decrease in clone abundance (death or differentiation)
- **Spatial migration**: Changes in UMAP position (transcriptional state transitions)
- **Temporal persistence**: Clone presence across multiple timepoints (memory formation)

**Patient Cohort Classification (Figure 7a Color Coding):**
- **Yellow (Control)**: NLS+PEM = Non-Long Survivors + PEM arm (patients 1, 4, 8, 9)
  - Control group receiving checkpoint blockade without experimental combination
- **Blue (Short-term)**: LTT+PEM (< mOS) = Long-Term Treatment + PEM, below median overall survival
  - Experimental arm patients with shorter survival (patients 7, 10, 12, 14, 18)
- **Red (Long-term)**: LTT+PEM (> mOS) = Long-Term Treatment + PEM, above median overall survival
  - Experimental arm patients with longer survival (patients 2, 3, 5, 13, 19, 20, 21)

**Clone Tracking Methodology (5-Step Process):**

**Step 1: Subpopulation Cell Identification**
- Extract all cells belonging to specific T cell subpopulation (e.g., Effector_CD8)
- Uses Seurat cluster assignments from annotation

**Step 2: Barcode-to-Clone Mapping**
- Map cell barcodes to their TCR CDR3β sequences
- Source: VDJ TCR sequencing data (`filtered_contig_annotations.csv`)
- Focus on TRB chain (β chain, more diverse than α chain)

**Step 3: Subpopulation-Specific Clone Identification**
- Extract unique CDR3β sequences from subpopulation cells
- These are clones enriched in this specific T cell subset

**Step 4: Clone Repertoire Filtering**
- Cross-reference with quantified clone repertoire data
- Filter for clones with detectable proportion in patient samples
- Ensures clones have sufficient abundance for reliable tracking

**Step 5: Complete Clone Cell Retrieval**
- Identify ALL cells bearing the filtered clones (not just those in focal subpopulation)
- Enables tracking clone migration across different T cell states
- Maps complete temporal trajectory of each clone

**Workflow:**

**Data Preparation (Lines 1-56):**
1. Load T cell Seurat object with UMAP embeddings
2. Filter for UF site patients and remove hyperactivated cluster 13
3. Define 14 T cell subpopulations with cluster mappings:
   - CD4+ subsets: Activated, Active, Naive, Memory
   - CD8+ subsets: Effector, Multiple Memory States, Stem-Like, Anergic, Naive
   - Specialized: Exhausted, Gamma-Delta, Proliferating, General CD8
4. Classify patients into survival cohorts:
   - Control: Patients 1, 4, 8, 9
   - Short-term survivors: Patients 7, 10, 12, 14, 18
   - Long-term survivors: Patients 2, 3, 5, 13, 19, 20, 21
5. Add `SurvivorGroup` metadata to Seurat object

**Function 1: Contour Density Tracking (Lines 57-262) - Figure 7a Style**

6. **Clone Identification Pipeline:**
   - Extract subpopulation cells from Seurat object
   - Map cell barcodes to TCR CDR3β clones
   - Filter for subpopulation-enriched clones
   - Intersect with quantified clone repertoire
   - Retrieve all cells bearing these clones across all timepoints

7. **UMAP Preparation:**
   - Extract UMAP coordinates for all cells
   - Set consistent x/y axis limits across all timepoints (enables comparison)
   - Maintain fixed aspect ratio for accurate spatial representation

8. **Visualization Generation (Per Timepoint):**
   - **Background**: Gray points showing all T cells (context)
   - **Highlighted Clones**: Color-coded by survival cohort
     - Yellow = Control cohort clones
     - Blue = Short-term survivor clones
     - Red = Long-term survivor clones
   - **Density Contours**: 2D kernel density estimation showing clone spatial distribution
   - **Coordinate Flip**: Match manuscript figure orientation

9. **Grid Assembly:**
   - Arrange timepoint plots in chronological order (C1, C2, C4, C6, C9, C17, C34)
   - Save as single PDF with all timepoints side-by-side
   - Enables visual tracking of clone dynamics over treatment course

**Function 2: Filled Density Tracking (Lines 264-488) - Enhanced Visualization**

10. Uses identical clone tracking methodology as Function 1
11. **Enhanced visualization** with filled density polygons instead of contour lines:
    - `stat_density_2d()` with filled polygons (bins=10, alpha=0.5)
    - More prominent visual representation of clone distribution
    - Better for presentations and highlighting spatial patterns
12. Same grid assembly and output as Function 1

**Main Execution Loop (Lines 492-515):**

13. Load input data files:
    - TCR clone-barcode mapping: `filtered_contig_annotations.csv`
    - Clone proportion table: `clonotype_df_proportion_all_T_cells.csv`

14. Iterate over all 14 T cell subpopulations:
    - Create output directory for each subpopulation
    - Execute clone tracking and visualization
    - Error handling for subpopulations with insufficient data

15. Generate two visualization types per subpopulation:
    - **Contour plots**: `{subpop}_clone_tracking_umap_all_patients_all_cells_fixed_axes.pdf` (Figure 7a style)
    - **Filled density plots**: `{subpop}_filled_density_plots_with_background_all_patients_all_cells_fixed_axes.pdf`

**Biological Interpretation:**

**Figure 7a**: CD8+ T Cell Clone Temporal Dynamics
- **Early Timepoints (C1-C2)**: Initial treatment response
  - Red contours (long-term survivors) may show different spatial patterns than blue (short-term)
  - Yellow (control) provides baseline reference
- **Middle Timepoints (C4-C6)**: Peak response phase
  - Clone expansion visible as larger, denser contours
  - Spatial repositioning indicates transcriptional state changes
- **Late Timepoints (C9-C34)**: Long-term immunity
  - Persistent contours indicate memory formation
  - Contraction or expansion reveals durability of response

**Clone Spatial Patterns:**
- **Localized contours**: Clones in specific transcriptional state
- **Dispersed contours**: Clones spanning multiple states (heterogeneity)
- **Overlapping cohorts**: Similar spatial distribution suggests shared biology
- **Separated cohorts**: Distinct patterns may explain survival differences

**Clinical Insights:**
- **Clone persistence**: Long-term survivors show sustained responses
- **Clone expansion dynamics**: Different expansion kinetics between cohorts
- **Subpopulation specificity**: Certain T cell subsets may drive therapeutic benefit
- **Spatial trajectories**: Clone migration patterns reveal differentiation pathways

**Key Outputs:**

**For Each T Cell Subpopulation (14 total):**
- Contour density plots (Figure 7a manuscript style)
- Filled density plots (enhanced visualization)
- Both showing 7-8 timepoints in chronological sequence

**Visualization Design:**
- Gray background: All T cells (provides anatomical context)
- Colored density contours/fills: Tracked clones (shows clone dynamics)
- Fixed axes across timepoints: Enables direct comparison
- Coordinate flip: Matches manuscript orientation

**Dependencies:**
- T cell Seurat object: `MK_T_Cells_seurat_obj.RDS`
- TCR data: `filtered_contig_annotations.csv` (VDJ sequencing, CDR3β sequences)
- Clone proportions: `clonotype_df_proportion_all_T_cells.csv`
- Libraries: Seurat, ggplot2, cowplot, dplyr

**Technical Notes:**
- Requires ≥3 points for density estimation (contingency for sparse clones)
- Uses 2D kernel density estimation for smooth contours
- Exports both contour and filled versions for flexibility
- Error handling for subpopulations with insufficient clone data

---

## Scripts Pending Documentation

- clonal_expansion.R
- get_metadata.R
- immune_checkpoint_comparison.R
- neoantigen_analysis.R
- optimal_transport.R
- short_term_long_term_analysis.R

---

*Last Updated: 2024-12-12*