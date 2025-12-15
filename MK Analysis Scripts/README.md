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

### ✅ get_metadata.R

**Purpose:** Utility script providing standardized timepoint nomenclature mapping and metadata sorting functions for clonal expansion analyses.

**Script Type:** Supporting utility (no figures/tables generated)

**Core Functionality:**

**1. Timepoint Standardization (`sample_order` vector):**
- Defines canonical timepoint ordering: S/Pre, C1, C2, C4, C6, C9, C18, C36
- Ensures consistent temporal ordering across all analyses
- Maps multiple nomenclature systems:
  - "S" = Screening/baseline (legacy nomenclature)
  - "Pre" = Pre-treatment baseline (current nomenclature)
  - "C1"-"C36" = Treatment cycles 1-36

**2. Metadata Sorting (`get_metadata_sorted()` function):**
- **Purpose**: Sort and standardize timepoint metadata with ordered factor levels
- **Process**:
  1. Extracts timepoint from sample origin by removing patient ID prefix
  2. Converts timepoint strings to ordered factor using `sample_order` levels
  3. Sorts metadata by timepoint chronological order
  4. Ensures consistent temporal ordering for downstream analyses
- **Use Cases**:
  - Preparing metadata for clonal expansion calculations
  - Ensuring chronological ordering in plots and tables
  - Standardizing timepoint labels across different data sources

**3. Boolean Conversion (`to_bool()` function):**
- **Purpose**: Convert string representations to boolean values
- **Mappings**:
  - "TRUE" or "True" or "true" → TRUE
  - "FALSE" or "False" or "false" → FALSE
  - All other inputs → FALSE (default)
- **Use Cases**: Parsing configuration files or user input

**Biological Context:**
This script supports clonal expansion analysis by standardizing timepoint nomenclature. Different data sources may use "S" vs "Pre" for baseline, or have varying timepoint label formats. Consistent ordering is critical for:
- Calculating temporal changes in clonal diversity (C1 vs Pre, C2 vs C1, etc.)
- Plotting time-series data with correct chronological sequence
- Merging datasets from different experimental batches

**Dependencies:**
- No external data files required
- Used by: clonal expansion pipeline scripts (Steps 1-4), other temporal analysis scripts

**Key Output:**
- Sorted metadata dataframes with standardized timepoint factor levels
- Enables accurate temporal comparisons and visualizations

---

### ✅ Clonal_Tracking/final_nomenclature/coi/optimal_transport.R

**Purpose:** Data preparation script for optimal transport analysis of T cell clones across timepoints - extracts source and target cell data with UMAP/PCA coordinates for Python-based optimal transport calculations.

**Script Type:** Data preparation utility (prepares inputs for `t_cell_optimal_transport.ipynb`)

**Core Concept - Optimal Transport Data Preparation:**
This R script identifies subpopulation-enriched TCR clones and organizes cell data into structured CSV files that enable optimal transport analysis in Python. It bridges the gap between Seurat single-cell analysis and mathematical optimal transport modeling.

**Clone Identification Pipeline (5 Steps):**

**Step 1: Subpopulation Cell Identification**
- Extract cells belonging to specific CD8+ T cell subpopulation (e.g., Effector_CD8)
- Uses Seurat cluster assignments from annotation

**Step 2: Barcode-to-Clone Mapping**
- Maps cell barcodes to TCR CDR3β amino acid sequences (clone identifiers)
- Source: VDJ TCR sequencing data
- Focuses on TRB chain (β chain, more diverse than α chain)

**Step 3: Subpopulation-Enriched Clone Identification**
- Extracts unique CDR3β sequences from subpopulation cells
- Identifies clones enriched in this specific T cell subset

**Step 4: Clone Repertoire Filtering**
- Cross-references with quantified clone proportion data
- Filters for clones with detectable abundance (proportion > 0) in patient samples
- Ensures only high-confidence clones are tracked

**Step 5: Complete Clone Cell Retrieval**
- Retrieves ALL cells bearing the filtered clones across ALL timepoints
- Not limited to focal subpopulation - enables tracking clone migration across T cell states
- Example: Effector_CD8 clones can be tracked as they transition to Memory states

**Workflow:**

**Data Loading and Cohort Definition:**
1. Load T cell Seurat object (UF site patients only)
2. Filter for CD8+ T cell clusters (1, 2, 3, 8, 10, 12, 14, 16, 17)
3. Remove patient-specific cluster 13 (hyperactivated/artifact)
4. Define patient survival cohorts:
   - Control: Patients 1, 4, 8, 9 (NLS+PEM)
   - Short-term: Patients 7, 10, 12, 14, 18 (LTT+PEM, below median OS)
   - Long-term: Patients 2, 3, 5, 13, 19, 20, 21 (LTT+PEM, above median OS)

**Main Function: `extractDataForOptimalTransport()`:**

5. Subset Seurat object by patient cohort and T cell subpopulation
6. Load TCR VDJ data and clone proportion tables
7. Extract UMAP coordinates (for visualization)
8. Extract PCA coordinates (for optimal transport distance calculations)
9. Execute 5-step clone identification pipeline
10. For each consecutive timepoint pair (Pre→C1, C1→C2, etc.):
    - **Source cells**: Highlighted clones at timepoint A (with UMAP/PCA)
    - **Target cells**: Highlighted clones at timepoint B (with UMAP/PCA)
    - **Gray cells**: Background cells at timepoint A (provides spatial context)
11. Save three CSV files per timepoint pair:
    - `source_cells_new.csv`: Source timepoint data
    - `target_cells_new.csv`: Target timepoint data
    - `gray_cells_new.csv`: Background cells for context
12. Special handling for last timepoint (no target available)

**Nested Loop Execution:**
13. Outer loop: 8 CD8+ T cell subpopulations
14. Inner loop: 3 patient survival cohorts
15. Total: 8 × 3 = 24 analysis units
16. Each unit generates CSV files for all available timepoint pairs

**CD8+ T Cell Subpopulations Analyzed:**
- **Effector_CD8** (Cluster 1): Terminally differentiated cytotoxic
- **Memory_Precursor_Effector_CD8** (Cluster 2): Transitional state
- **Exhausted_T** (Cluster 3): Dysfunctional, high checkpoint expression
- **Stem_Like_CD8** (Cluster 8): Self-renewing progenitors
- **Effector_Memory_CD8** (Cluster 10): Cytotoxic with memory
- **Central_Memory_CD8** (Cluster 12): Long-lived self-renewal
- **Proliferating_Effector** (Clusters 14, 16, 17): Actively dividing
- **All_CD8** (Combined): All CD8+ populations

**Output Structure:**
```
base_output_dir/
├── Effector_CD8/
│   ├── control/
│   │   ├── Timepoint_Pre/
│   │   │   ├── source_cells_new.csv
│   │   │   ├── target_cells_new.csv
│   │   │   └── gray_cells_new.csv
│   │   ├── Timepoint_C1/
│   │   └── ... (subsequent timepoints)
│   ├── short_term/
│   └── long_term/
└── ... (other subpopulations)
```

**CSV File Contents:**
Each CSV contains:
- **cell_id**: Cell barcode identifier
- **TimePoint**: Timepoint label (Pre, C1, C2, etc.)
- **SurvivorGroup**: Patient cohort (control, short_term, long_term)
- **seurat_clusters**: Cluster assignment
- **origin**: Sample ID
- **Patient**: Patient identifier
- **UMAP_1, UMAP_2**: 2D UMAP coordinates
- **PC_1 through PC_N**: All PCA dimensions

**Downstream Usage:**
The CSV files are loaded by `t_cell_optimal_transport.ipynb` which:
1. Computes optimal transport plan using Earth Mover's Distance
2. Generates arrow plots showing cell fate transitions
3. Creates density shift visualizations
4. Exports predecessor mapping for clone tracking analysis

**Biological Rationale:**
- **Optimal transport** predicts which source cells transform into which target cells
- **Clone tracking** enables following antigen-specific T cells through treatment
- **Subpopulation focus** reveals differentiation pathways (e.g., Effector → Memory)
- **Cohort comparison** identifies survival-associated clone dynamics

**Key Features:**
- Extracts both UMAP (visualization) and PCA (distance calculation) coordinates
- Handles sparse timepoints gracefully (not all patients have all timepoints)
- Maintains complete clone tracking (retrieves all clone-bearing cells, not just subpopulation)
- Organized output structure for systematic Python analysis

**Dependencies:**
- T cell Seurat object: `MK_T_Cells_seurat_obj.RDS`
- TCR VDJ data: `filtered_contig_annotations.csv` (barcode-to-CDR3β mapping)
- Clone proportions: `clonotype_df_proportion_all_T_cells.csv`
- Libraries: Seurat, dplyr, ggplot2

**Next Step:** Run `t_cell_optimal_transport.ipynb` to perform optimal transport analysis and generate visualizations.

---

### ✅ Clonal_Tracking/final_nomenclature/coi/t_cell_optimal_transport.ipynb

**Purpose:** Perform optimal transport analysis on T cell clones to predict cell fate transitions and differentiation patterns across treatment timepoints, generating movement visualizations and cluster distribution analyses.

**Manuscript Figures Generated:**
- **Figure 7b**: Cluster-level arrow plots showing T cell subpopulation movement patterns across timepoints
- **Figure 8a**: Stacked bar charts showing target cluster distributions for each source cluster of interest

**Core Concept - Optimal Transport for Cell Fate Prediction:**
Optimal transport provides a mathematical framework to predict which source cells (at timepoint A) transform into which target cells (at timepoint B) by solving the "Earth Mover's Distance" problem. This enables:
- **Cell fate prediction**: Identify most likely developmental trajectories
- **Spatial migration analysis**: Visualize transcriptional state transitions in UMAP space
- **Differentiation pattern quantification**: Measure cluster-to-cluster transformation probabilities
- **Cohort comparison**: Compare clone dynamics between survival groups

**Mathematical Framework:**

**1. Earth Mover's Distance (EMD):**
- Finds minimum-cost transformation from source distribution to target distribution
- Also known as Wasserstein-1 distance
- Formula: EMD = min Σ T[i,j] × C[i,j], where T is transport plan and C is cost matrix

**2. Cost Matrix:**
- Pairwise Euclidean distances between all source and target cells
- **Computed in PCA space** (high-dimensional, captures full transcriptional differences)
- Dimensions: N_source × N_target

**3. Transport Plan:**
- Matrix T where T[i,j] = probability mass moved from source cell i to target cell j
- Uniform distributions over source and target cells (equal weighting)
- Solver: `ot.emd()` (exact EMD solver from Python Optimal Transport library)

**4. Visualization Mapping:**
- Transport plan computed in PCA space (accurate distances)
- Arrows visualized in UMAP space (interpretable 2D representation)
- Displacement vectors: target_UMAP - source_UMAP for each cell

**Workflow:**

**PART 1: Data Loading and Preprocessing**

1. **Configuration Setup:**
   - Define 8 CD8+ T cell subpopulations to analyze
   - Define 3 patient cohorts (control, short_term, long_term)
   - Define 8 timepoints (Pre, C1, C2, C4, C6, C9, C18, C36)
   - Set cohort color mapping: yellow=control, blue=short_term, red=long_term
   - Define clusters of interest (COI): [1, 2, 3, 8, 10, 12, 14]

2. **Load Cluster Color Mapping:**
   - Read T_Cell_cluster_colors.csv for consistent cluster visualization
   - Maps cluster IDs to colors and cell type names

**PART 2: Main Function - `optimal_transport_visualization()`**

**Data Loading Phase:**

3. Load CSV files for specified subpopulation and cohort:
   - `source_cells_new.csv`: Highlighted clones at each timepoint
   - `gray_cells_new.csv`: Background cells for context
4. Identify available timepoint folders and extract timepoint labels
5. Gather data across all valid timepoints

**Axis Limits and Downsampling:**

6. Calculate global UMAP axis limits (x_min, x_max, y_min, y_max)
   - Ensures consistent visualization across all timepoints
7. Find minimum cell counts across timepoints for balanced downsampling
   - min_source_cells: Smallest number of highlighted clone cells
   - min_gray_cells: Smallest number of background cells
8. Downsample all timepoints to these minimum values (random_state=42)
   - Prevents timepoint bias in visualizations
   - Ensures equal representation across treatment course

**Optimal Transport Computation:**

9. For each consecutive timepoint pair (Pre→C1, C1→C2, etc.):
   
   **Step 9a: Extract PCA/UMAP Coordinates**
   - Identify PCA columns (PC_1, PC_2, ..., PC_N)
   - Extract full PCA coordinates for source and target cells
   - Extract UMAP coordinates for visualization
   
   **Step 9b: Define Uniform Mass Distributions**
   - Create uniform distribution over source cells: a = [1/N_source, ..., 1/N_source]
   - Create uniform distribution over target cells: b = [1/N_target, ..., 1/N_target]
   - Assumes equal weighting for all cells
   
   **Step 9c: Compute Cost Matrix**
   - Calculate pairwise Euclidean distances in PCA space
   - Cost matrix C: N_source × N_target
   - Higher dimensions capture full transcriptional landscape
   
   **Step 9d: Solve Optimal Transport**
   - Call `ot.emd(a, b, cost_matrix)` to compute transport plan T
   - T[i,j] = mass moved from source cell i to target cell j
   - Maximum iterations: 100,000 (complex optimization problem)
   
   **Step 9e: Extract Cell Fate Predictions**
   - For each source cell, find target cell with maximum transport probability
   - `target_indices = argmax(transport_plan, axis=1)`
   - Compute displacement vectors in UMAP space for visualization
   - Calculate Euclidean norms (arrow lengths)

**Single-Cell Arrow Analysis:**

10. Extract displacement vectors for downsampled source cells
11. Normalize displacement vectors to unit vectors
12. Scale by original arrow lengths (differential arrow lengths by distance)

**Cluster-Level Analysis:**

13. **Cluster Arrow Computation:**
    - Group source cells by Seurat cluster assignment
    - Lump related clusters: {16→14, 17→14}, {9→6, 18→6}
    - Calculate cluster centroids (median UMAP position)
    - Compute mean displacement vector per cluster
    - Weight arrow thickness by cluster proportion at source timepoint
    - Filter to clusters of interest (COI) only

14. **Aggregated Arrow Computation:**
    - Sum displacement vectors across all cluster arrows
    - Start point: median UMAP position of all source cells
    - Direction/magnitude: Vector sum of cluster-level displacements
    - Represents overall population movement

**Target Cluster Distribution Analysis (Figure 8a):**

15. **For each source cluster of interest (COI):**
    - Identify which target clusters COI cells map to via transport plan
    - Calculate fraction mapping to each target cluster
    - Formula: fraction[target_cluster] = count(COI→target) / total(COI)
    - Store distributions for stacked bar chart generation

**PART 3: Visualization Generation**

**Three PDF Types Generated Per Subpopulation-Cohort:**

**Visualization 1: Single-Cell Movement Plots**

16. Create multi-panel plot (one panel per timepoint)
17. For each timepoint panel:
    - Plot gray background cells (alpha=0.5, size=5)
    - Plot highlighted clone cells (cohort color, alpha=1, size=5)
    - Draw black arrows showing individual cell fate predictions
    - Arrow properties:
      - Direction: From source to predicted target location
      - Length: Proportional to displacement magnitude
      - Styling: head_width=0.05, linewidth=0.5, alpha=0.7
18. Save as `{subpop}_{cohort}_movement_plots_differential_arrow_lengths_equal_cells_using_PCs.pdf`

**Visualization 2: Cluster-Level Arrow Plots (Figure 7b)**

19. Create multi-panel plot (one panel per timepoint)
20. For each timepoint panel:
    - Plot gray background + highlighted clones (same as Viz 1)
    - Draw cluster-level arrows:
      - Start: Cluster centroid (median position)
      - Direction: Mean displacement of cluster cells
      - Thickness: Proportional to cluster proportion (linewidth = 1.0 + 8.0 × proportion)
      - Arrow size: head_width and head_length scale with proportion
      - Color: Black with black outline (cluster-specific colors commented out)
    - Only show arrows for clusters of interest
21. Save as `{subpop}_{cohort}_movement_plots_aggregated_cluster_arrows_equal_cells_using_PCs.pdf`
22. **This generates Figure 7b components**

**Visualization 3: Single Aggregated Arrow**

23. Create multi-panel plot (one panel per timepoint)
24. For each timepoint panel:
    - Plot gray background + highlighted clones
    - Draw single aggregated arrow:
      - Start: Global median of all source cells
      - Direction/magnitude: Sum of all cluster-level displacement vectors
      - Styling: Black, thick (linewidth=2.0), prominent (head_width=0.3)
25. Save as `{subpop}_{cohort}_movement_plots_aggregated_single_arrow_equal_cells_using_PCs.pdf`

**PART 4: Cluster Distribution Stacked Bar Charts (Figure 8a)**

**Target Cluster Grouping:**

26. For each source cluster of interest (COI):
    - Compile target cluster distributions across all timepoints and cohorts
    - Group non-COI target clusters into "Others" category
    - This simplifies visualization by focusing on COI-to-COI transitions

**Stacked Bar Chart Generation:**

27. Create large grid plot:
    - **Rows**: One row per source cluster of interest (7 clusters)
    - **Columns**: One column per timepoint (8 timepoints)
    - **Total**: 7 × 8 = 56 subplots

28. For each subplot (source cluster × timepoint):
    - X-axis: Three bars representing three cohorts (control, short_term, long_term)
    - Y-axis: Fraction mapping to each target cluster (stacked to 1.0)
    - **Stacking order**: Largest average fraction at bottom (most visible)
    - **Color coding**: Each target cluster has consistent color (from color mapping)
    - **Text labels**: Cell type names shown for fractions >5%

29. Add comprehensive labeling:
    - Column titles: Timepoint labels (Pre, C1, C2, etc.)
    - Row labels: Source cluster cell type names (on left margin)
    - Legend: Target cluster color-to-celltype mapping (right side)
    - Overall title: Indicates subpopulation being analyzed

30. Save as `{subpop}_target_cluster_distribution_coi_with_stack_labels.pdf`
31. **These stacked bar charts are Figure 8a**

**PART 5: Main Execution Loop**

32. **Nested loop structure:**
    - Outer loop: Iterate over 8 CD8+ T cell subpopulations
    - Inner loop: Iterate over 3 patient cohorts
    - Total: 8 × 3 = 24 optimal transport analyses

33. For each subpopulation-cohort combination:
    - Run `optimal_transport_visualization()` function
    - Collect target cluster distributions
    - Generate all three movement visualization types

34. **Aggregate distributions across cohorts:**
    - Combine distribution data from all three cohorts
    - Enable comparison in single stacked bar chart

35. Generate final stacked bar chart PDF for this subpopulation
36. Repeat for all 8 subpopulations

**Biological Interpretation:**

**Figure 7b - Cluster-Level Movement:**
- **Arrow direction**: Reveals dominant differentiation pathways
  - Rightward: Often represents maturation/differentiation
  - Leftward: May indicate dedifferentiation or state reversal
- **Arrow length**: Quantifies magnitude of transcriptional change
  - Longer arrows = Greater state transition
  - Shorter arrows = Stable transcriptional state
- **Arrow thickness**: Indicates cluster abundance
  - Thicker arrows = More cells in this cluster
  - Thinner arrows = Rare subpopulation
- **Cohort comparison**: Different colored plots reveal survival-associated patterns
  - Control vs short-term vs long-term differential dynamics

**Figure 8a - Target Cluster Distribution:**
- **Stacked bars show differentiation fate:**
  - Which target clusters do source COI cells become?
  - Example: Do Effector_CD8 (cluster 1) cells → Memory states or Exhausted states?
- **Temporal evolution:**
  - Columns show progression through treatment (Pre → C36)
  - Track how differentiation patterns change over time
- **Cohort stratification:**
  - Three bars per timepoint compare cohorts
  - Identify survival-associated differentiation programs
  - Example: Long-term survivors may show more Memory formation
- **"Others" category:**
  - Cells transitioning to non-COI clusters
  - Indicates differentiation out of CD8+ effector programs
  - Could represent cell death, state change, or technical dropout

**Key Differentiation Patterns:**
- **Effector → Memory**: Functional maturation and long-term immunity
- **Effector → Exhausted**: T cell dysfunction accumulation
- **Memory → Effector**: Reactivation and re-engagement
- **Stable clusters**: Cells maintaining transcriptional state
- **→ Others**: Differentiation to non-CD8 states or technical loss

**Clinical Insights:**
- **Long-term survivors**: May show preferential Memory formation
- **Short-term survivors**: Could exhibit Exhaustion accumulation
- **Control cohort**: Baseline differentiation dynamics without combination therapy
- **Treatment effects**: Changes in differentiation patterns induced by MK-3475 + MLA

**Technical Implementation:**

**Optimal Transport Parameters:**
- **Distance metric**: Euclidean distance in PCA space (all principal components)
- **Transport solver**: `ot.emd()` (exact Earth Mover's Distance)
- **Maximum iterations**: 100,000 (ensures convergence)
- **Mass distributions**: Uniform (equal weighting for all cells)

**Visualization Parameters:**
- **Downsampling**: Equal cell numbers across timepoints (prevents bias)
- **Fixed axes**: Consistent UMAP limits across all panels
- **Aspect ratio**: Equal (prevents distortion)
- **Color coding**: Cohort-specific (yellow/blue/red)
- **Arrow styling**: Differential lengths, proportion-weighted thickness

**Key Outputs:**

**For Each of 8 Subpopulations × 3 Cohorts = 24 Analysis Units:**

**Movement Visualizations:**
1. `{subpop}_{cohort}_movement_plots_differential_arrow_lengths_equal_cells_using_PCs.pdf`
   - Single-cell level arrows (individual fate predictions)
   - Multi-panel: One panel per timepoint
   - Black arrows showing predicted transitions

2. `{subpop}_{cohort}_movement_plots_aggregated_cluster_arrows_equal_cells_using_PCs.pdf`
   - Cluster-level arrows (average movement per subpopulation)
   - Multi-panel: One panel per timepoint
   - Arrow thickness proportional to cluster abundance
   - **Components of Figure 7b**

3. `{subpop}_{cohort}_movement_plots_aggregated_single_arrow_equal_cells_using_PCs.pdf`
   - Single aggregated arrow (overall population shift)
   - Multi-panel: One panel per timepoint
   - Large black arrow showing net movement

**Cluster Distribution Stacked Bar Charts (Figure 8a):**

4. `{subpop}_target_cluster_distribution_coi_with_stack_labels.pdf`
   - Grid: 7 source COI clusters × 8 timepoints = 56 subplots
   - Each subplot: 3 stacked bars (one per cohort)
   - Bars show fraction mapping to each target cluster
   - **Manuscript Figure 8a**

**Data Flow:**

```
optimal_transport.R (R)
    ↓  Prepares CSV files
t_cell_optimal_transport.ipynb (Python)
    ↓  Performs OT analysis
Movement visualizations + Stacked bar charts
    ↓  Manuscript figures
Figure 7b (cluster arrows) + Figure 8a (distributions)
```

**Clusters of Interest (COI) - Source Clusters:**
- **Cluster 1**: Effector_CD8 (terminally differentiated)
- **Cluster 2**: Memory_Precursor_Effector_CD8 (transitional)
- **Cluster 3**: Exhausted_T (dysfunctional)
- **Cluster 8**: Stem_Like_CD8 (self-renewing)
- **Cluster 10**: Effector_Memory_CD8 (cytotoxic memory)
- **Cluster 12**: Central_Memory_CD8 (long-lived)
- **Cluster 14**: Proliferating_Effector (actively dividing, includes 16, 17)

**Target Clusters:**
- All COI clusters (track COI-to-COI transitions)
- "Others" category (non-COI clusters, representing differentiation out of CD8+ effector programs)

**Cluster Lumping Rules:**
- Clusters 16, 17 → 14 (all Proliferating_Effector)
- Clusters 9, 18 → 6 (all Naive_CD4)
- Simplifies visualization while maintaining biological meaning

**Patient Cohorts:**
- **Control** (Yellow): Patients 1, 4, 8, 9 (NLS+PEM arm)
- **Short-term** (Blue): Patients 7, 10, 12, 14, 18 (LTT+PEM, below median OS)
- **Long-term** (Red): Patients 2, 3, 5, 13, 19, 20, 21 (LTT+PEM, above median OS)

**Statistical Considerations:**
- **Downsampling strategy**: Removes timepoint-dependent sampling bias
- **Uniform mass**: Treats all cells equally (no a priori weighting)
- **PCA-based distances**: Captures full transcriptional differences (not just UMAP projection)
- **Deterministic transport**: Exact EMD solver (not stochastic)
- **Fixed axes**: Enables direct visual comparison across timepoints

**Biological Interpretation:**

**Movement Patterns (Figure 7b):**
- **Directional consistency**: Similar arrows across cohorts suggest conserved biology
- **Directional divergence**: Different arrows indicate cohort-specific differentiation
- **Arrow magnitude**: Quantifies transcriptional state change
- **Temporal progression**: Sequential panels show treatment-induced evolution

**Distribution Patterns (Figure 8a):**
- **High fraction to Memory clusters**: Functional maturation, protective immunity
- **High fraction to Exhausted cluster**: T cell dysfunction, poor prognosis
- **High fraction to "Others"**: Loss of CD8+ effector phenotype
- **Cohort differences**: Differential stacking reveals survival-associated patterns
- **Temporal trends**: Columns show how fates evolve through treatment

**Clinical Relevance:**
- Identifies differentiation programs associated with therapeutic benefit
- Reveals how combination therapy (MK-3475 + MLA) alters T cell fate decisions
- Quantifies memory formation vs exhaustion accumulation
- Guides biomarker development and patient stratification

**Dependencies:**
- **Input CSV files**: Generated by `optimal_transport.R`
- **Python libraries**: pandas, numpy, matplotlib, ot (Python Optimal Transport), PdfPages
- **Cluster mapping**: `T_Cell_cluster_colors.csv`
- **Prior steps**: optimal_transport.R must be run first to generate CSV inputs

**Technical Notes:**
- Handles missing timepoints gracefully (skips if data absent)
- Error handling for failed OT convergence
- Fixed random seed (random_state=42) for reproducibility
- Memory-efficient: Processes one subpopulation-cohort at a time
- Vectorized NumPy operations for computational efficiency

---

### ✅ Clonal_Tracking/final_nomenclature/coi_central_memory_cd8/central_memory_cd8_t_cell_optimal_transport.ipynb

**Purpose:** Focused optimal transport analysis specifically for Central Memory CD8 T cells (cluster 12) to reveal differentiation patterns and fate transitions of this critical long-lived memory population.

**Manuscript Figure Generated:**
- **Figure 8b**: Stacked bar chart showing Central Memory CD8 target cluster distributions across timepoints and cohorts

**Relationship to Main Analysis:**
This notebook is a **focused version** of `t_cell_optimal_transport.ipynb`, analyzing only Central Memory CD8 cells instead of all 7 CD8+ subpopulations. The analysis pipeline is identical, but the cluster of interest (COI) is restricted to cluster 12 (Central Memory CD8).

**Key Difference from Main OT Notebook:**
- **Main notebook** (`t_cell_optimal_transport.ipynb`): Analyzes 7 source COI clusters → Figure 8a (7-row stacked bar chart)
- **This notebook**: Analyzes 1 source COI cluster (cluster 12 only) → Figure 8b (1-row stacked bar chart)

**Cluster of Interest:**
- **Cluster 12 - Central Memory CD8 T cells**:
  - Long-lived memory T cells with self-renewal capacity
  - Critical for sustained anti-tumor immunity
  - Express CD62L (lymph node homing) and CCR7
  - Low effector function but high proliferative potential
  - Can differentiate to Effector Memory or re-activate to Effectors

**Why Focus on Central Memory CD8?**
Central Memory CD8 T cells are particularly important for long-term therapeutic benefit:
- **Memory persistence**: Maintained throughout treatment course
- **Self-renewal**: Can replenish effector populations
- **Survival association**: Abundance correlates with patient outcomes
- **Differentiation plasticity**: Can transition to multiple effector/memory states
- **Clinical biomarker potential**: May predict long-term treatment response

**Mathematical Framework:**
Identical to main OT notebook:
- Optimal transport using Earth Mover's Distance
- PCA-based cost matrix (high-dimensional distances)
- UMAP-based visualization (2D arrows)
- Uniform mass distributions

**Workflow:**
1. Load CSV files for Central Memory CD8 clones (from `optimal_transport.R`)
2. Compute optimal transport plan for consecutive timepoint pairs
3. Generate three movement visualization types:
   - Single-cell arrows (individual fates)
   - Cluster-level arrows (weighted by proportion)
   - Aggregated arrow (overall shift)
4. **Analyze target cluster distributions** (key output for Figure 8b)
5. Generate stacked bar chart showing differentiation fates

**Target Cluster Analysis (Figure 8b Focus):**
For Central Memory CD8 cells at each timepoint, the stacked bar chart shows:
- **Fraction staying Central Memory (cluster 12)**: Self-renewal and maintenance
- **Fraction → Effector Memory (cluster 10)**: Functional maturation
- **Fraction → Effector CD8 (cluster 1)**: Re-activation and anti-tumor response
- **Fraction → Stem-Like CD8 (cluster 8)**: Dedifferentiation to progenitor state
- **Fraction → Others**: Differentiation out of CD8+ memory program

**Stacked Bar Chart Structure (Figure 8b):**
- **Grid**: 1 row (Central Memory only) × 8 timepoints (Pre through C36)
- **X-axis per subplot**: Three bars (control, short_term, long_term cohorts)
- **Y-axis**: Fraction mapping to each target cluster (stacked to 1.0)
- **Color coding**: Each target cluster has consistent color
- **Text labels**: Cell type names for fractions >5%

**Biological Interpretation:**

**Central Memory → Effector Transitions:**
- **High Central Memory self-renewal**: Sustained long-term immunity
- **Increased → Effector Memory**: Functional differentiation while maintaining memory
- **Increased → Effector CD8**: Re-activation and anti-tumor engagement
- **Increased → "Others"**: Loss of memory phenotype (exhaustion or terminal differentiation)

**Cohort Comparisons:**
- **Long-term survivors**: May show more Central Memory self-renewal
- **Short-term survivors**: May exhibit greater terminal differentiation or exhaustion
- **Control cohort**: Baseline differentiation dynamics without combination therapy

**Temporal Evolution:**
- **Early (Pre→C1→C2)**: Initial treatment response and differentiation
- **Middle (C4→C6)**: Peak effector re-activation
- **Late (C9→C36)**: Long-term memory maintenance patterns

**Clinical Insights:**
- Central Memory persistence predicts durable responses
- Differentiation patterns may stratify patients by prognosis
- Combination therapy effects on memory cell fate decisions
- Identifies mechanisms of long-term therapeutic benefit

**Key Outputs:**

**For Each of 8 Subpopulations × 3 Cohorts:**
1. Single-cell movement plots (PDF)
2. Cluster-level arrow plots (PDF)
3. Aggregated arrow plots (PDF)

**Main Figure Output:**
- `{subpop}_target_cluster_distribution_coi_with_stack_labels.pdf`
- **Figure 8b**: Central Memory CD8 target distribution stacked bar chart
- Shows differentiation fates across all timepoints and cohorts

**Configuration:**
- **clusters_of_interest**: [12] (Central Memory CD8 only)
- **Output directory**: `coi_central_memory_cd8/` (separate from main analysis)
- **All other parameters**: Identical to main OT notebook

**Dependencies:**
- **Input CSV files**: From `optimal_transport.R` (same as main analysis)
- **Cluster color mapping**: `T_Cell_cluster_colors.csv`
- **Python libraries**: pandas, numpy, matplotlib, ot, PdfPages

**Note:** This analysis uses the same input data as the main OT notebook but focuses visualization and interpretation specifically on Central Memory CD8 differentiation patterns, which are particularly clinically relevant for long-term therapeutic outcomes.

---

### ✅ immune_checkpoint_comparison.R

**Purpose:** Compare immune checkpoint gene expression across T cell populations between patient survival cohorts and treatment timepoints, visualizing spatial expression patterns and temporal dynamics.

**Manuscript Figure Generated:**
- **Figure 8c**: 2D density plots showing immune checkpoint expression patterns across timepoints and cohorts

**Core Concept - Immune Checkpoint Expression:**
Immune checkpoint molecules are regulatory proteins that modulate T cell activation and function:
- **Inhibitory checkpoints** (e.g., PD-1, CTLA-4, LAG3): Suppress T cell activity, prevent autoimmunity
- **Stimulatory checkpoints** (e.g., CD28, ICOS, 4-1BB): Enhance T cell activation and survival
- **Expression levels**: Correlate with T cell functional state (activated, exhausted, anergic)
- **Clinical relevance**: Therapeutic targets for checkpoint blockade therapy (anti-PD-1 = MK-3475)

**Checkpoint Genes Analyzed (21 total):**

**Inhibitory Checkpoints:**
- **PDCD1** (PD-1): Primary target of MK-3475 therapy
- **CD274** (PD-L1): PD-1 ligand, tumor escape mechanism
- **PDCD1LG2** (PD-L2): Alternative PD-1 ligand
- **CTLA4**: T cell activation inhibitor, second-generation checkpoint target
- **LAG3**: MHC class II binding, co-inhibitory
- **HAVCR2** (TIM-3): T cell exhaustion marker
- **TIGIT**: NK and T cell inhibitory receptor
- **BTLA**: B and T lymphocyte attenuator
- **VSIR** (VISTA): V-domain Ig suppressor of T cell activation
- **IDO1, IDO2**: Tryptophan catabolism, immunosuppressive

**Stimulatory Checkpoints:**
- **CD28**: Costimulatory molecule, essential for T cell activation
- **ICOS**: Inducible costimulator, T cell help and memory
- **TNFRSF4** (OX40): T cell survival and memory formation
- **TNFRSF9** (4-1BB): T cell expansion and persistence
- **TNFRSF18** (GITR): Regulatory T cell modulation

**Ligands:**
- **CD276** (B7-H3): Coinhibitory ligand
- **VTCN1** (B7-H4): Peripheral tolerance
- **NCR3LG1** (B7-H6): Inhibitory ligand
- **HHLA2** (B7-H7): Emerging checkpoint target
- **CD47**: "Don't eat me" signal, innate immunity

**Script Structure:**

**SECTION 1: FeaturePlot Generation with Data Export (Lines 1-251)**

**Workflow:**

1. **Define Patient Cohorts:**
   - Control: Patients 1, 4, 8, 9, 11 (NLS+PEM)
   - Short-term: Patients 7, 10, 12, 14, 18 (LTT+PEM, below median OS)
   - Long-term: Patients 2, 3, 5, 13, 19, 20, 21 (LTT+PEM, above median OS)

2. **Load and Filter Data:**
   - Load T cell Seurat object (all subpopulations)
   - Remove cluster 13 (patient-specific hyperactivation)
   - Filter for UF site patients matching survival data

3. **Add Survival Metadata:**
   - Map patients to survival groups (Control, ShortTerm, LongTerm)
   - Create ordered factor levels for consistent visualization

4. **Timepoint Selection and Downsampling:**
   - Focus on Pre, C1, C2 timepoints (early treatment response)
   - **Critical downsampling strategy**: For each survival group, ensure equal cells across timepoints
   - Formula: n_min = min(n_Pre, n_C1, n_C2) per group
   - Prevents timepoint-dependent sampling bias in expression comparisons

5. **Create Combined Metadata:**
   - Generate `Group_TimePoint` variable (e.g., "Control_Pre", "ShortTerm_C1")
   - Order: Control → ShortTerm → LongTerm, each with Pre → C1 → C2
   - Enables side-by-side cohort comparison at each timepoint

6. **Generate FeaturePlots:**
   - For each of 21 checkpoint genes:
     - Create 9-panel split FeaturePlot (3 cohorts × 3 timepoints)
     - Color scale: lightgrey (low) to red (high expression)
     - Quantile cutoffs: q10 (min) to q90 (max) for consistent scaling
     - Point ordering: High expression cells plotted on top (order=TRUE)
   - Arrange in grid: Gene labels (rows) × Group_TimePoint panels (columns)
   - Save as large PDF (width=20, height=2×21 genes)

7. **Export Underlying Data:**
   - Save FeaturePlot data to RDS file for downstream processing
   - Includes: UMAP coordinates, expression values, split labels
   - Structure: One row per cell per gene

**SECTION 2: 2D Density Visualization - Figure 8c (Lines 256-358)**

**Workflow:**

8. **Load FeaturePlot Data:**
   - Read RDS file generated in Section 1
   - Contains pre-computed UMAP coords and expression values

9. **Define `plot_2d_density()` Function:**
   
   **For each checkpoint gene:**
   
   **Step 9a: Expression Filtering**
   - Subset data for target gene
   - Identify top 10% expressing cells (90th percentile threshold)
   - These cells get density contours (high expressers)
   - All cells shown as gray background (spatial context)
   
   **Step 9b: Parse Metadata**
   - Extract TimePoint from split_label (Pre, C1, C2)
   - Extract Group from split_label (Control, ShortTerm, LongTerm)
   
   **Step 9c: Generate Three Subplots (Pre, C1, C2)**
   - For each timepoint panel:
     - **Base layer**: All cells in gray (background, alpha=0.5)
     - **Density contours**: Top 10% expressers, colored by cohort
       - Yellow contours: Control cohort high expressers
       - Blue contours: Short-term survivor high expressers
       - Red contours: Long-term survivor high expressers
     - **Weighting**: Density weighted by expression level (not just cell count)
     - **Bins**: Adaptive (3-10 contour lines based on cell number)
   
   **Step 9d: Panel Aesthetics**
   - Coordinate flip: Match manuscript orientation
   - Fixed aspect ratio: Prevent distortion
   - Minimal theme: Clean visualization
   - Title: Gene name and timepoint

10. **Generate Final PDF (Figure 8c):**
    - One page per checkpoint gene (21 pages total)
    - Each page: Three subplots (Pre, C1, C2) side-by-side
    - Dimensions: 15" wide × 5" tall per page
    - Progressive timepoints enable temporal comparison

**Biological Interpretation:**

**Checkpoint Expression Patterns:**
- **Dense contours**: High concentration of checkpoint-expressing cells
- **Dispersed contours**: Widespread low-level expression
- **Cohort-specific contours**: Differential expression between survival groups
- **Temporal changes**: Treatment-induced upregulation or downregulation

**Expected Patterns:**

**PD-1 (PDCD1) - Primary Therapeutic Target:**
- **Baseline (Pre)**: Higher in Control/Short-term (more exhaustion)
- **Post-treatment (C1-C2)**: May decrease in responders (checkpoint blockade efficacy)
- **Long-term survivors**: Could show lower sustained PD-1 (less exhaustion)

**CTLA-4 - Secondary Checkpoint:**
- **Activated T cells**: High expression during activation
- **Regulatory T cells**: Constitutively high (immune suppression)
- **Temporal dynamics**: May change with treatment

**LAG3, TIM-3, TIGIT - Exhaustion Markers:**
- **Co-expression**: Multiple checkpoints indicate severe exhaustion
- **Spatial clustering**: Exhausted cells may localize to specific UMAP regions
- **Cohort differences**: Short-term survivors may show more exhaustion markers

**CD28, ICOS, 4-1BB - Costimulatory:**
- **Activation markers**: High in functional T cells
- **Memory formation**: ICOS important for T follicular helper differentiation
- **4-1BB**: Correlates with T cell expansion and persistence

**Density Weighting Significance:**
- Contours weighted by expression level (not just presence/absence)
- Identifies regions with highest checkpoint expression intensity
- More biologically meaningful than simple cell counts

**Timepoint Progression (Pre → C1 → C2):**
- **Pre**: Baseline checkpoint landscape before treatment
- **C1**: Early response to checkpoint blockade (1 cycle)
- **C2**: Sustained response after 2 treatment cycles
- Progressive panels reveal treatment-induced dynamics

**Cohort Stratification:**
- **Control (Yellow)**: Baseline checkpoint expression patterns
- **Short-term (Blue)**: May show higher exhaustion markers
- **Long-term (Red)**: May show more balanced checkpoint expression (less exhaustion, maintained costimulation)

**Spatial Patterns in UMAP:**
- **Effector regions**: May show high PD-1, LAG3, TIM-3 (active but exhausted)
- **Memory regions**: May show CD28, ICOS (functional memory)
- **Exhausted regions**: Dense multi-checkpoint expression clusters
- **Naive regions**: Low checkpoint expression (quiescent)

**Clinical Relevance:**
- Reveals mechanisms of checkpoint blockade resistance (persistent high PD-1)
- Identifies additional therapeutic targets (LAG3, TIM-3, TIGIT combination therapy)
- Stratifies patients by checkpoint expression profiles
- Guides biomarker development for response prediction

**Key Outputs:**

**Section 1 Outputs:**
- `immune_checkpoint_timepoint_splitted_featureplots.pdf`:
  - Grid layout: 21 genes (rows) × 9 group-timepoint combinations (columns)
  - Standard FeaturePlots with expression color scale
- `featureplot_data.rds`:
  - Complete data for downstream density analysis
  - Contains UMAP coords, expression values, split labels

**Section 2 Outputs (Figure 8c):**
- `ImmuneCheckpoint_2D_Density_from_FeaturePlotData.pdf`:
  - **Manuscript Figure 8c**
  - 21 pages (one per checkpoint gene)
  - Each page: Three timepoint panels (Pre, C1, C2)
  - Expression-weighted density contours colored by survival cohort

**Visualization Design:**
- **Gray background**: All T cells (spatial context)
- **Colored contours**: Top 10% expressers by cohort
- **Weighted density**: Contour intensity reflects expression level, not just cell count
- **Adaptive bins**: 3-10 contour lines based on cell number (optimal visualization)

**Parameters:**
- **Timepoints analyzed**: Pre, C1, C2 (early response focus)
- **Expression cutoffs**: q10-q90 (removes extreme outliers)
- **Top expressers**: 90th percentile threshold for density contours
- **Density weighting**: By expression level (captures intensity)
- **Point size**: 0.4 for background, 0.6 for FeaturePlots

**Dependencies:**
- T cell Seurat object: `MK_T_Cells_seurat_obj.RDS`
- Survival data: `patient_survival_data.csv` (cohort assignments)
- Libraries: Seurat, dplyr, ggplot2, cowplot

**Technical Notes:**
- Downsampling ensures balanced comparison (equal n per group × timepoint)
- Set seed (123) for reproducibility
- Adaptive bin number prevents over/under-smoothing
- Coordinate flip for manuscript orientation consistency

---

*Last Updated: 2024-12-12*