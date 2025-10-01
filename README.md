# ProteusFold — Protein Folding Model (Notebook Documentation)

## File Directory of this Repository

### 1. Making atomic properties
- `make_atom_residue_info.py`

### 2. File for making unfolded protein
- `make_unfolded_protein_from_sequence.py`

### 3. Convert PDB files into CSV
- `make_dataframe_using_biopandas.py`

### 4. Combine unfolded and folded protein CSVs
- `combine_new_dataframe.py`

### 5. Combine all the dataframe CSVs into a single CSV
- `combine_all_dataframes.py`

### 6. Main notebook for ProteusFold
- `ProteusFold_final_version.ipynb`  
  *Please refer to the Table of Contents below for its internal structure.*


---

## Table of Contents — `ProteusFold_final_version.ipynb`

### Helper Functions
- [User's Choice](#users-choice)
- [Helper Function](#helper-function)
  - [Save checkpoint](#save-checkpoint)
  - [Load checkpoint](#load-checkpoint)
  - [Trigger download](#trigger-download)
- [Kabsch: Align folded to unfolded structure](#kabsch-align-folded-to-unfolded-structure) *(Run once — use the flag variable to run only once)*
- [CSV → PDB string function](#csv-to-pdb-string-function)
  - [Visualizing both unfolded and folded protein](#visualizing-unfolded-and-folded-protein)

---

### Tokenizations
- [Protein Tokenization](#protein-tokenization)
  - [Tokenization](#tokenization)
    - [Structural tokenization](#structural-tokenization)
      - [True coordinates](#true-coordinates)
      - [Detokenization: inverting tokens back](#detokenization-inverting-tokens-back)
      - [Structural encoder](#structural-encoder)
      - [Structural dataloader](#structural-dataloader)
      - [Declarations](#declarations)
      - [Train and testing](#train-and-testing)
      - [Final evaluation](#final-eval)
      - [Converting to structural tokens](#converting-to-structural-tokens)
      - [PCA](#pca)
    - [Atomic tokenization](#atomic-tokenization)
  - [Make sequence](#make-sequence)

---

### Translation (Main Model)
- [Translation (main model)](#translation-main-model)

---

### Atomic Properties Analysis
- [Atomic properties analysis](#atomic-properties-analysis)
  - [Phase 1: Collecting grad × input using embedding](#phase-1-collecting-grad-x-input-using-embedding)
  - [Recollecting the atomic properties](#recollecting-the-atomic-properties)
  - [Phase 2: Jacobian–Vector Product (JVP)](#phase-2-jacobian-vector-product-jvp)
  - [Visualizing grad × input (various)](#visualizing-grad-x-input-various)
  - [Visualizing grad × input of a single protein (interactive, py3Dmol)](#visualizing-grad-x-input-single-protein)

---

### Inferring & Converting Coordinates (Visualization)
- [Inferring and converting to original coordinates](#inferring-and-converting-to-original-coordinates) — various subsections performing this task
- [Visualizing protein](#visualizing-protein)

---

### Data Distribution
- [Data distribution](#data-distribution) — visualize chain and total sequence length

---

### Metric Calculation
- [Metrics calculation](#metrics-calculation)
  - [DockQ](#dockq)
  - [Per-chain metrics](#per-chain-metrics)
  - [Entire-protein metrics](#entire-protein-metrics)
  - [Analysing the metrics](#analysing-the-metrics)

---

### Explanation (Requires GPU for some steps)
- [Explanation](#explanation) *(Needs GPU for initial calculations; if precomputed, GPU not required)*
  - [Attention maps](#attention-maps)
  - [Visualizing attention maps (interactive, py3Dmol)](#visualizing-attention-maps)
  - [Gradient × Input](#gradient-x-input)
  - [Heat map of Gradient × Input](#heat-map-of-gradient-x-input)
  - [Visualizing Grad × Input (interactive, py3Dmol)](#visualizing-grad-x-input)
  - [Heat maps and other graphs of attention](#heat-maps-and-other-graphs-of-attention)
  - [Targeted heatmaps (focused on specific proteins)](#targeted-heatmaps)

---

### Location of Interest
- [Calculating anchor regions (locations of interest)](#calculating-anchor-regions)

---

### Positional Analysis
- [Positional analysis](#positional-analysis)
  - [Residue coverage of locations of interest](#residue-coverage-of-locations-of-interest)
  - [Distribution of the top-1 important residue](#distribution-of-the-top-1-important-residue)
  - [Distribution of importance across normalized position](#distribution-of-importance-across-normalized-position)
  - [Distribution of importance across raw position](#distribution-of-importance-across-raw-position)

---

### Comparison
- [Comparison](#comparison)

---

### Parameter Count
- [Parameter count](#parameter-count)

---

