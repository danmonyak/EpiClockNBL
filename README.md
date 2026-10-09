# EpiClockNBL

Code for **"Aggressive neuroblastomas start growing after infancy"** (Monyak et al., 2025).

This repository generates all figures and results in the paper. It estimates the mitotic age of neuroblastoma tumors from fluctuating CpG (fCpG) methylation in the TARGET cohort, validates it in an independent cohort (Henrich et al., GEO series GSE73515), and relates it to clinical variables, survival and gene expression.

Data preprocessing is done by Python modules and R Markdown files; figures and results are produced in Jupyter notebooks and R Markdown files.

## Contents

- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Setup](#setup)
- [Pipeline](#pipeline)
  - [1. Simulation](#1-simulation)
  - [2. TARGET Retrieval](#2-target-retrieval)
  - [3. Select fCpGs](#3-select-fcpgs)
  - [4. Process Supplementary Data](#4-process-supplementary-data)
  - [5. Gaussian Mixture Model](#5-gaussian-mixture-model)
  - [6. Main Analysis](#6-main-analysis)
  - [7. GSEA](#7-gsea)
- [Outputs and re-running](#outputs-and-re-running)

## Repository layout

| Path | Contents |
|---|---|
| `EpiClockNBL/` | Python package with the pipeline code that the notebooks call |
| `1. Simulation/` … `7. GSEA/` | Pipeline stages, run in order (see [Pipeline](#pipeline)) |
| `data/` | Small reference files (LUMP CpGs, Hallmark pathway names, protein-coding genes) |
| `config.json.example` | Template for the machine-specific `config.json` |

## Requirements

- Python 3.12 or later, with Jupyter
- R, with RStudio recommended for knitting the R Markdown files
- [GSEA desktop application](https://www.gsea-msigdb.org/gsea/downloads.jsp) (stage 7 only)
- Access to a SLURM cluster (stage 1 only)
- An external data directory with plenty of space; the TARGET download is large
- At least 16 GB of memory for the TARGET retrieval (stage 2)

### Python packages

- matplotlib
- numpy
- openpyxl
- pandas
- scikit-learn
- scipy
- seaborn
- tqdm

### R packages

**It is highly recommended to create two separate conda environments for R: one for most tasks, and one only for running rstan**, which can have installation and functionality issues.

General R environment:

- BiocFileCache
- DESeq2
- GSEABase
- GSVA
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19
- dplyr
- ggfortify
- ggplot2
- ggsci
- glmnet
- jsonlite
- kableExtra
- latex2exp
- msigdbr
- plyr
- remotes
- reshape2
- rmarkdown
- sesame
- sesameData
- survival
- tibble
- tidyr

A forked version of `TCGAbiolinks` is also needed; it is installed by a setup script in [stage 2](#2-target-retrieval).

rstan environment:

- bayesplot
- dplyr
- ggplot2
- jsonlite
- reshape2
- rmarkdown
- rstan

## Setup

Fork and clone this repository locally. Keep the clone's directory name as `EpiClockNBL`; the code uses it to locate the repository.

### 1. Python path

The repository directory must be on the Python path.

On Mac/Linux, use a **bash** shell to run all scripts and Jupyter notebooks (check with `echo $SHELL`). Run the following, replacing the template with the path to your clone:

```
repo_dir=/PATH/TO/REPO/PARENT/DIR/EpiClockNBL
echo "export PYTHONPATH=$PYTHONPATH:$repo_dir" >> ~/.bash_profile
```

On Windows, edit the system environment variables and append the path to your clone to the `PYTHONPATH` variable:

```
C:\PATH\TO\REPO\PARENT\DIR\EpiClockNBL
```

Always start Jupyter, and run Python scripts, from inside the repository directory (or one of its subdirectories).

### 2. R `repo_dir` variable

The R Markdown files expect a variable called `repo_dir`. **Repeat the following for each R environment.**

In R (preferably RStudio), run the following line and copy the path it prints:

```
file.path(Sys.getenv("R_HOME"), 'etc', 'Rprofile.site')
```

Append the following line to the file at that path (create the file if necessary), replacing the template with the path to your clone:

```
repo_dir <- '/PATH/TO/REPO/PARENT/DIR/EpiClockNBL'
```

or on Windows:

```
repo_dir <- 'C:\\PATH\\TO\\REPO\\PARENT\\DIR\\EpiClockNBL'
```

### 3. `config.json`

Copy `config.json.example`, rename the copy to `config.json`, and fill in:

- **official_indir** — path to a directory in an external file location (preferably Box) that can hold terabytes of data. All downloaded and intermediate data is stored here.
- **Figure_data_dir** — path to a directory in an external file location (preferably Box).
- **Windows** — `true` if working on a Windows machine, otherwise `false`.

On Windows, use double backslashes in the paths:

```
"official_indir": "C:\\PATH\\TO\\OFFICIAL_INDIR"
```

### 4. Data directories

Inside the *official_indir* directory, create two directories called `TARGET` and `Henrich`. `TARGET` holds the discovery cohort data and `Henrich` holds the validation cohort data.

## Pipeline

Run the stages in order. Each stage has a folder of the same name in the repository.

| Stage | What it does | How it is run |
|---|---|---|
| 1. Simulation | Simulates fCpG methylation in a growing tumor | SLURM jobs, then a notebook |
| 2. TARGET Retrieval | Downloads TARGET-NBL data and builds the clinical table | R Markdown, then a notebook |
| 3. Select fCpGs | Selects the fCpG (clock) sites | Notebook |
| 4. Process Supplementary Data | Prepares the validation cohort | Notebook |
| 5. Gaussian Mixture Model | Estimates the mitotic age of each tumor | Rscript (rstan environment) |
| 6. Main Analysis | Calendar ages, figures, survival | Notebooks and R Markdown |
| 7. GSEA | Gene set enrichment analysis | R Markdown and the GSEA application |

Stage 1 is independent of the others.

### 1. Simulation

The simulation needs a high-performance cluster because of its memory and compute requirements. *run_base_ensemble.py* simulates the fCpG states of a growing tumor up to a certain time point. The cells of the tumor are then divided equally into 200 parts, and each part is simulated separately for the remaining time by *run_split.py*.

From inside `1. Simulation`:

1. Create the directories that SLURM writes the job logs to:
   ```
   mkdir output_files error_files
   ```
2. Run the base simulation:
   ```
   sbatch run_base_ensemble.sh
   ```
3. When it has finished, schedule the *run_split.sh* jobs:
   ```
   sbatch job_scheduler_run_split.sh
   ```
4. Move the two output directories, `90_sites_NB_split_base` and `90_sites_NB_split_splitOutputs`, into a new directory called `sim_data` inside `1. Simulation`.
5. Open the notebook *Create Simulation Figures.ipynb* and run all cells. This recombines the data from the splits and creates the figures.

### 2. TARGET Retrieval

#### 2a. Data retrieval

1. From inside `2. TARGET Retrieval`, install the forked `TCGAbiolinks` package:
   ```
   bash setup.sh
   ```
   or on Windows, in PowerShell:
   ```
   .\setup.bat
   ```
2. **Windows only.** Link a virtual `P:` drive to the *official_indir* by running the following in PowerShell, replacing the template with the path to your *official_indir*:
   ```
   subst P: "C:\PATH\TO\OFFICIAL_INDIR"
   ```
   This works around a Windows path length limit in R that otherwise causes errors in *Data_Prep.Rmd*.
3. Open *Data_Prep.Rmd* and *Knit* the file. This retrieves the TARGET methylation, gene expression and clinical data (including the clinical supplement with MYCN status), saves it to *official_indir*/TARGET, and generates an HTML report.

   This usually takes less than 1 hour but can take a few hours. Use a machine with at least 16 GB of memory. With only 8 GB it can work, but it will take a few hours and the computer should not be used for anything else at the same time.

The TPM gene expression matrix is saved only as an `.rds` file at this point; the `.tsv` version is written in step 2c.

#### 2b. Data processing

Open the notebook *Data_Processing_Pipeline.ipynb* and run all cells to generate the annotated clinical table.

#### 2c. CIBERSORTx

1. From inside `2. TARGET Retrieval`, run:
   ```
   Rscript Save_Cibersort_Input.R
   ```
   This saves the TPM gene expression data as *cohort1.analysis_tumors.rnaseq_tpm.tsv* in *official_indir*/TARGET, restricted to the tumors in the analysis cohort that have gene expression data. On Windows, the virtual `P:` drive from step 2a must still be linked.
2. Upload *cohort1.analysis_tumors.rnaseq_tpm.tsv* to [CIBERSORTx](https://cibersortx.stanford.edu/) as the mixture file and run it with the LM22 signature matrix, B-mode batch correction, absolute mode and 500 permutations.
3. Save the results as *Cibersort_LM22_500perm_Analysis_Tumors.csv* in *official_indir*/TARGET.
4. **Windows only.** Unlink the virtual `P:` drive in PowerShell:
   ```
   subst P: /D
   ```

### 3. Select fCpGs

Open the notebook *Pipeline.ipynb* inside `3. Select fCpGs` and run all cells.

This selects the fCpG sites and saves the list to `3. Select fCpGs/outputs/NBL_Clock_CpGs.txt`. It also runs the sensitivity analysis, in which the tumors are randomly split into two halves and fCpGs are selected separately in each half (`outputs/split1` and `outputs/split2`).

### 4. Process Supplementary Data

1. Download the [series matrix file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73515/matrix/GSE73515_series_matrix.txt.gz) for the GSE73515 series into the *official_indir*/Henrich directory.
2. Decompress the file there, so that the directory contains *GSE73515_series_matrix.txt*.
3. Open the notebook *Pipeline.ipynb* inside `4. Process Supplementary Data` and run all cells.

### 5. Gaussian Mixture Model

A Gaussian mixture model is fit to each tumor's beta values in order to sample from the posterior of the mitotic age $\phi$. This requires the R library rstan, so activate the rstan environment first.

From inside `5. Gaussian Mixture Model`, run all four of the following. On Windows, run the `.bat` file of the same name instead of the `.sh` file.

TARGET cohort:
```
bash Run_GMM_STAN_TARGET.sh
```

Henrich cohort:
```
bash Run_GMM_STAN_Henrich.sh
```

Sensitivity analysis (the two halves from stage 3):

```
bash Run_GMM_STAN_TARGET_SPLIT1.sh
bash Run_GMM_STAN_TARGET_SPLIT2.sh
```

Each run takes a few hours. Progress is written to a file in the same folder, for example *TARGET.GMM_progress.txt* or *Henrich.GMM_progress.txt*, and the results are saved there as *TARGET.GMM_results.csv*, *Henrich.GMM_results.csv*, *SENSITIVITY_SPLIT1.TARGET.GMM_results.csv* and *SENSITIVITY_SPLIT2.TARGET.GMM_results.csv*.

### 6. Main Analysis

All files for this stage are inside `6. Main Analysis`. Run the steps in order.

#### a. Data processing

After all four GMM runs have finished, open the notebook *Pipeline.ipynb* and run all cells. This adds the mitotic ages to the clinical tables and saves the patient characteristics summary table.

#### b. Proliferation scores

1. Open *NBL_Proliferation_Scores.Rmd* in RStudio and click *Knit*. This calculates proliferation scores from the gene expression data.
2. Open *NBL_Predict_Prolieration_Score.Rmd* in RStudio and click *Knit*. This predicts proliferation scores from DNA methylation for the tumors without gene expression data.

#### c. Mitotic age and covariate analysis

Open the notebook *Make Figures.ipynb* and run all cells.

#### d. Tumor calendar age analysis

Open the notebook *Estimate Ages.ipynb* and run all cells. This requires the proliferation scores from step b.

#### e. Survival analysis

Open *NBL_survival.Rmd* in RStudio and click *Knit*.

#### f. PCA sensitivity analysis

Open the notebook *Sensitivity_Analysis_PCA.ipynb* and run all cells.

### 7. GSEA

All files for this stage are inside `7. GSEA`.

#### a. Data preprocessing

Open *GSEA_Data_Preprocessing.Rmd* in RStudio and knit the file. This preprocesses the gene expression and mitotic age data for use in GSEA. In the output HTML document, note the paths of the output `.cls` and `.gct` files.

#### b. Run GSEA

1. Download the GSEA software from https://www.gsea-msigdb.org/gsea/downloads.jsp
2. Open GSEA, navigate to *Load data*, and click *Browse for files*.
3. Load the `.cls` and `.gct` files.
4. Navigate to *Run GSEA* and select the following options.

   | Field | Value |
   |---|---|
   | Expression dataset | neuroblastoma_deseq_scale_factor |
   | Gene sets database | h.all.v2026.1.Hs.symbols.gmt |
   | Permutations | 1000 |
   | Phenotype labels | neuroblastoma_gmm_phi |
   | Collapse/Remap to gene symbols | Collapse |
   | Permutation type | gene_set |
   | Chip platform | Human_Ensembl_Gene_ID_MSigDB.v2026.1.Hs.chip |

   Under *Basic fields*:

   | Field | Value |
   |---|---|
   | Analysis name | original |
   | Metric for ranking genes | Spearman |
   | Save results in this folder | /PATH/TO/<official_indir>/TARGET/GSEA |

   Under *Advanced fields*:

   | Field | Value |
   |---|---|
   | Seed for permutation | 0 |

5. Select *Run* at the bottom of the page. GSEA should take about 90 seconds to run.

#### c. Create GSEA figure

1. Open *GSEA_Figure.Rmd* in RStudio.
2. Click the dropdown next to *Knit* and choose *Knit with Parameters*.
3. Look inside /PATH/TO/<official_indir>/TARGET/GSEA to find the GSEA output directory. Replace the *GSEA_output_dir_name* parameter with the name of that directory.
4. Look inside the GSEA output directory to find the names of the *pos* and *neg* files. They should be in this format: *gsea_report_for_neuroblastoma_gmm_phi_pos_1234567890.tsv* and *gsea_report_for_neuroblastoma_gmm_phi_neg_1234567890.tsv*. Replace the *pos_file* and *neg_file* parameters with the respective filenames.
5. Do not alter the other parameters.
6. Click *Knit*.

#### d. Sensitivity analysis

Repeat steps a to c for each half of the split, with the following changes (shown for `split1`; repeat with `split2`).

1. Knit *GSEA_Data_Preprocessing.Rmd* using *Knit with Parameters* and set:

   | Parameter | Value |
   |---|---|
   | output_gct_file | neuroblastoma_deseq_scale_factor.split1 |
   | output_cls_file | neuroblastoma_gmm_phi.split1 |
   | phi_var | phi.split1 |

2. When loading the files into GSEA, load *neuroblastoma_deseq_scale_factor.split1* and *neuroblastoma_gmm_phi.split1*.
3. Enter `split1` as the *Analysis name* under *Basic fields*.
4. When knitting *GSEA_Figure.Rmd*, enter `split1` for the *prefix* parameter.

## Outputs and re-running

- Figures are saved to `figures_revision` in the repository.
- Downloaded and intermediate data is saved under *official_indir*/TARGET and *official_indir*/Henrich.
- Several steps refuse to overwrite an existing output file if the new result differs from it, and raise an error instead. To regenerate such an output on purpose, delete the existing file first.
