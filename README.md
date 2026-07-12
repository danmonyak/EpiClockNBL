# Aggressive neuroblastomas start growing after infancy
# EpiClockNBL

## 1. Introduction

The code in this repository can be used to generate all figures and results found in Monyak et al. (2025). Preprocessing of some of the data sources is performed by Python scripts and R Markdown files, and generation of figures and results is done in Jupyter notebooks and R Markdown files.

Software requirements:
- Python 3
- R

## 2. Setup

Fork and clone this repository locally.

### Python

On Mac/Linux, use a **bash** shell to run all scripts and Jupyter notebooks. To see what shell is running, use ```echo $SHELL```. Run the following lines to append the path of the repository clone to the Python path:

```
repo_dir=/PATH/TO/REPO/PARENT/DIR/EpiClockNBL
echo "export PYTHONPATH=$PYTHONPATH:$repo_dir" >> ~/.bash_profile
```

(replacing the template with the appropriate path)

On Windows, edit system environmental variables to append to the PYTHONPATH variable:

```
C:\PATH\TO\REPO\PARENT\DIR\EpiClockNBL
```

### R

**It is highly recommended to create two separate conda environments for R, one for most tasks, and one for running rstan. Repeat the following for each environment:**

In R (preferably Rstudio), run the following line and copy the path outputted:

```
file.path(Sys.getenv("R_HOME"), 'etc', 'Rprofile.site')
```

Append the following line to the file at the path outputted (create the file if necessary):

```
repo_dir <- '/PATH/TO/REPO/PARENT/DIR/EpiClockNBL'
```

or if on Windows:

```
repo_dir <- 'C:\\PATH\\TO\\REPO\\PARENT\\DIR\\EpiClockNBL'
```

replacing the template above appropriately with the path to the local repository clone.

#### Dependencies (general R environment)

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
- reshape2
- sesame
- sesameData
- survival

#### Dependencies (rstan environment)

- bayesplot
- dplyr
- ggplot2
- jsonlite
- reshape2
- rstan

### Path variables

Copy `config.json.example` and rename to `config.json`. Then open the file and insert appropriate paths for the following attributes:
- **official_indir** — Path to a directory in an external file location (preferably Box) that can hold terabytes of data
- **Figure_data_dir** — Path to a directory in an external file location (preferably Box) that can hold terabytes of data
- **Windows** — *true* if working on a Windows machine, other *false*

**Note**: If on Windows, use double backslashes ("\\") in the path, i.e.:

```
"official_indir": "C:\\PATH\\TO\\OFFICIAL_INDIR"
```

### Directories

Inside the directory entered for *official_indir*, create directories called *TARGET* and *Henrich*. Henrich will hold the validation cohort data.

## 2. Supplementary Data Download

1. Download the [series matrix file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73515/matrix/GSE73515_series_matrix.txt.gz) for the GSE73515 series into the *official_indir*/Henrich directory.
2. Decompress/extract the file into the Henrich directory.

## 3. Pipeline

### 1. Simulation

### 2. TARGET Data Retrieval

First install the `TCGAbiolinks` package by running the following:
```
bash setup.sh
```

or if on Windows, run in powershell:

```
.\setup.bat
```

Next, if on Windows, run the following in powershell:

```
subst P: "C:\PATH\TO\OFFICIAL_INDIR"
```

replacing the template above appropriately with the path to the local repository clone. This links a virtual *P* drive to the *official_indir*, which solves a Windows-specific path length limit in R, that otherwise would cause errors when running *Data_Prep.Rmd*.


Then, to retrieve the TCGA data and generate the HTML output, open *Data_Prep.Rmd* inside *2. TARGET Retrieval* and *Knit* the file.

It could take up to a few hours to run, though it will likely take less than 1 hour. This script should be run on a machine of at least 16 GB of memory. If only 8 GB of memory is available, it can work but it will take a few hours and the computer should not be used for anything else at the same time.

If on Windows, run the following in powershell in order to unlink the virtual *P* drive:
```
subst P: /D
```

#### 2b. Data Processing

To generate the annotated clinical table, open the notebook Data_Processing_Pipeline.ipynb inside "2. TARGET Data Retrieval" and run all cells.

### 3. Select fCpGs

Open the notebook Pipeline.ipynb inside "3. Select fCpGs" and run all cells.

### 4. Process Supplementary Data

Open the notebook Pipeline.ipynb inside "3. Select fCpGs" and run all cells.

### 5. Gaussian Mixture Model

We will fit a GMM to the tumors' beta values, in order to sample from the posterior of the mitotic age $\phi$. This requires the r library rstan, and it is likely a good idea to create a new conda virtual environment specifically for this purpose, as rstan can have installation and functionality issues. Once you have activated the virtual environment in terminal, navigate to "5. Gaussian Mixture Model" and run the following in terminal:

TARGET cohort:
```
bash Run_GMM_STAN_TARGET.sh
```

Henrich cohort:
```
bash Run_GMM_STAN_Henrich.sh
```

For each cohort, the program will take a few hours to run, and you can check on the progress in these files: TARGET.GMM_progress.txt and Henrich.GMM_progress.txt.

### 6. Main Analysis

#### **a. Data processing**
After the GMM script has been run for both cohorts, open the notebook Pipeline.ipynb inside "6. Analysis" and run all cells.

#### **b. Tumor calendar age analysis**
Open the notebook *Make Figures.ipynb* inside "6. Analysis" and run all cells.

#### **c. Correlation and Miscellaneous Analysis**

Open the notebook *Estimate Ages.ipynb* inside "6. Analysis" and run all cells.

#### **d. Survival Analysis**

1. Open *Analysis/NBL_survival.Rmd* file in Rstudio.
2. Click *Knit*.

### 7. GSEA

**Data preprocessing**

Open *GSEA/GSEA_Data_Preprocessing.Rmd* file in Rstudio. Knit the file to preprocess the gene expression and mitotic age data for use in GSEA. In the output HTML document, note the paths of the output .cls and .gct files.

**Run GSEA**

1. Download the GSEA software from https://www.gsea-msigdb.org/gsea/downloads.jsp
2. Open GSEA, navigate to *Load data*, click *Browse for files*
3. Load the .cls and .gct files.
<!-- 4. Load the *h.all.v2026.1.Hs.symbols.gmt* file inside the *GSEA* folder in the git repository.  -->
5. Navigate to *Run GSEA*. Select the following options:

&emsp;&emsp;&emsp;&emsp;- Expression dataset: neuroblastoma_deseq_scale_factor<br>
&emsp;&emsp;&emsp;&emsp;- Gene sets databse: h.all.v2026.1.Hs.symbols.gmt<br>
&emsp;&emsp;&emsp;&emsp;- Permutations: 1000<br>
&emsp;&emsp;&emsp;&emsp;- Phenotype labels: neuroblastoma_gmm_phi<br>
&emsp;&emsp;&emsp;&emsp;- Collapse/Remap to gene symbols: Collapse<br>
&emsp;&emsp;&emsp;&emsp;- Permutation type: gene_set<br>
&emsp;&emsp;&emsp;&emsp;- Chip platform: Human_Ensembl_Gene_ID_MSigDB.v2026.1.Hs.chip<br>

&emsp;&emsp;&emsp;Under *Basic fields*:<br>
&emsp;&emsp;&emsp;&emsp;- Analysis name: original<br>
&emsp;&emsp;&emsp;&emsp;- Metric for ranking genes: Spearman<br>
&emsp;&emsp;&emsp;&emsp;- Save results in this folder: /PATH/TO/<official_indir>/TARGET/GSEA<br>

&emsp;&emsp;&emsp;Under *Advanced fields*:<br>
&emsp;&emsp;&emsp;&emsp;- Seed for permutation: 0<br>

6. Select *Run* at the bottom of the page. GSEA should take ~90 seconds to run.

**Create GSEA figure**

1. Open *GSEA/GSEA_Figure.Rmd* file in Rstudio.
2. Click the dropdown next to Knit and choose *Knit with Parameters*.
3. Navigate to /PATH/TO/<official_indir>/TARGET/GSEA to view the GSEA output directory. Replace the *GSEA_output_dir_name* parameter with the GSEA output directory name.
4. Navigate inside the to find the name of the *pos* and *neg* files. They should be in this format: *gsea_report_for_neuroblastoma_gmm_phi_pos_1234567890.tsv* and *gsea_report_for_neuroblastoma_gmm_phi_neg_1234567890.tsv*. Replace the *pos_file* and *neg_file* parameters with the respective filenames.
5. Do not alter the other parameters.
6. Click *Knit*.

**Sensitivity analysis**

Perform the same pipeline with a few changes.<br><br>
1. When running the *GSEA_Data_Preprocessing.Rmd* file, use *Knit with Parameters* and use the following:<br>
&emsp;&emsp;&emsp;&emsp;- output_gct_file: "neuroblastoma_deseq_scale_factor.split1"<br>
&emsp;&emsp;&emsp;&emsp;- output_cls_file: "neuroblastoma_gmm_phi.split1"<br>
&emsp;&emsp;&emsp;&emsp;- phi_var: "phi.split1"<br>

&emsp;&emsp;(repeat with "split2")

2. When loading the files into GSEA, load neuroblastoma_deseq_scale_factor.split1 and neuroblastoma_gmm_phi.split1 (repeat with "split2").
3. Enter "split1" or "split2" under *Analysis name* under *Basic fields*.
4. When Knitting the GSEA_Figure.Rmd file, enter "split1" (or "split2") for the *prefix* parameter.

