# Aggressive neuroblastomas start growing after infancy
# EpiClockNBL

## 1. Introduction

The code in this repository can be used to generate all figures and results found in Monyak et al. (2025). Preprocessing of some of the data sources is performed by Python scripts and R Markdown files, and generation of figures and results is done in Jupyter notebooks and R Markdown files.

Software requirements:
- Python 3
- R

## 2. Setup

Fork and clone this repository locally as normal.

### Python

Use a **bash** shell to run all scripts and Jupyter notebooks. To see what shell is running, use ```echo $SHELL```. Run the following line to append the path of the repository clone to the Python path:

```
repo_dir=/PATH/TO/REPO/PARENT/DIR/EpiClockNBL
echo "export PYTHONPATH=$PYTHONPATH:$repo_dir" >> ~/.bash_profile
```

### R

In your R environment, preferably Rstudio, run the following line and copy the path outputted:

```
file.path(Sys.getenv("R_HOME"), 'etc', 'Rprofile.site')
```

Append the following line to the file at the path outputted (create the file if necessary):

```
repo_dir <- '/PATH/TO/REPO/PARENT/DIR/EpiClockNBL'
```

replacing the code above appropriately with the path to the local repository clone.

####

Dependencies

- jsonlite
- BiocFileCache
- plyr

### Path variables

Copy `config.json.example` and rename to `config.json`. Then open the file and insert appropriate paths for the following attributes:
- **official_indir** — Path to a directory in an external file location (preferably Box) that can hold terabytes of data
- **Figure_data_dir** — Path to a directory in an external file location (preferably Box) that can hold terabytes of data

## 2. Supplementary Data Retrieval



## 3. Pipeline

### 2. TARGET Data Retrieval

First install the `TCGAbiolinks` package by running the following:
```
bash setup.sh
```

Then, to retrieve the TCGA data and generate the HTML output, run the following:
```
bash Run_Data_Prep.sh
```

It could take up to a few hours to run, though it will likely take less than 1 hour. This script should be run on a machine of at least 16 GB of memory.

### 3. Select fCpGs

Open the notebook Pipeline.ipynb inside "3. Select fCpGs" and run all cells.
