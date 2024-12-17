# Blood risk biomarkers of leprosy patients

A brief description of what this project does and who it's for.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Requirements

For this project, you will need to use the following programming languages:
- [R 4.4.2](https://www.r-project.org/)
- [Python 3.12.8](https://www.python.org/)

The code within this repository is used on a Windows 11 device. It is expected that they work on other computers, however some paths and other device specific variables might need to be adjusted to that operating system

## Installation

Instructions on how to install and set up the project.

```bash
# Clone the repository
git clone https://github.com/IamAart/Leprosy_RNA.git

# Navigate to the project directory
cd Leprosy_RNA
```

Instructions on how to install all python packages in a virtual environment

```bash
# Create a virtual environment
python -m venv "./venv"

# Install dependencies
pip install -r "requirements.txt
```

Instructions on how to install all R packages

```bash
# install the requiRements package
source r_install.sh
```

## Usage

Within this repository, two main tasks are performed, namely Differential Gene Expression (DGE) analysis and using differentially expressed genes with Random Forest model. This is perfomed in order to find the best biomarker to use as identification of a risk of an individual progressing with leprosy

### Differential gene expression analysis

1. Set different variables in the `.env` file
    - `DGE_METHOD`, you wish to use to `"LimmaVoom", "DESeq2" or "EdgeR"` inside the 
    - `CUT_OFF`, you wish to use as p-value cut off that differentially expressed genes should be below
    - `LOG2_OFF`, you wish to use a fold change cut off that differentially expressed genes should be above or below -LOG2_OFF. Currently, 1.5 is used, which means a fold change difference of 50% between conditions
    - `BOOL_MDS_PLOT`, `"TRUE"` or `"FALSE"` to indicate whether MDS plots should be made during the analysis
2. Run the command `R -e 'source("./DGE_Analysis.R")
3. See the results in the folder `"./data/<DGE_METHOD>"` and the MDS plots (if set to TRUE) in the folder `"./data/Plots"`

### Random Forest
<!-- TODO: add this and plots -->

## Contributors

- Aart Rosmolen (s2548526)

## Contact

Aart Rosmolen - [a.rosmolen@umail.leidenuniv.nl](mailto:a.rosmolen@umail.leidenuniv.nl)

Project Link: [Leprosy_RNA](https://github.com/IamAart/Leprosy_RNA)