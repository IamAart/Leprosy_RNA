# Prediction of leprosy susceptibility from an RNA-seq population study

This repository is created in order to analyse data originating from a field trial, in which whole blood samples were collected from households of an individual that is diagnosed with leprosy.

This data will be used in this study to find genes that are differentially expressed, comparing progressors and non-progressors at the timepoint where neither showed any clinical symptoms. With these differentially expressed genes, a subgroup of genes will be identified with the use of machine learning models that can accurately predict the early development of leprosy. 

Multiple differential gene expression (DGE) analysis methods (DESeq2, EdgeR, LimmaVoom) are used on RNA-Seq data between progressors without symptoms and household contacts. Subsequently, the predictive performance of four different machine learning models (Random Forest (RF), Support Vector Machine (SVM), Linear Logistic Regression (LLR), and K-Nearest-Neighbours (kNN)) along with the performance of two different feature selection methods (Chi-Squared and Recursive Feature Elimination) was compared to find the optimal predictor to identify genes leprosy susceptibility genes from the DGE results.

## Disclaimer

The results shown in the corresponding paper may differ from those produced by the current version of the code in the repository. This discrepancy arises from the non-deterministic nature of the random forest classifier, as well as changes introduced during the integration of different code components, specifically the differential gene expression analysis and the machine learning analysis. 

## Table of Contents

- [Papers](#papers)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Papers
This repository is used for two different papers. 
<!-- TODO make links working -->
1. The first paper is based on the [Bachelor Thesis](link) written by Aart Rosmolen and Emma van Eerde. 
2. The second paper is based on the [Research](link) performed by Matheus Rogério Almeida et al. 

## Requirements

For this project, you will need to use the following programming languages:
- [R 4.4.2](https://www.r-project.org/)
- [Python 3.12.8](https://www.python.org/)

Additionally, it was performed on the following hardware specifications:
- **Operating System**: CentOS 7, Linux with the use of Elrepo 7
- **CPU**: Intel Xeon E5-2630v3 16 cores @ 2.40GHz (32 threads)
- **RAM**: 1.5 TB

## Installation

Instructions on how to install and set up the project.

```bash
# Clone the repository
git clone https://github.com/IamAart/Leprosy_RNA.git

# Navigate to the project directory
cd Leprosy_RNA
```

#### Instructions on how to install all python packages in a virtual environment

```bash
# Create a virtual environment
python -m venv "./venv"

# Install dependencies
pip install -r "requirements.txt
```

#### Instructions on how to install all R packages

```bash
# install the requiRements package
source r_install.sh
```
This installation of R can sometimes go wrong. In case this happens, you should manually make sure that all the libraries from the files `./R_BiocManager.txt` and `./R_requirements.txt` are installed/imported.

## Usage

Within this repository, two main tasks are performed, namely Differential Gene Expression (DGE) analysis and using differentially expressed genes with different machine learning models. How to use the code with these 2 different tasks will be explained in the next two sections

### Differential gene expression analysis

1. Set different variables in the `.env` file
    - `DGE_METHOD`, you wish to use to `"LimmaVoom", "DESeq2" or "EdgeR"` inside the 
    - `CUT_OFF`, you wish to use as p-value cut off that differentially expressed genes should be below
    - `LOG2_OFF`, you wish to use a fold change cut off that differentially expressed genes should be above or below -LOG2_OFF. Currently, 1.5 is used, which means a fold change difference of 50% between conditions
    - `BOOL_MDS_PLOT`, `"TRUE"` or `"FALSE"` to indicate whether MDS plots should be made during the analysis
2. Run the command 
    ```bash 
    R -e 'source("./DGE_Analysis.R")`
    ```
3. See the results in the folder `"./data/<DGE_METHOD>"` and the MDS plots (if set to TRUE) in the folder `"./data/Plots"`
4. From the results of the 3 DGE methods, you are able to generate a venn diagram file. These venndiagrams are based on each biotype (Coding, Non-coding or both) The command for this would be: 
    ```bash
    cd "./src"
    python "./venn_diagram.py"
    ```

### Machine Learning
If you have performed all the steps above and you have the venn diagram files, you are able to run the following commands:
1. In order to perform all the machine learning options: 
    ```bash
    cd "./src"
    python "./machine_learning.py"
    ```
2. To analyse the results that are saved in the folder `./data/Predictions/`, you can run the following command:
    ```bash
    cd "./src"
    python "./analyse_predictions"
    ```
    
    This file consists of multiple different parts, which can be targetted based on what is commented or not. Therefore, it is important to uncomment certain aspects before you run this command. What part you need to uncomment can be seen in the file itself.
3. In order to run the subset of the intersection between the two different methods, shown in the flowchart `./data/flowchart.png`, you will need to adjust the `.env` file. The changes are:
    - `SUBSET_MODEL` needs to be changed to one of the 4 different models used ("RF", "SVM", "LLR", "kNN")
    - `SUBSET_MIN` needs to be a string of an integer containing the minimum amount of genes 
    - `SUBSET_MAX` needs to be a string of an integer containing the maximum amount of genes 
    
    After setting these changes, you can redo `step 1 & 2` in order to run and analyse the intersection subset
4. To generate some visuals with the results you can make use of the file `plots.R`. This can be run with the following command:
    ```bash 
    R -e 'source("./plots.R")`
    ```
    The plots will be generated in the folder `./data/Plots`

## Contributors
- Aart Rosmolen
- Emma van Eerde 
- Matheus Rogério Almeida 

## Contact
If you have any questions or remarks regarding this repository, you can contact us:
- [Aart Rosmolen](mailto:arosmolen2500@gmail.com)
- [Emma van Eerde](mailto:emmarosavaneerde@gmail.com)
- [Matheus Rogério Almeida](mailto:m.rogerio_almeida@lumc.nl)

The two research papers can be found at 
<!-- TODO: Add these links -->
1. [Bachelor Thesis](link) written by Aart Rosmolen and Emma van Eerde. 
2. [Research](link) performed by Matheus Rogério Almeida et al. 

The project repository can be found at https://github.com/IamAart/Leprosy_RNA.