# ProGeo-neo-v3.0
ProGeo-neo-v3.0: a ProteoGenomics-based Neoantigen Prediction Pipeline in Noncoding Regions

## Overview of ProGeo-neo-v3.0 pipeline
![](pipeline.png)
## 1. Introduction


## 2. Installation
We provide two methods of installation.
### 2.1 Installation according to the User's Manual 
Running environment: ProGeo-neo-v3.0 requires a Linux operation system (centos7) with Python (V3.8), R (V4.0), Perl (V5.16) and Java (V1.7) installed.
Run the following codes before getting started.
```
cd ProGeo-neo-v3.0
bash start.sh
```

## 3.	Usage
### 3.1 command
#### Module1: Noncoding somatic variant calling and HLA typing
```
python model_wes.py control_name_R1.fastq.gz control_name_R2.fastq.gz case_name_R1.fastq.gz case_name_R2.fastq.gz
# eg:
python model_wes.py con_R1.fastq.gz con_R2.fastq.gz case_R1.fastq.gz case_R2.fastq.gz

python model_rnaseq.py control_name_R1.fastq.gz control_name_R2.fastq.gz case_name_R1.fastq.gz case_name_R2.fastq.gz
# eg:
python model_rnaseq.py con_R1.fastq.gz con_R2.fastq.gz case_R1.fastq.gz case_R2.fastq.gz

```
The results of RNAseq data preprocessing, call mutation and HLA typing are in the "rna_result", "mut_result" and "HLAtype" folders, respectively.

#### Module2: Generation of tumor-specific variant peptides
```
python model_nonconding_mutpep.py control_name_R1.fastq.gz case_name_R1.fastq.gz
# eg:
python model_nonconding_mutpep.py con_R1.fastq.gz case_R1.fastq.gz
```
The results of the generated mutant peptides are under the "mut_result" folder.

#### Module3: Database construction and variant peptide identification 
```
python model_prediction.py
```
The result files will be stored under the "preneo/" folder.

#### Module4: Neoantigen prediction and selection
```
python model_MS.py
```
**notes:**

Neoantigen predictions and filtering results will be stored in the "mass/" directory.
#### Module5: Database construction and variant peptide identification 
```
python model_filter.py
```
The result files will be stored under the "neofilter/" folder.
#### Module6: Database construction and variant peptide identification 
```
python model_mrna_seq.py
```
The result files will be stored under the "mrna_seq/" folder.



For detailed software installation and usage of ProGeo-neo-v3.0, please read the User's Manual.
For coding neoantigen detection, we previously established a proteogenomic pipeline ProGeo-neo, ProGeo-neo source code and documentation are available at https://github.com/kbvstmd/ProGeo-neo.
For Non-coding neoantigen detection, we previously established a proteogenomic pipeline PGNneo, PGNneo source code and documentation are available at https://github.com/tanxiaoxiu/PGNneo.
