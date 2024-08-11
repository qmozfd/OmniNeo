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
### 2.2 Docker
We provide a docker image (https://github.com/qmozfd/ProGeo-neo-v3.0) which contains all package dependencies. You need to install docker on your system in advance. Download the Dockerfile form (https://github.com/qmozfd/ProGeo-neo-v3.0), then the command docker build  qmozfd/ProGeo-neo-v3.0:v1 will pull the ProGeo-neo-v3.0 image into your local machine.
```
docker build -t ProGeo-neo-v3.0 .
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
The results of RNAseq data preprocessing, call mutation and HLA typing are in the “rna_result”, “mut_result” and “HLAtype” folders, respectively.

#### Module2: Generation of tumor-specific variant peptides
```
python model_nonconding_mutpep.py control_name_R1.fastq.gz case_name_R1.fastq.gz
# eg:
python model_nonconding_mutpep.py con_R1.fastq.gz case_R1.fastq.gz
```
The results of the generated mutant peptides are under the “mut_result” folder.

#### Module3: Database construction and variant peptide identification 
```
python model_prediction.py
```
The result files will be stored under the “ms_resultmqpar/combined/txt” folder.

#### Module4: Neoantigen prediction and selection
```
python model_MS.py
```
**notes:**
Input the HLA types predicted in 3.1 or other types that the user interested in when the system prompts:
"please input an HLA class I allele like 'HLA-A02:01' or multiple alleles like 'HLA-A02:01,HLA-B15:01,HLA-C01:02':".

Neoantigen predictions and filtering results will be stored in the “preneo” directory.
#### Module5: Database construction and variant peptide identification 
```
python model_filter.py
```
The result files will be stored under the “ms_resultmqpar/combined/txt” folder.
#### Module6: Database construction and variant peptide identification 
```
python model_mrna_seq.py
```
The result files will be stored under the “ms_resultmqpar/combined/txt” folder.

### 3.2 GUI
In addition to the command line, we also provide a GUI to run the tool. Users can select a working directory by clicking on "Choose Directory", submit data by clicking on "Choose File" and start analysis by clicking on "Start Analysis". The log popup will show you the status of each step. You need to run the following command to bring up the GUI.
```
python gui.py
```
![](GUI.jpg)

For detailed software installation and usage of ProGeo-neo-v3.0, please read the User's Manual.
For coding neoantigen detection, we previously established a proteogenomic pipeline ProGeo-neo, ProGeo-neo source code and documentation are available at https://github.com/kbvstmd/ProGeo-neo.
