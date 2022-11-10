
# XAEM: a novel method for isoform quantification


## 1. Introduction
This vignette shows how to use XAEM to quantify isoform expression from RNA-seq data for **multiple samples**.

Software requirements for XAEM:

- R version 3.3.0 or later with installed packages: foreach and doParallel
- C++11 compliant compiler (g++ >= 4.7)
- XAEM is currently tested in Linux OS environment
#### Reference transcriptome: 
XAEM requires a fasta file of transcript sequences as reference. XAEM supports all kinds of reference and annotation for any species. In our paper we use the UCSC hg19 reference. You can download the reference of transcripts here: [transcripts.fa.gz](https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/transcripts.fa.gz)

#### X matrix (design matrix) :  
X matrix is an essential object for bias correction and isoform quantification (see our paper for more details). For users working on human, the X matrix can be downloaded here: [X_matrix.RData](https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/X_matrix.RData). For other species the X matrix will be added soon.

## 2. Download and installation
If you use the binary verion of XAEM (recommended):

- Download the latest binary version from XAEM website:
```sh
wget https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/XAEM-binary-0.1.1.tar.gz
```
- Uncompress to folder
```sh
tar -xzvf XAEM-binary-0.1.1.tar.gz
```
- Move to the XAEM_home directory and do configuration for XAEM
```sh
cd XAEM-binary-0.1.1
bash configure.sh
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/XAEM-binary-0.1.1/lib:$LD_LIBRARY_PATH
export PATH=/path/to/XAEM-binary-0.1.1/bin:$PATH
```
If you want to build XAEM from sources:

Download XAEM from XAEM website and move to XAEM_home directory
```sh
wget https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/XAEM-source-0.1.1.tar_.gz
unzip XAEM-source-0.1.1.zip
cd XAEM-source-0.1.1
bash configure.sh
```
XAEM requires information of flags from Sailfish including DFETCH_BOOST, DBOOST_ROOT, DTBB_INSTALL_DIR and DCMAKE_INSTALL_PREFIX. Please refer to the Sailfish website for more details of these flags.
Do installation by the following command:
```sh
DBOOST_ROOT=/path/to/boostDir/ DTBB_INSTALL_DIR=/path/to/tbbDir/ DCMAKE_INSTALL_PREFIX=/path/to/expectedBuildDir bash install.sh
After the installation is finished, remember to add the paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
export LD_LIBRARY_PATH=/path/to/expectedBuildDir/lib:$LD_LIBRARY_PATH
export PATH=/path/to/expectedBuildDir/bin:$PATH
```
Do not forget to replace "/path/to/" by your local path.

## 3. Preparation for the annotation reference
### 3.1 Indexing transcripts
Using TxIndexer to index the transcript sequences in the reference file (transcripts.fa). For example:
```sh
wget https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/transcripts.fa.gz
gunzip transcripts.fa.gz
TxIndexer -t /path/to/transcripts.fa -o /path/to/TxIndexer_idx
```
### 3.2 Construction of the X matrix (design matrix)
This step constructs the X matrix required by the XAEM pipeline. For users working in human the X can be downloaded here: [X_matrix.RData](https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/X_matrix.RData). It's recommended to make a project folder and put X matrix in that folder, e.g. /path/to/XAEM_project. The command is:
```sh
mkdir /path/to/XAEM_project
wget https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/X_matrix.RData -P /path/to/XAEM_project
```
The steps to construct your own design matrix are:

- Generate simulated RNA-seq data using the R package polyester

```sh
#### R package polyester and Biostrings are required
Rscript XAEM_home/R/genPolyesterSimulation.R /path/to/transcripts.fa /path/to/design_matrix
```
- Run GenTC to generate Transcript Cluster (TC) using simulated data. GenTC will generate an eqClass.txt file as the input for next step.
```sh
GenTC -i /path/to/TxIndexer_idx -l IU -1 /path/to/design_matrix/sample_01_1.fasta -2 /path/to/design_matrix/sample_01_2.fasta -p 8 -o /path/to/design_matrix
```
- Create the design matrix using the eqClass.txt from last step. The design matrix will be saved in X_matrix.RData. "H=0.025" is the threshold to filter false positive neighbors in each X matrix. (Please see our paper Section 2.2.1)
```sh
Rscript XAEM_home/R/buildCRP.R in=/path/to/design_matrix/eqClass.txt out=/path/to/design_matrix/X_matrix.RData H=0.025
```
## 4. XAEM: step by step instruction and explanation
Suppose we already created a working directory “XAEM_project” (/path/to/XAEM_project/) for quantification of transcripts, and the X_matrix.RData is saved in this folder.
### 4.1 Generating the equivalence class table
The command to generate equivalence class table for each sample is similar to "sailfish quant".  For example, we want to run XAEM for sample1 and sample2 with 8 cpus:
```sh
XAEM -i /path/to/TxIndexer_idx -l IU -1 s1_read1.fasta -2 s1_read2.fasta -p 8 -o /path/to/XAEM_project/eqc_sample1
XAEM -i /path/to/TxIndexer_idx -l IU -1 s2_read1.fasta -2 s2_read2.fasta -p 8 -o /path/to/XAEM_project/eqc_sample2
```
If the data is compressed in gz format. We can combine with gunzip for decompression on-fly:
```sh
XAEM -i /path/to/TxIndexer_idx -l IU -1 <(gunzip -c s1_read1.gz) -2 <(gunzip -c s1_read2.gz) -p 4 -o /path/to/XAEM_project/eqc_sample1
XAEM -i /path/to/TxIndexer_idx -l IU -1 <(gunzip -c s2_read1.gz) -2 <(gunzip -c s2_read2.gz) -p 4 -o /path/to/XAEM_project/eqc_sample2
```
### 4.2 Creating Y count matrix

After running XAEM there will be the output of the equivalence class table for **multiple samples**. We then create the Y count matrix. For example, if we want to run XAEM parallelly using 8 cores, the command is:
```sh
Rscript Create_count_matrix.R workdir=/path/to/XAEM_project core=8
```
### 4.3 Updating the X matrix and isoform expression using AEM algorithm
When finish the construction of Y count matrix, we use the AEM algorithm to update the X matrix. The updated X matrix is then used to estimate the transcript (isoform) expression. The command is:
```sh
Rscript AEM_update_X_beta.R workdir=/path/to/XAEM_project core=8 design.matrix=X_matrix.RData isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData merge.paralogs=FALSE isoform.method=average remove.ycount=TRUE
```
The output in this step will be saved in XAEM_isoform_expression.RData, which is the TPM value and raw read counts of multiple samples.

**Parameter setting**
- workdir: the path to working directory
- core: the number of cpu cores for parallel computing
- design.matrix: the path to the design matrix
- isoform.out (default=XAEM_isoform_expression.RData):  the output contains the estimated expression of individual transcripts, **where the paralogs are split into separate isoforms**. This file contains two objects: isoform_count and isoform_tpm for estimated counts and normalized values (TPM). The expression of the individual isoforms is calculated with the corresponding setting of parameter “isoform.method” below.
- isoform.method (default=average):  to report the expression of the individual members of a paralog as the average or total expression of the paralog set (value=average/total).
- paralog.out (default=XAEM_paralog_expression.RData): the output contains the **estimated expression of merged paralogs**. This file consists of two objects: XAEM_count and XAEM_tpm  for the estimated counts and normalized values (TPM). The standard error of the estimate is supplied in object XAEM_se stored in *.standard_error.RData.
- merge.paralogs (default=TRUE) (*****): the parameter to turn on/off (value=TRUE/FALSE) the paralog merging in XAEM. Please see the details of how to use this parameter in the note at the end of this section.
-remove.ycount (default=TRUE): to clean all data of Ycount after use.

**Note(*)**:In XAEM pipeline we provide this parameter (merge.paralog) to merge or not merge the paralogs within the updated X matrix (please see XAEM paper Section 2.2.3 and Section 2.3).  **Turning on (default) the paralog merging step produces a more accurate estimation. Turning off this step can produce the same sets of isoforms between different projects.**

## 5. A complete run of XAEM by copy and paste
This section presents a tutorial to run XAEM pipeline with a toy example. Suppose that input data contain two RNA-seq samples and server supplies 4 CPUs for computation. We can test XAEM by just copy and paste of the example commands.

- Download the binary file of XAEM
```sh
mkdir tmp_test
cd tmp_test
wget https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/XAEM-binary-0.1.1.tar.gz
tar -xzvf XAEM-binary-0.1.1.tar.gz
cd XAEM-binary-0.1.1
bash configure.sh
export LD_LIBRARY_PATH=/path/to/XAEM-binary-0.1.1/lib:$LD_LIBRARY_PATH
export PATH=/path/to/XAEM-binary-0.1.1/bin:$PATH
```
- Download fasta file and index it
```sh
wget https://github.com/WenjiangDeng/XAEM/releases/download/v0.1.1/transcripts.fa.gz
gunzip transcripts.fa.gz
TxIndexer -t transcripts.fa -o TxIndexer_idx
```
- Download the X matrix and RNA-seq data of sample1 and sample2
```sh
## Download input RNA-seq samples
# Create a XAEM project to save the data
mkdir XAEM_project
cd XAEM_project

# Download the RNA-seq data
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/XAEM_datasources/sample1_read1.fasta.gz
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/XAEM_datasources/sample1_read2.fasta.gz
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/XAEM_datasources/sample2_read1.fasta.gz
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/XAEM_datasources/sample2_read2.fasta.gz
cd ..
```
- Generate the eqclass table and Y count matrix
```sh
XAEM -i TxIndexer_idx -l IU -1 <(gunzip -c XAEM_project/sample1_read1.fasta.gz) -2 <(gunzip -c XAEM_project/sample1_read2.fasta.gz) -p 4 -o XAEM_project/eqc_sample1
XAEM -i TxIndexer_idx -l IU -1 <(gunzip -c XAEM_project/sample2_read1.fasta.gz) -2 <(gunzip -c XAEM_project/sample2_read2.fasta.gz) -p 4 -o XAEM_project/eqc_sample2
## R packages foreach and doParallel are required

Rscript $PWD/XAEM-binary-0.1.1/R/Create_count_matrix.R workdir=$PWD/XAEM_project core=$CPUNUM design.matrix=$PWD/X_matrix.RData
```
- Estimate isoform expression using AEM algorithm
```sh
Rscript $PWD/XAEM-binary-0.1.1/R/AEM_update_X_beta.R workdir=$PWD/XAEM_project core=$CPUNUM design.matrix=$PWD/X_matrix.RData isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData
```
The outputs are stored in the folder of “XAEM_project” including XAEM_isoform_expression.RData and XAEM_paralog_expression.RData.

## 6. Dataset for differential-expression (DE) analysis
In XAEM paper we have used the RNA-seq data from the breast cancer cell line (MDA-MB-231) for DE analysis. Since the original data is generated by our collaborators in Mayo Clinic and not published yet, we provide the equivalence class table by running the read-alignment tool Rapmap, which is the same mapper of Salmon and totally independent from XAEM algorithm. We also prepare the R scripts and the guide to replicate the DE analysis results in the paper. The dataset can be downloaded from Google drive [using the link here](https://drive.google.com/drive/folders/1vPGo4fpl7NC4_qzZCWsqK7yNRUahZf9p?usp=sharing).

Reference: [Alternating EM algorithm for a bilinear model in isoform quantification from RNA-seq data. Bioinformatics. 2020 Feb 1;36(3):805-812.](https://academic.oup.com/bioinformatics/article/36/3/805/5545974)

Please also visit our webiste at http://fafner.meb.ki.se/biostatwiki/ for more bioinformatics tools.
