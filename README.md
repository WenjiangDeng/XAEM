
# XAEM: a novel method for isoform quantification


## 1. Introduction
This document shows how to use XAEM to quantify isoform expression for multiple samples from RNA-seq reads.


Software requirements for XAEM:

- R version 3.3.0 or later with installed packages: foreach and doParallel
- C++11 compliant compiler (g++ >= 4.7)
- Reference transcriptome: XAEM requires a fasta file of transcript sequences as reference. XAEM supports all kinds of reference and annotation for any species. In our paper XAEM uses the UCSC hg19 reference. You can download the reference of transcripts here: [transcripts.fa.gz](http://fafner.meb.ki.se/biostatwiki/2018_XAEM/transcripts.fa.gz)

#### X matrix (design matrix) :  
X matrix is an essential object for bias correction and isoform quantification (see our paper for more details). For users working on human the X matrix can be downloaded here: [X_matrix.RData](https://github.com/WenjiangDeng/XAEM/raw/master/X_matrix.RData). For other species the X matrix will be added soon.

## 2. Download and installation
If you use the binary verion of XAEM (recommended):

- Download the latest binary version from XAEM website:
```sh
wget https://github.com/WenjiangDeng/XAEM/raw/master/XAEM-binary-0.1.0.tar.gz
```
- Uncompress to folder
```sh
tar -xzvf XAEM-binary-0.1.0.tar.gz
```
- Move to the XAEM_home directory and do configuration for XAEM
```sh
cd XAEM-binary-0.1.0
bash configure.sh
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/XAEM-binary-0.1.0/lib:$LD_LIBRARY_PATH
export PATH=/path/to/XAEM-binary-0.1.0/bin:$PATH
```
If you want to build XAEM from sources:

Download XAEM from XAEM website and move to XAEM_home directory
```sh
wget https://github.com/WenjiangDeng/XAEM/raw/master/XAEM-source-0.1.0.zip
unzip XAEM-source-0.1.0.zip
cd XAEM-source-0.1.0
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

## 3. Index
Index the reference file
```sh
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/transcripts.fa.gz
gunzip transcripts.fa.gz
TxIndexer -t /path/to/transcripts.fa -o /path/to/TxIndexer_idx
```
## 4. XAEM: step by step instruction and explanation
### 4.1 Construction of the X matrix (design matrix)
This step constructs the X matrix required by the XAEM pipeline. For users working in human the X can be downloaded here: [X_matrix.RData](https://github.com/WenjiangDeng/XAEM/raw/master/X_matrix.RData). It's recommended to make a project folder and put the file in that folder, e.g. /path/to/XAEM_project. The command is:
```sh
mkdir /path/to/XAEM_project
wget https://github.com/WenjiangDeng/XAEM/raw/master/X_matrix.RData -P /path/to/XAEM_project
```
The steps to construct the design matrix are:

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

### 4.2 Generating the equivalence class table
The command to generate equivalence class table for each sample is similar to "sailfish quant".  For example, we want to run XAEM for sample1 and sample2 with 4 cpus:
```sh
XAEM -i /path/to/TxIndexer_idx -l IU -1 s1_read1.fasta -2 s1_read2.fasta -p 4 -o /path/to/XAEM_project/eqc_sample1
XAEM -i /path/to/TxIndexer_idx -l IU -1 s2_read1.fasta -2 s2_read2.fasta -p 4 -o /path/to/XAEM_project/eqc_sample2
```
If the data is compressed in gz format. We can combine with gunzip for decompression on-fly:
```sh
XAEM -i /path/to/TxIndexer_idx -l IU -1 <(gunzip -c s1_read1.gz) -2 <(gunzip -c s1_read2.gz) -p 4 -o /path/to/XAEM_project/eqc_sample1
XAEM -i /path/to/TxIndexer_idx -l IU -1 <(gunzip -c s2_read1.gz) -2 <(gunzip -c s2_read2.gz) -p 4 -o /path/to/XAEM_project/eqc_sample2
```
### 4.3 Creating Y count matrix

After running XAEM there will be the output of equivalence class table for multiple samples. We then create the Y count matrix. For example, if we want to run XAEM parallelly using 8 cores, the command is:
```sh
Rscript Create_count_matrix.R workdir=/path/to/XAEM_project core=8
```
### 4.4 Updating the X matrix and isoform expression using AEM algorithm
When finish the construction of Y count matrix we then use the AEM algorithm to update the X matrix. The updated X matrix is then used to estimate the isoform expression. The command is:
```sh
Rscript AEM_update_X_beta.R workdir=/path/to/XAEM_project core=8
```
The output in this step will be saved in XAEM_isoform_expression.RData, which is the TPM value and raw read counts of multiple samples.

**Note**: In XAEM pipeline we provide an extra step to merge the paralogs within the updated X matrix (please see our paper Section 2.2.3 and Section 2.3). The new X matrix is then used to estimate the final isoform expression. The paralog merging step produces more accurate estimation but can yield different sets of isoforms between different projects. If you want to run XAEM with this step in your project, you can simply add the "Merge=TRUE" parameter and run the command as below:
```sh
Rscript AEM_update_X_beta.R workdir=/path/to/XAEM_project core=8 Merge=TRUE
```
## 5. A complete run of XAEM by copy and paste
This section shows the tutorial to run XAEM pipeline. We can test XAEM by just copy and paste of the example commands.

- Download the binary file of XAEM
```sh
mkdir tmp_test
cd tmp_test
wget https://github.com/WenjiangDeng/XAEM/raw/master/XAEM-binary-0.1.0.tar.gz
tar -xzvf XAEM-binary-0.1.0.tar.gz
cd XAEM-binary-0.1.0
bash configure.sh
export LD_LIBRARY_PATH=/path/to/XAEM-binary-0.1.0/lib:$LD_LIBRARY_PATH
export PATH=/path/to/XAEM-binary-0.1.0/bin:$PATH
```
- Download fasta file and index it
```sh
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/transcripts.fa.gz
gunzip transcripts.fa.gz
TxIndexer -t transcripts.fa -o TxIndexer_idx
```
- Download the X matrix and RNA-seq data of sample1 and sample2
```sh
mkdir XAEM_project
cd XAEM_project
wget https://github.com/WenjiangDeng/XAEM/raw/master/X_matrix.RData
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample1_read1.fasta.gz
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample1_read2.fasta.gz
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample2_read1.fasta.gz
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample2_read2.fasta.gz
cd ..
```
- Generate the eqclass table and Y count matrix
```sh
XAEM -i TxIndexer_idx -l IU -1 <(gunzip -c XAEM_project/sample1_read1.fasta.gz) -2 <(gunzip -c XAEM_project/sample1_read2.fasta.gz) -p 4 -o XAEM_project/eqc_sample1
XAEM -i TxIndexer_idx -l IU -1 <(gunzip -c XAEM_project/sample2_read1.fasta.gz) -2 <(gunzip -c XAEM_project/sample2_read2.fasta.gz) -p 4 -o XAEM_project/eqc_sample2
## R packages foreach and doParallel are required

Rscript Create_count_matrix.R workdir=$PWD/XAEM_project core=8
```
- Estimate isoform expression using AEM algorithm
```sh
Rscript AEM_update_X_beta.R workdir=$PWD/XAEM_project core=8
cd XAEM_project
```


Reference: tba

Please also visit our webiste at http://fafner.meb.ki.se/biostatwiki/ for more bioinformatics tools : )
