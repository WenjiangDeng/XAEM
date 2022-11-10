## Take the workdir and core arguments
workdir=NULL
design.matrix="X_matrix.RData"
core = 8 #default

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="workdir") workdir=as.character(res[2])
	if (res[1]=="core") core=as.numeric(res[2])
  if (res[1]=="design.matrix") design.matrix=as.character(res[2])
}

cat("\n Create_count_matrix.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n workdir: ",workdir)
cat("\n core: ",core)
cat("\n design.matrix: ",design.matrix)
cat("\n ----------------------------------------------------- ")


source("/path/to/R/Rsource.R")


options(stringsAsFactors=FALSE)
setwd(workdir)
load(design.matrix)

if(!dir.exists("Ycount")) dir.create("Ycount")
setwd(paste(workdir,"/Ycount",sep=""))
flist = list.files(workdir,pattern="eqClass.txt",recursive=TRUE,full.names = TRUE)
tx_length = flist[1]
#initialization for parallel computing
library(foreach)
library(doParallel)
registerDoParallel(cores=core)

res=foreach(id = 1:length(flist),.combine=c) %dopar% {
  y = NULL
  y = crpcount(flist[id])
  samplename = paste(y$samplename,".RData",sep="")

  save(y,file=samplename)
  return(flist[id])
}
res=NULL

flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
npat = sapply(CRP,nrow)   # number of occupancy patterns per clusterloc2 = which(npat>1) 
loc2 = which(npat>1) 
CCRP1 = CCRP[loc2]
Y=NULL
for(id in 1:length(flist)){
 cat("Merging results from sample ",flist[id],' ...\n')
 load(flist[id])
 if(id==1){
  Y = y[[1]][loc2]
 }
 if(id>1){
  y1 = y[[1]][loc2]
  for(i in 1:length(Y)){
   Y[[i]] = cbind(Y[[i]],sample1=y1[[i]][,'sample1'])
   }
  }
 }



samplename1 = NULL
for(id in 1:length(flist))
 {
  s.1 = strsplit(flist[id],"/")[[1]]
  s = s.1[length(s.1)]
  s = gsub(".RData","",s)
  samplename1 = c(samplename1,s)
 }

 setwd(workdir)

 for(i in 1:length(Y))
 {
 y2 = Y[[i]]
 xloc = which(colnames(y2) != "sample1")
 y2.1 = y2[,-xloc]
 y3 = cbind(CCRP1[[i]],y2.1)
 Y[[i]] = as.matrix(y3)
 }

save(Y,samplename1,tx_length,file='Ycount.RData')

cat("\n...Done...\n")
