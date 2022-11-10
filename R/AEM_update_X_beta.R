## 27 Nov 2019 / Nghia:
# add standard error for the estimates
## 04 June 2019 / Nghia:
# add parameter: isoform.method=average/total to report the expression of the individual members of a paralog i) average (default) or ii) total from the paralog set
## 01 Apr 2019 / Wenjiang:
# add "merge.paralogs" parameter to turn on/off the paralog merging in XAEM. The default is off, which will generate the same set of isoforms between different projects. To turn it on, just add "merge.paralogs=TRUE"
# Example of command: Rscript buildCRP.R in=eqClass.txt isoform.out=X_matrix.RData merge.paralogs=TRUE

## Take the workdir and core arguments
workdir=NULL
core = 8 #default
merge.paralogs = TRUE ## default is to combine paralogs in the updated X to obtain the best performance
fout="XAEM_isoform_expression.RData"
foutr="XAEM_paralog_expression.RData"
design.matrix="X_matrix.RData"
isoform.method="average" #  "average" or "total"
remove.ycount=TRUE

saveSubset=FALSE #save singleton and paralogs 
noBiasResult=FALSE #export the results without bias correction
foutr_noBias="XAEM_noBiasCor.RData"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="workdir") workdir=as.character(res[2])
	if (res[1]=="core") core=as.numeric(res[2])
  if (res[1]=="design.matrix") design.matrix=as.character(res[2])
  if (res[1]=="isoform.out") fout=as.character(res[2])
  if (res[1]=="paralog.out") foutr=as.character(res[2])
	if (res[1]=="merge.paralogs") merge.paralogs=as.logical(res[2])
  if (res[1]=="isoform.method") isoform.method=as.character(res[2])
  if (res[1]=="remove.ycount") remove.ycount=as.logical(res[2])
}

cat("\n AEM_update_X_beta.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n workdir: ",workdir)
cat("\n core: ",core)
cat("\n design.matrix: ",design.matrix)
cat("\n isoform.out: ",fout)
cat("\n paralog.out: ",foutr)
cat("\n merge.paralogs: ",merge.paralogs)
cat("\n isoform.method: ",isoform.method)
cat("\n remove.ycount: ",remove.ycount)
cat("\n ----------------------------------------------------- ")


source("/path/to/R/Rsource.R")


options(stringsAsFactors=FALSE)
setwd(workdir)

#load input data
load(design.matrix)
load("Ycount.RData")
#set parallel
library(foreach)
library(doParallel)
ncores = detectCores()
nc = min(ncores,core)     # use 8 or 16 as needed!!
cl <- makePSOCKcluster(nc)   #
registerDoParallel(cl)

##### start from here
X.y= Y


fun = function(crpdat, maxiter.X=5, modify=TRUE){
  #xloc = grep('N', colnames(crpdat))
  xloc = which(colnames(crpdat) != "sample1")
  X0 = matrix(crpdat[,xloc], ncol=length(xloc))
  Ymat = crpdat[,-xloc]
  est = AEM(X0, Ymat, maxiter.X=maxiter.X, modify=modify)
  return(est)
}
X.y= Y

EST <- foreach(i= 1:length(X.y)) %dopar% fun(X.y[[i]], maxiter.X=5)

names(EST) = names(X.y)

#save(EST,file="Updated_X.Rdata")

####
#### Do not run add paralog step 
####
if(!merge.paralogs)
{
x.all = list()
X.y=Y
for(i in 1:length(X.y))
{
 x1=x2=NULL
 x1 = EST[[i]]$X
 x.y =  X.y[[i]]
 #xloc = grep('N', colnames(x.y))
 xloc = which(colnames(x.y) != "sample1")
 colnames(x1) = colnames(x.y)[xloc]
 x.all[[i]] = x1
}
}

####
#### run add paralog step with X from EST result
####
if(merge.paralogs)
{
x.all = list()
X.y=Y
for(i in 1:length(X.y))
{
 x1=x2=NULL
 x1 = EST[[i]]$X
 x.y =  X.y[[i]]
 #xloc = grep('N', colnames(x.y))
 xloc = which(colnames(x.y) != "sample1")
 colnames(x1) = colnames(x.y)[xloc]
 x.all[[i]] = try(ccrpfun(x1),silent = TRUE)
}
run.err=NULL
for(i in 1:length(X.y))
 if(class(x.all[[i]])== "try-error")
  run.err=c(run.err,i)

 for(i in run.err)
 {
  x1=x2=NULL
  x1 = EST[[i]]$X
  x.y =  X.y[[i]]
  #xloc = grep('N', colnames(x.y))
  xloc = which(colnames(x.y) != "sample1")
  colnames(x1) = colnames(x.y)[xloc]
  x.all[[i]] = ccrpfun(x1,clim=5) #clim shoule be smaller, such as 10
 }
}
#save(x.all, file="One_more_Collapsing_X.RData")

cat("\n...Estimation using AEM algorithm...\n")
X.y=Y
beta.all = list()# use the new X to calculate new beta
se.all = list()# keep standard error
for(i in 1:length(X.y))
{
 x2=NULL
 x.y =  X.y[[i]]
 #xloc = grep('N', colnames(x.y))
 xloc = which(colnames(x.y) != "sample1")
 y = x.y[,-xloc]
 x2 = x.all[[i]]
 beta1 = foreach(j=1:ncol(y)) %dopar% EM(x2,y[,j])
 beta2 = Reduce(rbind,beta1)
 rownames(beta2) = NULL
 beta.all = c(beta.all,list(beta2))

 if (merge.paralogs){ #compute standard error only if using merge.paralogs=TRUE
   s2=getSE(x2,beta2,y)
   se.all = c(se.all,list(s2))   
 }
}
if (saveSubset) save(beta.all,file="Beta_final_paralog.Rdata")

##### process singletons
## singletons
estfun = function(mat){
  CRP.y = mat$crpcount
  #sum(!(names(CRP.y)==names(CRP))) 
  ## cluster info
  npat = sapply(CRP,nrow)   # number of occupancy patterns per cluster
  table(npat)
  loc1 = which(npat==1)  # clusters with 1 pattern

  TC1 = sapply(CRP.y[loc1],function(x) x[1, ncol(x)])
   names(TC1)=names(CRP.y[loc1])
  ## single tx
  est.all =  c(TC1)
  return(est.all)
}

flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
result_est=NULL
for(id in 1:length(flist)){ # call crpcount()
  load(flist[id])
  est1 = estfun(y)# estimation step
  result_est = cbind(result_est,est1)
  #cat("sample ",id,'\n')
}

#samplename1 = gsub(pattern = "pq", replacement = "N", x = samplename1)
colnames(result_est)=samplename1

if (saveSubset) save(result_est,file="Est_result_Singletons.RData")

seq1 = Reduce(cbind,beta.all);seq1=t(seq1)
XAEM_count = rbind(result_est,seq1)

if (merge.paralogs){ #export standard error only if using merge.paralogs=TRUE
  se1=Reduce(cbind,se.all);se1=t(se1)
  result_se=sqrt(result_est) #use sqrt(ycount) for singletons
  XAEM_se=rbind(result_se,se1)

  foutr_se=foutr
  foutr_se=gsub(".RData","",foutr_se)
  foutr_se=paste(foutr_se,".standard_error.RData",sep="")
  save(XAEM_se, file=foutr_se)
}
##### done with the estimation

### collect the list of txnames from count data, more than 1 transcripts if the isoform is a paralog
txList = sapply(rownames(XAEM_count),function(x){return(unlist(strsplit(x," ")))})
txNum=sapply(txList,length)
### get median txlength for paralog
txLen=sapply(txList, function(x){
  pick=names(txlength) %in% x
  return(median(txlength[pick]))
})
### compute TPM
#normalize count to length
isoform_lenNorm=apply(XAEM_count,2,function(x)return(x/txLen))
libsize_lenNorm=apply(isoform_lenNorm,2,sum)
XAEM_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
XAEM_tpm=t(XAEM_tpm)
#keep information of raw output of XAEM
save(XAEM_count, XAEM_tpm,file=foutr)

## expand isoforms can not separated from CRP
paralogID=which(txNum >1)
paralog_count=XAEM_count[paralogID,]
paralog_tpm=XAEM_tpm[paralogID,]
#get mapping ID
matchID=sapply(c(1:length(paralogID)), function(x) rep(paralogID[x],txNum[paralogID][x]))
matchID=unlist(matchID)
#get isoform names
matchNames=sapply(c(1:length(paralogID)), function(x) unlist(strsplit(names(paralogID[x])," ")))
matchNames=unlist(matchNames)
#update to data XAEM_count
expandDat=matrix(0,ncol = ncol(XAEM_count), nrow = sum(txNum[paralogID]))
expandDat=XAEM_count[matchID,] #if (isoform.method=="total"):the counts of isoform members are equal to the count of paralog
if (isoform.method=="average"){ #the counts of isoform members are equal to the average count of paralog
  expandDat_txnum=txNum[match(names(matchID),names(txNum))]
  expandDat2=apply(cbind(expandDat_txnum,expandDat),1,function(x) x[-1]/x[1])
  expandDat=t(expandDat2)
}
rownames(expandDat)=matchNames
isoform_count=XAEM_count
isoform_count=rbind(isoform_count,expandDat)
isoform_count=isoform_count[-paralogID,]
##recalculate TPM
txLen2=txlength[match(rownames(isoform_count),names(txlength))]
isoform_lenNorm=apply(isoform_count,2,function(x)return(x/txLen2))
libsize_lenNorm=apply(isoform_lenNorm,2,sum)
isoform_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
isoform_tpm=t(isoform_tpm)

#export to file
pick=names(txlength) %in% rownames(isoform_count)
pick=which(!pick)
if(length(pick)>0){
  newDat=matrix(0,ncol = ncol(isoform_count), nrow = length(pick))
  rownames(newDat)=names(txlength)[pick]
  colnames(newDat)=colnames(isoform_count)
  isoform_count=rbind(isoform_count,newDat)
  isoform_tpm=rbind(isoform_tpm,newDat)
}
isoform_count=isoform_count[names(txlength),]
isoform_tpm=isoform_tpm[names(txlength),]

save(isoform_count,isoform_tpm,file=fout)


###### No bias correction 
if (noBiasResult){
  cat("\n Get results with no bias correction")
  #paralogs
  X.y= Y
  beta.all = list()
  for(i in 1:length(X.y))
  {
   x2=NULL
   x.y =  X.y[[i]]
   #xloc = grep('N', colnames(x.y))
   xloc = which(colnames(x.y) != "sample1")
   y = x.y[,-xloc]
   x2 = x.y[,xloc,drop=FALSE] ## X is not updated using AEM algorithm
    beta1 = foreach(j=1:ncol(y)) %dopar% EM(x2,y[,j])
   beta2 = Reduce(rbind,beta1)
   rownames(beta2) = NULL
  }
  #save(beta.all,file="Beta_final_paralog.RData")
  ## singletons
  estfun = function(mat){
    CRP.y = mat$crpcount
    sum(!(names(CRP.y)==names(CRP))) 
    ## cluster info
    npat = sapply(CRP,nrow)   # number of occupancy patterns per cluster
    table(npat)
    loc1 = which(npat==1)  # clusters with 1 pattern

    TC1 = sapply(CRP.y[loc1],function(x) x[1, ncol(x)])
     names(TC1)=names(CRP.y[loc1])
    ## single tx
    est.all =  c(TC1)
    return(est.all)
  }

  flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
  result_est=NULL
  for(id in 1:length(flist)){ # call crpcount()
    load(flist[id])
    est1 = estfun(y)# estimation step
    result_est = cbind(result_est,est1)
  }
  colnames(result_est)=samplename1
  seq1 = Reduce(cbind,beta.all);seq1=t(seq1)
  XAEM_count = rbind(result_est,seq1)

  ### collect the list of txnames from count data, more than 1 transcripts if the isoform is a paralog
  txList = sapply(rownames(XAEM_count),function(x){return(unlist(strsplit(x," ")))})
  txNum=sapply(txList,length)
  ### get median txlength for paralog
  txLen=sapply(txList, function(x){
    pick=names(txlength) %in% x
    return(median(txlength[pick]))
  })
  ### compute TPM
  #normalize count to length
  isoform_lenNorm=apply(XAEM_count,2,function(x)return(x/txLen))
  libsize_lenNorm=apply(isoform_lenNorm,2,sum)
  XAEM_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
  XAEM_tpm=t(XAEM_tpm)
  #keep information of raw output of XAEM
  save(XAEM_count, XAEM_tpm,file=foutr_noBias)
} #done for no bias


### clean Ycount
if (remove.ycount){
  system("rm -rf Ycount")
  system("rm -rf Ycount.RData")
}

cat("\n...Done...\n")
