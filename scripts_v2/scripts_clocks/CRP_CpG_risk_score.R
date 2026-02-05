# CRP_CpG_risk_score by M. Wielscher et al 2022
# REF: https://www.nature.com/articles/s41467-022-29792-6
# REF: https://github.com/Mwielscher/EWAS_CRP

library(data.table)
rm(list=ls())
score_dir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/CRP_CpG_risk_score"
out_dir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE"
experimental_name <- "epimode"

cpgs_with_weights <- read.csv(paste0(score_dir,"/1151CpGs_with_effectMetaAnalysis.csv"), header = T)
colnames(cpgs_with_weights) <- c("CpG_Site","Coefficient")
cpgs_with_weights$Predictor <- "CRP_CpG_risk_score"
cpgs_with_weights$"Mean_Beta_Value" <- NA

## FOR THE USER - PLACE INPUT FILES HERE - CAN CHANGE readRDS or read.csv if you so wish 
load(file=paste0(out_dir,"/out_script0/",experimental_name,"_bVals_detectionPvalFlt.RData"))
data <- bVals
dim(data)

row_names <- data$ID_REF
data <- as.matrix(data[,-1])
rownames(data) <- row_names
head(data)
sexageinfo <- NULL
covariates <- NULL

## Start to Process Files 

message("1. Loading data") 

message("1.1 Loading Methylation data - rows to be CpGs and columns to be individuals") 
rm(cpgs)
cpgs <- cpgs_with_weights
length(cpgs$CpG_Site) # 1511
length(unique(cpgs$CpG_Site)) # 1511
## Check if Data needs to be Transposed

message("2. Quality Control and data Preparation") 

message("2.1 Checking if Row Names are CpG Sites") 

if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-t(data) 
}

message("2.2 Subsetting CpG sites to those required for Predictor Calculation") 

## Subset CpG sites to those present on list for predictors 

coef=data[intersect(rownames(data), cpgs$CpG_Site),]

message("2.3 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from Training Sample") 

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site

coef <- if(all(cpgs$CpG_Site %in% rownames(coef))) { message("All sites present"); coef } else if(nrow(coef)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef=rbind(coef,mat)
} 
all(cpgs$CpG_Site %in% rownames(coef))

message("2.4 Convert NA Values to Mean for each Probe") 

## Convert NAs to Mean Value for all individuals across each probe 
table(is.na(coef))

na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))
table(is.na(coef))

message("2.5 Checking if Beta Values are Present") 

## Conversion to Beta Values if M Values are likely present  

m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}

message("3. Calculating the Predictors") 

loop = unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp, na.rm = T)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = colSums(tmp2)
  }
} 
out$ID <- row.names(out)
out <- out[,c(ncol(out),1:(ncol(out)-1))] 

## combine Sex and Age information
if(is.null(sexageinfo)){
  out$'True Age' <- NA 
  out$Sex <- NA 
} else { 
  ids = out$ID
  sexageinfo = sexageinfo[match(ids, sexageinfo$ID),] 
  out$'True Age' <- sexageinfo$Age
  out$Sex <- sexageinfo$Sex
  message("4. Sex and Age Info Added")
}

## Correct order of output file 
head(out)
out <- out[,c(1,ncol(out),(ncol(out)-1), 2:c(ncol(out)-2))]


## combine covariates  information
if(is.null(covariates)){
  NULL 
} else { 
  ids = out$ID
  covariates = covariates[match(ids, covariates$ID),] 
  out <- merge(out, covariates, by = "ID")
  for(i in names(out)[5:10]){ 
    out[,i] <- resid(lm(as.formula(paste(paste(bt(i),"~"),paste(bt(paste(c(names(covariates)[-which(names(covariates) %in% "ID")]))), collapse = "+"))), na.action = na.exclude, data = out))
  }
  out <- out[,c(1:10)]
  message("5. Covariates")
}

out$Sex <- NULL
out$`True Age` <- NULL
#write.csv(out,paste0(out_dir,"/",experimental_name,"_CRP_CpG_risk_score_output.csv"))
save(out, file=paste0(out_dir,"/out_clocks/",experimental_name,"_CRP_CpG_risk_score_output.RData"))


