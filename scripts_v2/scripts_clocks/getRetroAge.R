# kmk : This script is based on REF Ndhlovu et al 2024 https://pmc.ncbi.nlm.nih.gov/articles/PMC11464121/#sec24
# Model coefficients were downloaded from : https://zenodo.org/records/11099870 on 7th November 2024

##### Load clock coefficients #####
rm(list=ls())

### Process samples run on EPICv1.0 ###

# compatible with EPICv1.0
load("/home/PERSONALE/katarzyn.kwiatkowsk2/RetroelementAge/Retroelement_AgeV1coefficients_withMeanBetaValue.RData")
cpgs <- model_coeffs_v1
int <- intercept_v1

##### Load data for the cohort to be analyzed #####

# SampleSheet with chrAge
SampleSheet <- fread("/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/SampleSheetEpimode.csv")

# Filtered and normalized B-values
out_dir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE"
experimental_name <- "epimode"
load(file=paste0(out_dir,"/out_script0/",experimental_name,"_bVals_detectionPvalFlt.RData"))
dim(bVals)

##### Prepare Data #####
row_names <- bVals$ID_REF
data <- as.matrix(bVals[,-1])
rownames(data) <- row_names
head(data)

# Check if rows are CpGs
if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-t(data) 
}

# Subset CpG sites to those present on list for predictors
length(intersect(rownames(data), cpgs$CpG_Site)) # 1205 of 1317 (91%) CpGs are present in dataset 
coef=data[intersect(rownames(data), cpgs$CpG_Site),] 

# Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site

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

# Input random NAs substituting them with probe mean value across all individuals
table(is.na(coef))

na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))
table(is.na(coef))

# Check if dataset is composed of B-values. if not - convert them to B-values
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}

##### Calculate Retro-Age #####

loop = unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs$CpG_Site[cpgs$Predictor %in% i]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=(colSums(tmp_coef$Coefficient*tmp, na.rm = T))+int
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient+int
    out[colnames(coef),i] = colSums(tmp2)
  }
} 

out$ID <- row.names(out)
out <- out[,c(ncol(out),1:(ncol(out)-1))] 

##### Add data on ChrAge, Sex and Group #####
ids = out$ID
sexageinfo = SampleSheet[match(ids, SampleSheet$SampleID),] 
out$Age <- sexageinfo$Age
out$Female <- sexageinfo$Female

# Reorder columns
head(out)
out <- out[,c(1,ncol(out),(ncol(out)-1), 2:c(ncol(out)-2))]
head(out)
outV1 <- out
colnames(outV1) <- sub("V1","",colnames(outV1))

##### Save outcome and finish up #####
save(out, file=paste0(out_dir,"/out_clocks/",experimental_name,"_RetroAge_output.RData"))

