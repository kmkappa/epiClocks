
# https://www.nature.com/articles/s41514-022-00085-y
# https://github.com/rsinghlab/AltumAge
# https://rdrr.io/bioc/MEAT/man/BMIQcalibration.html


library(readr)
library(tibble)
library(data.table)
library(SummarizedExperiment)
library(MEAT)
rm(list=ls())

infile <- "epimode"

# create csv with beta values as for horvath but with previous datMiniAnnotation3.csv file
dat0 <- fread("/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/epimode_bVals_detectionPvalFlt.csv")
dat0 <- as.data.frame(dat0)

datMiniAnnotation=read.csv("/home/PERSONALE/katarzyn.kwiatkowsk2/filtering_IlluminaArrays/datMiniAnnotation3.csv")

match1=match(datMiniAnnotation[,1], dat0[,1] )
match1[1:20]
dat0Reduced=dat0[match1,]
dat0Reduced[,1]=as.character(datMiniAnnotation[,1])

options(stringsAsFactors=F)
set.seed('12345')
#out.csv = paste0("/home/gs66/Desktop/epicPainNet_ClocksUp/out_script0/",infile,"_horvath.csv")

input=dat0Reduced

colnames(input)[1]='ProbeID'
ann=datMiniAnnotation
cpgs=input$ProbeID
check=is.element(ann$Name,cpgs)
table(check)

miss.cpg=ann$Name[!check]
input.subject=colnames(input)[-1]
nsubject=dim(input)[2]-1
nmiss.cpg=length(miss.cpg)
add=data.frame(matrix(data=NA,nrow=nmiss.cpg,ncol=dim(input)[2]))
names(add)=names(input)
add[,1]=miss.cpg

output=rbind(input,add)
output=subset(output,ProbeID %in% ann$Name)
cat('check my new input dimension\n')
print(dim(output))
dim(ann)[1]==dim(output)[1]

head(output[,1:3])
output$permu=sample(dim(output)[1])
output=output[order(output$permu),]
output$permu<-NULL
head(output[,1:3])
#write.table(output,out.csv,sep=',',row.names=F,quote=F)


# Load phenotypes
pheno_data <- read_csv("/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/SampleSheetEpimode.csv")
head(pheno_data)
pheno_data <- pheno_data[,c("SampleID","Age","Female")]
pheno_data$SampleID <- as.character(pheno_data$SampleID)
pdata <- column_to_rownames(pheno_data, "SampleID")
head(pdata)

# Load matrix of beta-values
input_data <- data.frame(output)
rownames(input_data) <- NULL

colnames(input_data)
colnames(input_data)[1] <- "ProbeID"

head(input_data)
input_data <- column_to_rownames(input_data, "ProbeID")

t_input_data <- transpose(input_data)
colnames(t_input_data) <- rownames(input_data)
t_input_data$SampleID <- colnames(input_data)

t_input_data[1:3,1:3]
counts <- column_to_rownames(t_input_data, "SampleID")
counts[1:3,1:3]

all(rownames(pdata)==colnames(t(counts)))
pdata <- pdata[colnames(t(counts)),]
all(rownames(pdata)==colnames(t(counts)))
# Create a SummarizedExperiment object to coordinate phenotypes and methylation into one object.
se_data <- SummarizedExperiment(assays=list(beta=t(counts)),
                                     colData=pdata)


# Run clean_beta() to clean the beta-matrix
se_clean <- se_data # no cleaning # 30084 cpg
#se_clean <- clean_beta(SE = se_data, version = "MEAT2.0") # 18748 cpg
#se_clean <- clean_beta(SE = se_data, version = "MEAT") # 19401 cpg

# Run BMIQcalibration() to calibrate the clean beta-matrix
se_calibrated <- se_clean # no calibration
#se_calibrated <- BMIQcalibration(SE = se_clean, version = "MEAT2.0")
#se_calibrated <- BMIQcalibration(SE = se_clean, version = "MEAT")

assay(se_calibrated)
class(assay(se_calibrated))
t_assay <- as.data.frame(t(assay(se_calibrated)))
t_assay[1:3,1:3]
t_assay$SampleID <- row.names(t_assay)

t_assay_with_age_sex <- merge(pheno_data[,c("SampleID","Age","Female")],t_assay, by = "SampleID")
t_assay_with_age_sex[1:3,1:3]
rownames(t_assay_with_age_sex) <- t_assay_with_age_sex$SampleID
t_assay_with_age_sex$SampleID <- NULL
t_assay_with_age_sex[1:3,1:3]
colnames(t_assay_with_age_sex)[1] <- "age"
colnames(t_assay_with_age_sex)[2] <- "female"
t_assay_with_age_sex[1:3,1:3]
dim(t_assay_with_age_sex)

#write_csv(t_assay_with_age,"/home/gs66/Desktop/BRIC_vaccineHBV/out_horvath_october2023/BMIQcalibrated_beta_values_meat2.csv")
#write_csv(t_assay_with_age,"/home/gs66/Desktop/epicPainNet_ClocksUp/out_clocks_2023/BMIQcalibrated_beta_values_meat1.csv")
#write_csv(t_assay_with_age,"/home/gs66/Desktop/epicPainNet_ClocksUp/out_clocks_2023/AltumAge_output_preds_noBMIQ.csv")
write.csv(t_assay_with_age_sex,"/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/input_file_for_AltumAge_noBMIQ.csv")

# go to AltumAge.py script to calculate the predictions #
# use this last piece to create finalAltumAge_BricRoma.RData :
altumage <- read.csv("/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/epimode_AltumAge_output_preds_noBMIQ.csv", header = F)
colnames(altumage)[1] <- "AltumAge"
t_assay_with_age_sex$AltumAge <- altumage$AltumAge
t_assay_with_age_sex$SampleID <- rownames(t_assay_with_age_sex)
out <- t_assay_with_age_sex[,c("SampleID","AltumAge")]
save(out, file="/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_clocks/epimode_finalAltumAge_epicPainNet.RData")


