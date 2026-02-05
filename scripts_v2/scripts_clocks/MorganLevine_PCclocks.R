##### PC-Clocks #####
# REF: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9586209/#SD2
# https://github.com/MorganLevineLab/PC-Clocks

#scp -r -J katarzyn.kwiatkowsk2@137.204.50.15 /Users/kk/Documents/pc/Home/PC-Clocks katarzyn.kwiatkowsk2@137.204.48.210:/home/PERSONALE/katarzyn.kwiatkowsk2/

rm(list=ls())

# Load the Functions you need to calculate the PCClocks (you need to change the path to the directory 
#       where you installed the code)

clocksDir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/PC-Clocks/" # dir with clocks

source(paste(clocksDir, "run_calcPCClocks.R", sep = ""))
source(paste(clocksDir, "run_calcPCClocks_Accel.R", sep = ""))

# Load my data
library(data.table)
library(dplyr)
out_dir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE"
experimental_name <- "epimode"
load(file=paste0(out_dir,"/out_script0/",experimental_name,"_bVals_detectionPvalFlt.RData"))
data <- bVals
dim(data)

datMeth2 <-as.matrix(t(data[,-1]))
colnames(datMeth2) <- data$ID_REF
rownames(datMeth2)

ss <-  read.csv(paste0(out_dir,"/SampleSheetEpimode.csv"))
head(ss)
dim(ss)
datPheno2 <- ss[,c("SampleID","Age","Female")]
datPheno2 <- datPheno2[match(rownames(datMeth2),as.character(datPheno2$SampleID)),]
rownames(datMeth2) == as.character(datPheno2$SampleID)

library(tibble)
datPheno2 <- as_data_frame(datPheno2)
class(datPheno2)

dim(datMeth2)
dim(datPheno2)
# IMPORTANT FORMATTING NOTE: If you are not using the example methylation and Pheno data, you will need to have specific
#     formatting. Please ensure that your Methylation dataframe/ matrix is of the methylation beta values and row names
#     are sample names, and column names are CpGs.
#     For the pheno data, ensure that the data frame/ matrix has rows as samples, and columns as whatever phenotype
#     variables you have/ wish. This can also include the original CpG clocks if you used the online Horvath calculator
#     as well. HOWEVER, the pheno matrix MUST have a column named "Age", and a column named "Female" (capital required),
#     especially if you want to calculate GrimAge and its components. Female = 1 is F and Female = 0 is M.
#
#     If you don't have this information, you can also just set it so Females is all 1 (all samples labeled female) and
#     all the same Age. Just know that this won't be accurate for PCGrimAge or components, and that you can't run
#     the acceleration calculations with calcPCClock_Accel.
#
#     The code below is going to ask if you would like to check the order of your datPheno and datMeth samples to ensure
#     they line up. For this to work, you will need to type the column name of datPheno with the names of the samples or 
#     'skip'.

# Get the PC Clocks values and the PC Clock Acceleration values #SampleID
out <- calcPCClocks(path_to_PCClocks_directory = clocksDir,
                                datMeth = datMeth2,
                                datPheno = datPheno2)

save(out, file=paste0(out_dir,"/out_clocks/",experimental_name,"_PCclocks_output.RData"))

head(out)







