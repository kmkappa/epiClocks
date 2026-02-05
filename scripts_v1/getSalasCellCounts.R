setwd("/home/PERSONALE/francesca.ferraresi6/Deconvolution/data")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("FlowSorted.Blood.EPIC", force=TRUE)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ExperimentHub", force=TRUE)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", force=TRUE)

## Step 1: Extract the mix samples
library(FlowSorted.Blood.EPIC)
library(ExperimentHub)

FlowSorted.Blood.EPIC <- libraryDataGet('FlowSorted.Blood.EPIC')
str(FlowSorted.Blood.EPIC)
head(FlowSorted.Blood.EPIC)

#############STEP2############
#Crea il tuo RGset con script0

library(minfi)
library(limma)
library(RColorBrewer)

pal <- brewer.pal(8,"Dark2")
experiment_label <- "cart_norm"

getwd()

baseDir<-"/home/PERSONALE/francesca.ferraresi6/idatEPIC_CART/EPIC_CART/idat"
list.files(baseDir, recursive = TRUE)

targets <- read.metharray.sheet(baseDir)

targets$Basename
targets[1:4,]

str(targets)

RGset <- read.metharray.exp(targets = targets, extended=T, force = T)
RGset

sampleNames(RGset) <- targets$Sample_Name
targets$Group <- as.factor(targets$Group)
targets$Age <- as.numeric(targets$Age)
targets$Sex <- as.factor(targets$Sex)
#targets$Outcome <- as.factor(targets$Outcome)
targets

head(RGset)

## Step 3: use your favorite package for deconvolution.

load("FlowSorted.BloodExtended.EPIC.RData")

load("IDOLOptimizedCpGsBloodExtended.rda")## Deconvolute the target data set 450K or EPIC blood DNA methylation.
## We recommend ONLY the IDOL method, the automatic method can lead to severe
## biases.


## We recommend using Noob processMethod = "preprocessNoob" in minfi for the
## target and reference datasets.
## Cell types included are "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
## "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", and "Treg"
## Use estimateCellCounts2 from FlowSorted.Blood.EPIC.
## Do not run with limited RAM the normalization step requires a big amount
## of memory resources. Use the parameters as specified below for
## reproducibility.
#

 prop_extCAR <- estimateCellCounts2(RGset,

     compositeCellType = "BloodExtended",

     processMethod = "preprocessNoob",

     probeSelect = "IDOL",

     cellTypes = c(

         "Bas", "Bmem", "Bnv",

         "CD4mem", "CD4nv",

         "CD8mem", "CD8nv", "Eos",

         "Mono", "Neu", "NK", "Treg"),

CustomCpGs =if(RGset@annotation[1]=="IlluminaHumanMethylationEPIC"){

IDOLOptimizedCpGsBloodExtended}else{IDOLOptimizedCpGsBloodExtended450k})

 perc_extCAR<-round(prop_extCAR$prop*100,1)

head(perc_extCAR)

save.image("Deconv12cell_CARTall.RData")

###Boxplot Deconvoluzione

##da lista a dataframe
df <- as.data.frame((prop_extCAR$prop))
str(df)
head(df)

write.csv(df, "Deconvoluzione12cell_CARTall.csv")

# BOXPLOT
library(ggplot2)

# Trasposta per avere i tipi cellulari come variabili e aggiungi una colonna per i tipi cellulari
df_transposed <- as.data.frame(t(df))
df_transposed$CellType <- rownames(df_transposed)

# Formato lungo
df_long <- reshape2::melt(df_transposed, id.vars = "CellType")
head(df_long)

# Rinomina le colonne per avere nomi più significativi
colnames(df_long) <- c("CellType", "Patient", "Value")

# Crea il boxplot con ggplot2
library(ggplot2)
ggplot(df_long, aes(x = CellType, y = Value)) +
  geom_boxplot() +
  labs(title = "Boxplot dei Tipi Cellulari", x = "Tipi Cellulari", y = "Predicted")

ggsave("boxplot_deconvoluzione12cell_CARINGR.pdf", plot = last_plot(), device = "pdf")






