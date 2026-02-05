rm(list=ls())
library(ggpubr)
library(plyr)
##########
# load dict of variables to analyze
data_dict <- read.csv("/Users/k/data/Epimode/out_clocks/DNAmAgeCalcProject_DataDict_ToAnalyze.csv")

data_dict$CurrentName <- gsub("\\s","\\.",data_dict$CurrentName)
data_dict$CurrentName <- sub("\\%",".",data_dict$CurrentName)
data_dict$CurrentName <- sub("\\:","\\.",data_dict$CurrentName)
data_dict$CurrentName <- sub("\\(","\\.",data_dict$CurrentName)
data_dict$CurrentName <- sub("\\)","\\.",data_dict$CurrentName)
data_dict$CurrentName <- sub("\\.\\.","\\.",data_dict$CurrentName)

dnam_trait <- data_dict$CurrentName[data_dict$ToAnalyze==1 & data_dict$FieldType=="DNAm_Trait"]
dnam_tissue_immcells <- data_dict$CurrentName[data_dict$ToAnalyze==1 & (data_dict$FieldType=="DNAm_Tissue_ImmuneCells" | data_dict$FieldType=="DNAm_Tissue_BrainCells")]
dnam_stem_cell_div <- data_dict$CurrentName[data_dict$ToAnalyze==1 & data_dict$FieldType=="DNAm_StemCellDivisions"]
dnam_grimAge_comps <- data_dict$CurrentName[data_dict$ToAnalyze==1 & data_dict$FieldType=="DNAm_GrimAgeComponent"]
dnam_grimAge2_comps <- data_dict$CurrentName[data_dict$ToAnalyze==1 & data_dict$FieldType=="DNAm_GrimAge2Component"]
dnam_fitAge_comps <- data_dict$CurrentName[data_dict$ToAnalyze==1 & data_dict$FieldType=="DNAm_FitAgeComponent"]
dnam_ages <- c(data_dict$CurrentName[data_dict$ToAnalyze==1 & (data_dict$FieldType=="DNAm_Age" ) ],"PoAmProjector","DunedinPACE", "AltumAge")
# check :
sum(data_dict$ToAnalyze==1)+3 == length(dnam_ages)+length(dnam_fitAge_comps)+length(dnam_grimAge_comps)+length(dnam_stem_cell_div)+length(dnam_tissue_immcells)+length(dnam_trait)+length(dnam_grimAge2_comps)

# load Horvath clock output
df <- read.csv("/Users/k/data/Epimode/out_clocks/DNAmAgeCalcProject_15786_Results.csv")
colnames(df)[1] <- "SampleID"
colnames(df) <- gsub("\\s","\\.",colnames(df))
colnames(df) <- sub("\\%","",colnames(df))
colnames(df) <- sub("\\:","\\.",colnames(df))
colnames(df) <- sub("\\.\\.","\\.",colnames(df))
df$Sex <- as.factor(revalue(as.character(df$Sex), c("male"="M","female"="F")))
df$Subject <- sub("_.*","", df$SampleID)
df$Group <- gsub("[0-9]","", df$Subject)
df$Timepoint <- sub(".*_","", df$SampleID)
unique(df$Subject)

# check :
all(data_dict$CurrentName[data_dict$ToAnalyze==1] %in% colnames(df))

# ss <-  read.csv("/Users/k/data/Epimode/SampleSheetEpimode.csv")
# ss$SampleID <- as.factor(ss$SampleID)
# ss$Batch <- as.factor(ss$Sample_Plate)
# ss$Group <- as.factor(ss$Group)
# ss$Timepoint <- as.factor(ss$Timepoint)
# head(ss)

ss <-  read.csv("/Users/k/data/Epimode/Epimode_MOH_phenotypeInT2.csv") # considering only MOH group with the clinical phenotype assigned in T2
ss$Subject <- as.factor(ss$Subject)
ss$Female <- NULL
ss$Overuse <- as.factor(ss$Overuse)
ss$Type <- as.factor(ss$Type)
ss$MHD <- as.numeric(ss$MHD)
head(ss)

# select only moh samples
df_all <- df
df <- df[df$Group=="MOH",] # 80 samples
unique(df$Subject)

# merge with phenotype data
df <- merge(df, ss, by = "Subject") # 74 samples (missing T2 data for : "MOH16" "MOH17" "MOH19" "MOH21" "MOH24")

##### load DunedinPACE scores #####
dir_clocks <- "/Users/k/data/Epimode/out_clocks/"
exp_name <- "epimode"
# created by PoAmProjector.R ; named "out"
load(paste0(dir_clocks,exp_name,"_DunedinPaceScore.RData"))
head(out)
score <- as.data.frame(out)
score$SampleID <- as.factor(rownames(score))
colnames(score)[1] <- "DunedinPACE"
score$SampleID[1:5]
df <- merge(score, df, by = "SampleID")

##### load PoAmProjector scores #####
# created by PoAmProjector.R ; named "out"
load(paste0(dir_clocks,exp_name,"_PoAmProjectorScore.RData"))
head(out)

score <- as.data.frame(out)
score$SampleID <- as.factor(rownames(score))
colnames(score)[1] <- "PoAmProjector"
score$SampleID[1:5]
df <- merge(score, df, by = "SampleID")

##### load CRP_CpG_risk scores #####
# created by CRP_CpG_risk_score.R
load(paste0(dir_clocks,exp_name,"_CRP_CpG_risk_score_output.RData"))
head(out)

score <- as.data.frame(out)
score$SampleID <- as.factor(rownames(score))
colnames(score)
score$ID <- NULL
colnames(score)[1] <- "CRP_CpG_risk_score"
score$SampleID[1:5]

df <- merge(score, df, by = "SampleID")

##### load EpiScores #####
# created by getEpiScores.R
load(paste0(dir_clocks,exp_name,"_EpiScores_output.RData"))
head(out)

score <- as.data.frame(out)
score$SampleID <- as.factor(rownames(score))
score$SampleID 
colnames(score) <- gsub("\\s","\\.",colnames(score))
colnames(score) <- gsub("-","",colnames(score))
colnames(score) <- gsub("\\/","",colnames(score))
colnames(score) <- gsub("\\.x","",colnames(score))
colnames(score) <- sub("\\%",".",colnames(score))
colnames(score) <- sub("\\:","\\.",colnames(score))
colnames(score) <- sub("\\(","\\.",colnames(score))
colnames(score) <- sub("\\)","\\.",colnames(score))
colnames(score) <- sub("\\.\\.","\\.",colnames(score))
colnames(score)[which(colnames(score) %in% c("Body.Mass.Index","Body.Fat.","HDL.Cholesterol","Waist.Hip.Ratio","Epigenetic.Age.Zhang."))] <- paste0(colnames(score)[which(colnames(score) %in% c("Body.Mass.Index","Body.Fat.","HDL.Cholesterol","Waist.Hip.Ratio","Epigenetic.Age.Zhang."))], "_2")
colnames(score) <- sub("6Ckine","Ckine6",colnames(score))

score$SampleID[1:5]

episcores <- colnames(score)[which(colnames(score) != c("ID","SampleID"))] # 109 proteins + 7 additional traits
episcores <- c(episcores, "CRP_CpG_risk_score")
episcores

df <- merge(score, df, by = "SampleID")
colnames(df) <- gsub("\\s","\\.",colnames(df))
colnames(df)<- gsub("-","",colnames(df))
colnames(df) <- gsub("\\/","",colnames(df))
colnames(df) <- gsub("\\.x","",colnames(df))
colnames(df) <- sub("\\%",".",colnames(df))
colnames(df) <- sub("\\:","\\.",colnames(df))
colnames(df) <- sub("\\(","\\.",colnames(df))
colnames(df) <- sub("\\)","\\.",colnames(df))
colnames(df) <- sub("\\.\\.","\\.",colnames(df))
# check :
all(episcores %in% colnames(df))
episcores[which(!(episcores %in% colnames(df)))]

##### load MorganLevine's PCclocks #####
# created by MorganLevine_PCclocks.R
load(paste0(dir_clocks,exp_name,"_PCclocks_output.RData"))
head(out)

score <- as.data.frame(out[,c(-2,-3)])
pcclocks <- colnames(score)[which(colnames(score) != c("SampleID"))]
pcclocks # 14 clocks
score$SampleID[1:5]

df <- merge(score, df, by = "SampleID")

##### load AltumAge #####
load(paste0(dir_clocks,exp_name,"_finalAltumAge.RData"))
head(out)

score <- as.data.frame(out)
df <- merge(score, df, by = "SampleID")

table(df$Timepoint)
table(df$Timepoint, df$Overuse,useNA = "ifany")

# save(df, file = "/Users/k/data/Epimode/out_clocks/allEpigeneticClocks.RData" )
# save(dnam_ages,dnam_grimAge_comps,dnam_grimAge2_comps,dnam_fitAge_comps,dnam_tissue_immcells,dnam_trait,dnam_stem_cell_div, CRP_CpG_risk_score, episcores, pcclocks, altumage,
#      file = "/Users/k/data/Epimode/out_clocks/variablesY.RData" )
