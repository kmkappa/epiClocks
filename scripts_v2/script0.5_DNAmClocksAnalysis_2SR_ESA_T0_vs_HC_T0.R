rm(list=ls())
library(ggpubr)
library(plyr)
##########
# load data and function
load("/Users/k/data/Epimode/DNAmClocks/out_clocks/allEpigeneticClocks.RData")
load("/Users/k/data/Epimode/DNAmClocks/out_clocks/variablesY.RData")
source("/Users/k/data/Epimode/DNAmClocks/fun_get_2SR_onVarX.R")

##### 2SR model comparing MOH_T0_Episodic vs MOH_T0_Chronic correcting for sex #####
out_dir <- "/Users/k/data/Epimode/DNAmClocks/out_2SR_ESA_T0_vs_HC_T0"
if (file.exists(out_dir)){}else{dir.create(out_dir)}
# keep original df & select specific analysis_name to use for subsequent computations
df$Group <- paste0(df$Group,"_", df$Timepoint)
df_ori<- df;  analysis_name <- "epimode_ESA_vs_HC";

# select phenotypes of interests
to_sel <- c("ESA_T0","HC_T0")

df <- df[df$Group %in% to_sel,]
table(df$Group)
df$Group <- as.factor(df$Group)
any(is.na(df$Group))

##### DNAmAge_Ages #####
variablesY_subset <- "dnam_ages"
variablesY <- eval(parse(text=variablesY_subset))
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=F)

##### DNAmAge_GrimAgeComps #####
variablesY_subset <- "dnam_grimAge_comps"
variablesY <- eval(parse(text=variablesY_subset))
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=T)

##### DNAmAge_GrimAge2Comps #####
variablesY_subset <- "dnam_grimAge2_comps"
variablesY <- eval(parse(text=variablesY_subset))
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=F)

##### DNAmAge_FitAgeComps #####
variablesY_subset <- "dnam_fitAge_comps"
variablesY <- eval(parse(text=variablesY_subset))
variablesY <- variablesY[!(variablesY %in% c("DNAmGrip_noAge","DNAmGrip_wAge","DNAmFEV1_noAge"))] # impossible to create model for these variables - they have too many NAs that substituted outliers
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=F)

##### PCclocks #####
variablesY_subset <- "pcclocks"
variablesY <- eval(parse(text=variablesY_subset))
variablesY <- variablesY[variablesY!="PCLeptin"] # there are 24 PCLeptin NAs (that substituted outliers) preventing us to create model
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=F)

##### DNAmAge_tissue_immcells #####
variablesY_subset <- "dnam_tissue_immcells"
variablesY <- eval(parse(text=variablesY_subset))
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=T)

##### DNAmAge_Traits #####
variablesY_subset <- "dnam_trait"
variablesY <- eval(parse(text=variablesY_subset))
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=F)

##### DNAmAge_StemCellDiv #####
variablesY_subset <- "dnam_stem_cell_div"
variablesY <- eval(parse(text=variablesY_subset))
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=F)

##### EpiScores #####
variablesY_subset <- "episcores"
variablesY <- eval(parse(text=variablesY_subset))
VarX <- "Age"
VarCov1 <- "Sex"
VarCov2 <- NA
refGroup <- "HC_T0"
caseGroup1 = "ESA_T0"
out_df_name <- paste0("2SR_",analysis_name,"_",variablesY_subset)
get_2SR_onVarX(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=T)



