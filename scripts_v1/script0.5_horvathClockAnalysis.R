library(ggpubr)

##########

# load Horvath clock output
df <- read.csv("/home/gs66/Desktop/BRIC_vaccineHBV/out_horvath/hbv_norm_bVals_noob_RAW_horvath.output.csv")
# load sample sheet
ss <-  read.csv("/home/gs66/Desktop/BRIC_vaccineHBV/SampleSheet.csv")
# merge both data frames
df_merged <- merge(df, ss[,c("SampleID","Group", "Sex", "Batch","BcellCount")], by = "SampleID")
# select DNAm biomarkers of interest
horvathClockVariables <- c("DNAmAge","DNAmAgeHannum","DNAmAgeSkinBloodClock","DNAmPhenoAge","DNAmGrimAge", "DNAmTL",
                           "DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPAI1",  "DNAmTIMP1",
                           "DNAmPACKYRS","CD8T","CD4T","CD8.naive","CD4.naive","CD8pCD28nCD45RAn","NK","DunedinPoAm")

# 2SR model comparing NR vs R correcting for sex
out_dir <- "/home/gs66/Desktop/BRIC_vaccineHBV/out_horvath/out_2SR"
VarX <- "Age"
VarCov <- "Sex"
variablesY <- horvathClockVariables
refGroup <- "R"
caseGroup1 = "NR"
out_df_name <- "2SR_HBV_RvsNR"
get_2SR_onVarX(df_merged, VarX , VarCov, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot=F, res_lineplot=F, boxplot=F)

##### FUNCTIONS : Residual 2SR model #####
# test :
#var_i <- "DNAmAge"; refGroup <- "R"; caseGroup1 = "NR"; VarX <- "Age"; variablesY <- horvathClockVariables;

# model correcting for biological sex
get_2SR_onVarX <- function(df_noOutliers, VarX , VarCov, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot, res_lineplot, boxplot){
  AvrgGr1 <- c(); AvrgCTRL <- c();
  MedianGr1 <- c(); MedianCTRL <- c();
  Pval_Gr1_CTRL <- c();
  for(var_i in variablesY){
    if(is.na(VarCov)){
      lmCTRL <- lm(eval(parse(text = var_i)) ~ eval(parse(text = VarX)), data = df_noOutliers[df_noOutliers$Group== refGroup,])
    }
    else{lmCTRL <- lm(eval(parse(text = var_i)) ~ eval(parse(text = VarX)) + eval(parse(text = VarCov)), data = df_noOutliers[df_noOutliers$Group== refGroup,])}
    
    predCTRL <- predict(lmCTRL, df_noOutliers[df_noOutliers$Group==refGroup,])
    predGr1 <- predict(lmCTRL, df_noOutliers[df_noOutliers$Group==caseGroup1,])
    resCTRL <- df_noOutliers[df_noOutliers$Group==refGroup, var_i] - predCTRL
    resGr1 <- df_noOutliers[df_noOutliers$Group==caseGroup1, var_i] - predGr1
    df_noOutliers$residuals[df_noOutliers$Group==refGroup] <- resCTRL
    df_noOutliers$residuals[df_noOutliers$Group==caseGroup1] <- resGr1
    
    out <- t.test(df_noOutliers$residuals[df_noOutliers$Group==caseGroup1],df_noOutliers$residuals[df_noOutliers$Group==refGroup])
    Pval_Gr1_CTRL <- c(Pval_Gr1_CTRL,out$p.value); AvrgGr1 <- c(AvrgGr1,out$estimate[1]); AvrgCTRL <- c(AvrgCTRL,out$estimate[2]);
    MedianGr1 <- c(MedianGr1,median(df_noOutliers[df_noOutliers$Group==caseGroup1, var_i],na.rm = T));
    MedianCTRL <- c(MedianCTRL,median(df_noOutliers[df_noOutliers$Group==refGroup, var_i],na.rm = T));
    
    df_noOutliers$residuals <- NA
    df_noOutliers$residuals[df_noOutliers$Group==refGroup] <- resCTRL
    df_noOutliers$residuals[df_noOutliers$Group==caseGroup1] <- resGr1
    
    # Visualizations
    # age_lineplot
    if(age_lineplot==T){
      i <- which(variablesY == var_i)
      pdf(paste0(out_dir,"/age_lineplot_",out_df_name, "_", as.character(var_i),".pdf"), width = 7, height = 6)
      p <- qplot(Age, eval(parse(text = var_i)) , data= df_noOutliers, color = Group, size=I(1.5))+
        geom_smooth(method="lm", se=T) +
        labs(x="ChrAge" , y = paste0("",as.character(var_i)),
             subtitle = paste0("p-value ",caseGroup1," vs ",refGroup," = ", as.character(format(round(Pval_Gr1_CTRL[i],3),nsmall=3))) ) +
        ggtitle(paste0("ChrAge and ",as.character(var_i)," by Phenotype Group")) +
        scale_color_discrete(name = "Phenotype") +
        scale_linetype_manual(values=c("solid"))
      print(p)
      dev.off()
    }
    
    # res_lineplot
    if(res_lineplot==T){
      i <- which(variablesY == var_i)
      pdf(paste0(out_dir,"/lineplot_",out_df_name, "_", as.character(var_i),".pdf"), width = 7, height = 6)
      p <- qplot(Age, residuals , data= df_noOutliers, color = Group, size=I(1.5))+
        geom_smooth(method="lm", se=T) +
        labs(x="ChrAge" , y = paste0("residuals_",as.character(var_i)),
             subtitle = paste0("p-value ",caseGroup1," vs ",refGroup," = ", as.character(format(round(Pval_Gr1_CTRL[i],3),nsmall=3))) ) +
        ggtitle(paste0("ChrAge and res_",as.character(var_i)," by Phenotype Group")) +
        scale_color_discrete(name = "Phenotype") +
        scale_linetype_manual(values=c("solid"))
      print(p)
      dev.off()
    }
    
    #boxplot
    if(boxplot==T){
      i <- which(variablesY == var_i)
      pdf(paste0(out_dir,"/boxplot_",out_df_name, "_", as.character(var_i),".pdf"), width = 4, height = 6)
      p <- ggboxplot(df_noOutliers[(df_noOutliers$Group==refGroup | df_noOutliers$Group==caseGroup1),], y = 'residuals' , x= 'Group', add= 'jitter',color = 'Group',
                     #order = c("AMC", "GGD"),
                     main= paste0(as.character(var_i),"\n",
                                  "p-value ",caseGroup1," vs ",refGroup," = ", as.character(format(round(Pval_Gr1_CTRL[i],3),nsmall=3))),
                     xlab="Phenotype Group", ylab= paste0("residuals_",as.character(var_i)))
      print(p)
      dev.off()
    }
    
  }
  model_df <- data.frame(var = variablesY, AvrgGr1 = AvrgGr1, AvrgCTRL = AvrgCTRL,
                         Pval_Gr1_CTRL = Pval_Gr1_CTRL,
                         MedianGr1 = MedianGr1, MedianCTRL = MedianCTRL)
  model_df$adj_Pval_Gr1_CTRL <- p.adjust(model_df$Pval_Gr1_CTRL, method="BH")
  write.csv(model_df,paste0(out_dir,"/",as.character(out_df_name),".csv"), row.names = FALSE)
}

