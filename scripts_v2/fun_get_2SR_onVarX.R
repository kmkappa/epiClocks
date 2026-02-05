##### FUNCTIONS : Residual 2SR model for 2 groups #####
#df <- df_ori; var_i <- "DNAmAge"; refGroup <- "R"; caseGroup1 = "NR"; VarX <- "Age"; variablesY <- dnam_ages; VarCov1 <- "Sex";VarCov2 <- NA;
df <- df; var_i <- "tnsc"; refGroup <- "Overuse0"; caseGroup1 = "Overuse1"; VarX <- "Age"; variablesY <- dnam_stem_cell_div; VarCov1 <- "Sex";VarCov2 <- NA;
df <- df; var_i <- "DNAmTL"; refGroup <- "Episodic"; caseGroup1 = "Chronic"; VarX <- "Age"; variablesY <- dnam_ages; VarCov1 <- "Sex";VarCov2 <- NA;

# model correcting for biological sex and complications
# update : preliminar step removing outliers - substituting them with NAs

get_2SR_onVarX <- function(df, VarX , VarCov1, VarCov2, variablesY, out_df_name, out_dir, refGroup, caseGroup1, age_lineplot, res_lineplot, boxplot){
  AvrgGr1 <- c(); AvrgCTRL <- c();
  MedianGr1 <- c(); MedianCTRL <- c();
  Pval_Gr1_CTRL <- c();
  num_tot <- c(); num_na_pre <- c(); num_na_post <- c();
  df_noOutliers <- df
  
  for(var_i in variablesY){
    print(var_i)
    if(sd(df[,c(var_i)])==0){
      num_na_pre <- c(num_na_pre, NA)
      num_na_post <- c(num_na_post, NA)
      num_tot <- c(num_tot, NA)
      Pval_Gr1_CTRL <- c(Pval_Gr1_CTRL,NA)
      AvrgGr1 <- c(AvrgGr1,NA)
      AvrgCTRL <- c(AvrgCTRL,NA)
      MedianGr1 <- c(MedianGr1,NA)
      MedianCTRL <- c(MedianCTRL,NA)
    }
    else{
      num_na_before_outliers <- sum(is.na(df[,c(var_i)]))

      # remove outliers
      if(refGroup=="Overuse0" | refGroup=="Episodic" ){
        quantiles <- quantile(df[(df$Group %in% c(refGroup,caseGroup1)),c(var_i)], probs = c(.25, .75), na.rm=T); 
        range <- 1.5 * IQR(df[(df$Group %in% c(refGroup,caseGroup1)),c(var_i)])
        df_noOutliers[(df_noOutliers$Group %in% c(refGroup,caseGroup1)),c(var_i)][!(df_noOutliers[(df_noOutliers$Group %in% c(refGroup, caseGroup1)),c(var_i)] > (quantiles[1] - range) & df_noOutliers[(df_noOutliers$Group %in% c(refGroup, caseGroup1)),c(var_i)] < (quantiles[2] + range))] <- NA
      }
      else{
        quantiles <- quantile(df[(df$Group %in% c(refGroup)),c(var_i)], probs = c(.25, .75), na.rm=T); 
        range <- 1.5 * IQR(df[(df$Group %in% c(refGroup)),c(var_i)])
        df_noOutliers[(df_noOutliers$Group %in% c(refGroup)),c(var_i)][!(df_noOutliers[(df_noOutliers$Group %in% c(refGroup)),c(var_i)] > (quantiles[1] - range) & df_noOutliers[(df_noOutliers$Group %in% c(refGroup)),c(var_i)] < (quantiles[2] + range))] <- NA
        
        quantiles <- quantile(df[(df$Group %in% c(caseGroup1)),c(var_i)], probs = c(.25, .75), na.rm=T); 
        range <- 1.5 * IQR(df[(df$Group %in% c(caseGroup1)),c(var_i)])
        df_noOutliers[(df_noOutliers$Group %in% c(caseGroup1)),c(var_i)][!(df_noOutliers[(df_noOutliers$Group %in% c(caseGroup1)),c(var_i)] > (quantiles[1] - range) & df_noOutliers[(df_noOutliers$Group %in% c(caseGroup1)),c(var_i)] < (quantiles[2] + range))] <- NA
      }
      
      num_na_with_outliers <- sum(is.na(df_noOutliers[,c(var_i)]))
      num_na_pre <- c(num_na_pre, num_na_before_outliers)
      num_na_post <- c(num_na_post, num_na_with_outliers)
      num_tot <- c(num_tot,nrow(df_noOutliers))

      # create linear model on ref group
      if(is.na(VarCov1) | is.na(VarCov2)){
        if(is.na(VarCov1) & is.na(VarCov2)){
          # model without covariates
          lmCTRL <- lm(eval(parse(text = var_i)) ~ eval(parse(text = VarX)), data = df_noOutliers[df_noOutliers$Group== refGroup,])
        }
        else if(is.na(VarCov1) & !(is.na(VarCov2))){
          # model only with VarCov2 covariates
          lmCTRL <- lm(eval(parse(text = var_i)) ~ eval(parse(text = VarX)) + eval(parse(text = VarCov2)), data = df_noOutliers[df_noOutliers$Group== refGroup,])
        }
        else if(is.na(VarCov2) & !(is.na(VarCov1))){
          # model only with VarCov1 covariates
          lmCTRL <- lm(eval(parse(text = var_i)) ~ eval(parse(text = VarX)) + eval(parse(text = VarCov1)), data = df_noOutliers[df_noOutliers$Group== refGroup,])
        }
      }
      else{
        # model with 2 covariates
        lmCTRL <- lm(eval(parse(text = var_i)) ~ eval(parse(text = VarX)) + eval(parse(text = VarCov1)) + eval(parse(text = VarCov2)), data = df_noOutliers[df_noOutliers$Group== refGroup,])
      }
      
      predCTRL <- predict(lmCTRL, df_noOutliers[df_noOutliers$Group==refGroup,])
      predGr1 <- predict(lmCTRL, df_noOutliers[df_noOutliers$Group==caseGroup1,])
      resCTRL <- df_noOutliers[df_noOutliers$Group==refGroup, var_i] - predCTRL
      resGr1 <- df_noOutliers[df_noOutliers$Group==caseGroup1, var_i] - predGr1
      df_noOutliers$residuals <- NA
      df_noOutliers$residuals[df_noOutliers$Group==refGroup] <- resCTRL
      df_noOutliers$residuals[df_noOutliers$Group==caseGroup1] <- resGr1
      
      out <- t.test(df_noOutliers$residuals[df_noOutliers$Group==caseGroup1],df_noOutliers$residuals[df_noOutliers$Group==refGroup])
      Pval_Gr1_CTRL <- c(Pval_Gr1_CTRL,out$p.value); AvrgGr1 <- c(AvrgGr1,out$estimate[1]); AvrgCTRL <- c(AvrgCTRL,out$estimate[2]);
      MedianGr1 <- c(MedianGr1,median(df_noOutliers[df_noOutliers$Group==caseGroup1, var_i],na.rm = T));
      MedianCTRL <- c(MedianCTRL,median(df_noOutliers[df_noOutliers$Group==refGroup, var_i],na.rm = T));
      
      # Visualizations
      # age_lineplot
      if(age_lineplot==T){
        if (file.exists(paste0(out_dir,"/age_lineplots"))){}else{dir.create(paste0(out_dir,"/age_lineplots"))}
        i <- which(variablesY == var_i)
        pdf(paste0(out_dir,"/age_lineplots/age_lineplot_",out_df_name, "_", as.character(var_i),".pdf"), width = 7, height = 6)
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
        if (file.exists(paste0(out_dir,"/res_lineplots"))){}else{dir.create(paste0(out_dir,"/res_lineplots"))}
        i <- which(variablesY == var_i)
        pdf(paste0(out_dir,"/res_lineplots/lineplot_",out_df_name, "_", as.character(var_i),".pdf"), width = 7, height = 6)
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
        if (file.exists(paste0(out_dir,"/boxplots"))){}else{dir.create(paste0(out_dir,"/boxplots"))}
        i <- which(variablesY == var_i)
        pdf(paste0(out_dir,"/boxplots/boxplot_",out_df_name, "_", as.character(var_i),".pdf"), width = 4, height = 6)
        p <- ggboxplot(df_noOutliers[(df_noOutliers$Group==refGroup | df_noOutliers$Group==caseGroup1),], y = 'residuals' , x= 'Group', add= 'jitter',color = 'Group',
                       order = c(as.character(caseGroup1), as.character(refGroup)),
                       main= paste0(as.character(var_i),"\n",
                                    "p-value ",caseGroup1," vs ",refGroup," = ", as.character(format(round(Pval_Gr1_CTRL[i],3),nsmall=3))),
                       xlab="Phenotype Group", ylab= paste0("residuals_",as.character(var_i))) + theme(plot.title = element_text(size=12))
        print(p)
        dev.off()
      }
    }
  }
  model_df <- data.frame(var = variablesY, Num_Tot= num_tot, Num_NA_pre = num_na_pre, Num_NA_post = num_na_post,
                         AvrgGr1 = AvrgGr1, AvrgCTRL = AvrgCTRL,
                         MedianGr1 = MedianGr1, MedianCTRL = MedianCTRL,
                         Pval_Gr1_CTRL = Pval_Gr1_CTRL)
  for(m in p.adjust.methods[-length(p.adjust.methods)]){
    model_df <- cbind(model_df, p.adjust(model_df$Pval_Gr1_CTRL, method=m))
    colnames(model_df)[ncol(model_df)] <- paste0("adj_Pval_Gr1_CTRL_",m)
  }
  #p.adjust.methods=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
  
  write.csv(model_df,paste0(out_dir,"/",as.character(out_df_name),".csv"), row.names = FALSE)
}

