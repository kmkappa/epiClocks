
rm(list=ls())
################################
library(minfi)
library(limma)
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")

out_dir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/"
ifelse(!dir.exists(file.path(out_dir)),
        dir.create(file.path(out_dir)),
        "Directory Exists")
setwd(out_dir)
#scp -J katarzyn.kwiatkowsk2@137.204.50.15 /Users/k/data/Epimode katarzyn.kwiatkowsk2@137.204.48.210:/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/

# DIR with sampleSheet with samples to analyze
ssDir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE"
# DIR with all samples to analyze
baseDir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/idat"
a <- list.files(baseDir, recursive = TRUE, full.names=T)
a
experiment_label <- "epimode"
targets <- read.metharray.sheet(ssDir, patter = "SampleSheetEpimo")
targets
dim(targets) # 142 15
targets[1:3,]

RGset <- read.metharray.exp(targets = targets, extended=T, force = T)
RGset

sampleNames(RGset) <- targets$SampleID
#save(RGset,file= paste0(experiment_label,"_RGset.RData"))

# Calculate detection p-values
detP <- detectionP(RGset)
head(detP)
# remove poor quality samples # here: 0 removed, all kept
keep <- colMeans(detP) < 0.05
table(keep) # 142 TRUE
RGset <- RGset[,keep]
RGset
targets <- targets[keep,]
# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP) # 866238     142

#####################################
### Try different Normalizations
####################################

# ## Funnorm Normalization
# mSetSq <- preprocessFunnorm(RGset) 
# # visualise what the data looks like before and after normalization
# pdf(paste(experiment_label,"_densityPlots_beforeAndAfterNormalization_funnorm.pdf",sep=""))
# par(mfrow=c(1,2))
# densityPlot(RGset, sampGroups=targets$Group,main="Raw", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# densityPlot(getBeta(mSetSq), sampGroups=targets$Group,
#             main="Normalized", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# dev.off()

# ## SWAN Normalization
# mSetSq <- preprocessSWAN(RGset) 
# # visualise what the data looks like before and after normalization
# pdf(paste(experiment_label,"_densityPlots_beforeAndAfterNormalization_swan.pdf",sep=""))
# par(mfrow=c(1,2))
# densityPlot(RGset, sampGroups=targets$Group,main="Raw", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# densityPlot(getBeta(mSetSq), sampGroups=targets$Group,
#             main="Normalized", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# dev.off()

# ## Noob Normalization
# mSetSq <- preprocessNoob(RGset) 
# # visualise what the data looks like before and after normalization
# pdf(paste(experiment_label,"_densityPlots_beforeAndAfterNormalization_noob.pdf",sep=""))
# par(mfrow=c(1,2))
# densityPlot(RGset, sampGroups=targets$Group,main="Raw", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# densityPlot(getBeta(mSetSq), sampGroups=targets$Group,
#             main="Normalized", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# dev.off()

# ## Illumina Normalization
# mSetSq <- preprocessIllumina(RGset) 
# # visualise what the data looks like before and after normalization
# pdf(paste(experiment_label,"_densityPlots_beforeAndAfterNormalization_illumina.pdf",sep=""))
# par(mfrow=c(1,2))
# densityPlot(RGset, sampGroups=targets$Group,main="Raw", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# densityPlot(getBeta(mSetSq), sampGroups=targets$Group,
#             main="Normalized", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# dev.off()

# ## Quantile Normalization
# mSetSq <- preprocessQuantile(RGset) 
# # visualise what the data looks like before and after normalization
# pdf(paste(experiment_label,"_densityPlots_beforeAndAfterNormalization_quantile.pdf",sep=""))
# par(mfrow=c(1,2))
# densityPlot(RGset, sampGroups=targets$Group,main="Raw", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# densityPlot(getBeta(mSetSq), sampGroups=targets$Group,
#             main="Normalized", legend=FALSE)
# #legend("top", legend = levels(factor(targets$Group)), 
# #       text.col=brewer.pal(8,"Dark2"))
# dev.off()

# scp -J katarzyn.kwiatkowsk2@137.204.50.15 katarzyn.kwiatkowsk2@137.204.48.210:/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/epimode_densityPlots_beforeAndAfterNormalization* /Users/k/data/Epimode/out_script0/

#####################################
### The Winner Normalization
####################################

## Quantile Normalization
mSetSq <- preprocessQuantile(RGset) 
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])
dim(bVals) # 856199     142
bVals <- data.frame(rownames(bVals), bVals)
colnames(bVals)[1] <- "ID_REF"
# Save normalized b-values to use for epigenetic clocks
write.csv(bVals,file =paste0(experiment_label,"_bVals_detectionPvalFlt.csv"), row.names = F)
save(bVals,file =paste0(experiment_label,"_bVals_detectionPvalFlt.RData"))

# # Download to local pc
# scp -r -J katarzyn.kwiatkowsk2@137.204.50.15 katarzyn.kwiatkowsk2@137.204.48.210:/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/epimode_bVals_Raw.RData /Users/k/data/Epimode/DNAmClocks/out_script0/
# scp -r -J katarzyn.kwiatkowsk2@137.204.50.15 katarzyn.kwiatkowsk2@137.204.48.210:/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/epimode_bVals_Raw_detectionPvalFlt.RData /Users/k/data/Epimode/DNAmClocks/out_script0/
# scp -r -J katarzyn.kwiatkowsk2@137.204.50.15 katarzyn.kwiatkowsk2@137.204.48.210:/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/epimode_bVals_detectionPvalFlt.RData /Users/k/data/Epimode/DNAmClocks/out_script0/


