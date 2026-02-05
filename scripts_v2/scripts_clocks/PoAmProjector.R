#' Generate the Dunedin Methylation Pace of Aging Scores!
#'
#' \code{PoAmProjector} returns the Dunedin Pace of Aging Methylation Scores
#'
#' @param betas A numeric matrix containing the percent-methylation for each probe.  Missing data should be 'NA's.  The rows should be probes, with the probe ID as the row name, and the columns should be samples, with sample names as the column name.
#' @param proportionOfProbesRequired (default: 0.8).  This value specificies the threshold for missing data (see description for more details on how missing data is handled)
#' @return A list of mPoA values.  There will be one element in the list for each mPoA model.  Each element will consist of a numeric vector with mPoA values.  The names of the values in the vector will be the sample names from the 'betas' matrix.
#' @details This function returns the Dunedin Methylation Pace of Aging scores for methylation data generated from either the Illumina 450K array or the Illumina EPIC array.  The Age38 score is the one described in the eLife paper (2020).  The Age45 score is one that has been trained on data based on 3 waves of collection (26, 38, and 45).  The manuscript is currently in preparation, but has been shown to be more accurate than the Age38 score.
#' Missing data handled in two different ways (and the threshold for both is set by the 'proportionOfProbesRequired' parameter).  First, if a sample is missing data for more probes than the threshold, the sample will get an NA back for a score.  If a particular probe is missing fewer samples than the threshold, then missing data is set to the mean in the provided 'betas' matrix.  If a probe is missing more samples than the threshold, then all samples in the 'betas' matrix have their value replaced with the mean of the training data for that particular model.
#' Because of how we handle missing data, it is reccomended that entire cohorts be run at once as a large 'betas' matrix.
#' @examples
#' PoAmProjector(betas)


rm(list=ls())
library(data.table)
source("/home/PERSONALE/katarzyn.kwiatkowsk2/DunedinPACE/PoAmProjector.R") # load function calculating score
load(file="/home/PERSONALE/katarzyn.kwiatkowsk2/DunedinPACE/sysdata.rda") # loads model

# Load my data
out_dir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE"
experimental_name <- "epimode"

load(file=paste0(out_dir,"/out_script0/",experimental_name,"_bVals_detectionPvalFlt.RData"))
data <- bVals
head(data)[1:3,1:3]
dim(data)

betas <-as.matrix(data[,-1])
rownames(betas) <- data$ID_REF

out <- PoAmProjector(betas,proportionOfProbesRequired=0.8)$DunedinPoAm
out <- data.frame(out)
colnames(out)[1] <- "PoAmProjector"
save(out,file=paste0(out_dir,"/out_clocks/",experimental_name,"_PoAmProjectorScore.RData",sep=""))


##### NEW DunedinPACE score #####

# REF : https://github.com/danbelsky/DunedinPACE
# R in conda env cart2

rm(list=ls())
library(data.table)

load(file="/home/PERSONALE/katarzyn.kwiatkowsk2/DunedinPACE/sysdata.rda") # loads model

# Load my data
out_dir <- "/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE"
experimental_name <- "epimode"

load(file=paste0(out_dir,"/out_script0/",experimental_name,"_bVals_detectionPvalFlt.RData"))
data <- bVals
head(data)[1:3,1:3]
dim(data)

betas <-as.matrix(data[,-1])
rownames(betas) <- data$ID_REF

library("DunedinPACE")
mPACE <- PACEProjector(betas)
out <- as.data.frame(mPACE)

save(out, file=paste0(out_dir,"/out_clocks/",experimental_name,"_DunedinPaceScore.RData"))
