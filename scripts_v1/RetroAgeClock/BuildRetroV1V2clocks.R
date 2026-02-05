
library(data.table)
library(caret)
library(glmnet)

# Set file paths to DNA methylation dataset for RetroelementV1 CpGs from 12,670 people age 12-100 years old
file_path <- 'RetroAgeV1MethylationDataset.csv'
metadata_path <- 'RetroMetadata.csv'


# Read the CSV files using fread
data <- fread(file_path, header = TRUE, data.table = FALSE)
metadata <- read.csv(metadata_path, header = TRUE)

# Set row names to first column and remove it
rownames(data) <- data[, 1]
data <- data[, -1]

# Transpose DNA methylation dataset
methylation_data <- t(data)


# Extract chronological ages from metadata
ages <- metadata$Decimal.Chronological.Age


# Create a data partition
set.seed(42) # Set seed for reproducibility
train_indices <- createDataPartition(ages, p = 0.8, list = FALSE)

# Split the DNA methylation retroelement dataset into training and validation sets
beta_noob <- data
beta_noob_train <- beta_noob[, train_indices]
beta_noob_validation <- beta_noob[, -train_indices]

train_targets <- ages[train_indices]
validation_targets <- ages[-train_indices]

datMeth_train <- t(beta_noob_train)
datMeth_validation <- t(beta_noob_validation)

#Build Retroelement Clock utilizing glmnet 
library(doParallel)
registerDoParallel(8)

cv = cv.glmnet(datMeth_train, train_targets, nfolds=10,alpha=0.5, family="gaussian", parallel=TRUE)
fit = glmnet(datMeth_train, train_targets, family="gaussian", alpha=0.5, nlambda=100)
plot(cv)


#Plot the training and test chronological age and retroelement predictions 
plot(train_targets,predict(fit,datMeth_train,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(train_targets,predict(fit,datMeth_train,s = cv$lambda.min))

plot(validation_targets,predict(fit,datMeth_validation,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(validation_targets,predict(fit,datMeth_validation,s = cv$lambda.min))


#Save coefficents
tmp_coeffs <- coef(cv, s = "lambda.min")
predictors <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)


write.csv(predictors, file="Retroelement_AgeV1_coefficents.csv", quote=F)
sessionInfo()














library(data.table)
library(caret)
library(glmnet)

# Set file paths to DNA methylation dataset for RetroelementV2 CpGs from 12,670 people age 12-100 years old
file_path <- 'RetroAgeV2MethylationDataset.csv'
metadata_path <- 'RetroMetadata.csv'


# Read the CSV files using fread
data <- fread(file_path, header = TRUE, data.table = FALSE)
metadata <- read.csv(metadata_path, header = TRUE)

# Set row names to first column and remove it
rownames(data) <- data[, 1]
data <- data[, -1]

# Transpose DNA methylation dataset
methylation_data <- t(data)


# Extract chronological ages from metadata
ages <- metadata$Decimal.Chronological.Age


# Create a data partition
set.seed(42) # Set seed for reproducibility
train_indices <- createDataPartition(ages, p = 0.8, list = FALSE)

# Split the DNA methylation retroelement dataset into training and validation sets
beta_noob <- data
beta_noob_train <- beta_noob[, train_indices]
beta_noob_validation <- beta_noob[, -train_indices]

train_targets <- ages[train_indices]
validation_targets <- ages[-train_indices]

datMeth_train <- t(beta_noob_train)
datMeth_validation <- t(beta_noob_validation)

#Build Retroelement Clock utilizing glmnet 
library(doParallel)
registerDoParallel(8)

cv = cv.glmnet(datMeth_train, train_targets, nfolds=10,alpha=0.5, family="gaussian", parallel=TRUE)
fit = glmnet(datMeth_train, train_targets, family="gaussian", alpha=0.5, nlambda=100)
plot(cv)


#Plot the training and test chronological age and retroelement predictions 
plot(train_targets,predict(fit,datMeth_train,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(train_targets,predict(fit,datMeth_train,s = cv$lambda.min))

plot(validation_targets,predict(fit,datMeth_validation,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(validation_targets,predict(fit,datMeth_validation,s = cv$lambda.min))


#Save coefficents
tmp_coeffs <- coef(cv, s = "lambda.min")
predictors <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)


write.csv(predictors, file="Retroelement_AgeV2_coefficents.csv", quote=F)
sessionInfo()

