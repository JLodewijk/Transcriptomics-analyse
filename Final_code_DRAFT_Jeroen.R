#########################################
# Name: Raquel Manzano & Jeroen Lodewijk#
# Student numbers: & 930203525030       #
# Project: Bioinformatics               #
#########################################

#############
# Packages #
#############

# install.packages("car")
#
# install.packages("class")
#
# install.packages("leaps")
#
# install.packages("MASS")

####################
#### Libraries #####
####################

# Library used for doing a variance inflation factors test.
library(car)

# Library used for performing k-Nearest Neighbour Classification.
library(class)

# Library used for regression subset selection
library(leaps)

# Library used for predictions.
library(MASS)

# Library used for the plotting of the PCA.
library(plotrix)

# Library used for performing a randdomForest.
library(randomForest)

# Library is used for the creation of trees.
library(tree)

# Used for performing a svm.
library(e1071)

# Used for the glmnet for lasso and ridge linear regression
library(gbm)

# Lasso adn the ridge methods
library(glmnet)

####################
#### Functions #####
####################

#' GlmPredictionErrorRate
#' Perform a prediction on the test data and calculate the error rate of that prediction.
#' @param glm.fit: Fitted generalized linear model that has been trained by a trainings set.
#' @param tissue.1: First tissue you want to use for prediction.
#' @param tissue.2: Second tissue you want to use for prediction.
#' @param test.data: A test dataset, that doesn't contain any trainings data.
#' @param type.prediction: Type of prediction that is beinf performed by the predict(). Default on response.
#'
#' @return Table of predictions and the error rate of the predictions against the test data.
GlmPredictionErrorRate <-
  function(glm.fit,
           tissue.1,
           tissue.2,
           test.data,
           type.prediction = "response") {
    # Make a predicion of the glm on the test data.
    glm.probs <-
      predict(glm.fit, test.data, type = type.prediction)
    
    # Create replicates of tissue.1 based on the number of observation the test data has.
    glm.pred <- rep(tissue.1, nrow(test.data))
    
    # Any predicitions higher than 0.5 will be set as tissue.2
    glm.pred[glm.probs > 0.5] <- tissue.2
    
    # Give the results of the prediction vs the test data.
    print(table(glm.pred, test.data$tissue))
    return(1 - mean(glm.pred == test.data$tissue))
  }

#' QdaLdaPredictionError
#' Gives the prediction made by a lda or qda on the 
#' @param model: qda or lda fit thas has been performed.
#' @param test.data: A test dataset, that doesn't contain any trainings data.
#'
#' @return Table of predictions and the error rate of the predictions against the test data.
QdaLdaPredictionError <- function(model, test.data){
  class <- predict(model, test.set)$class
  print(table(class, test.set$tissue))
  return(1 - mean(class == test.set$tissue))
}


#' PredictivePerformanceLm
#' Check how well the predictive performance is of the linear model.
#' @param y: the y used in the creation of the linear model.
#' @param data.set: the dataset that has been used for the creation of the linear model.
#' @param training.data: the trainings data indexes.
#' @param lm.training: the linear model, where the trainings data has been apllied to.
#'
#' @return R-squeare and fractions of variability explained by the model.
PredictivePerformanceLm <-
  function(y, data.set, training.data, lm.training) {
    # Create a test dataset to check how much of the variability is explained by the current model.
    # Take all the entries that do not belong to the trainings dataset.
    test <- -training.data
    
    # Create a prediction using the lm of the trainings data on the test data.
    test.pred <- predict(lm.training, newdata = data.set[test,])
    
    # See how well it predicts the y.
    test.y    <- data.set[test, y]
    
    # Calculate the Sum of Squares total, residual, regression and total. In order to know the fraction of variability explained by the model.
    SS.total      <- sum((test.y - mean(test.y)) ^ 2)
    SS.residual   <- sum((test.y - test.pred) ^ 2)
    SS.regression <- sum((test.pred - mean(test.y)) ^ 2)
    
    # Calculate the R-square (NOT the fraction of variability explained by the model)
    cat("R-square = ", 1 - SS.residual / SS.total)
    cat("\nFraction of variability explained by the model = ",
        SS.regression / SS.total)
  }

#' train.selection
#' Creation of a train variable, chan choose on within a range of a percentage.
#' @param percentage: percentage of the data that is going to be used for the training of the data
#' @param data: data itself where the training selection will be taken from.
#' @param seed: A fixed seed, that is taken to make the sample taken reproducable.
#' @param random.seed: determines if a random seed is even set.
#'
#' @return train: vector of indexes that contain the training selection.
train.selection <-
  function(percentage,
           data,
           random.seed = TRUE,
           seed = as.numeric(1)) {
    if (random.seed == TRUE) {
      # Set a fixed seed to make the data reproducable.
      set.seed(seed)
    }
    
    # Generate a vector of indexes based on the sample and the seed.
    train <- sample(1:nrow(data), round(nrow(data) * percentage))
    return (train)
  }

#' error.rate.prediction.trees
#' Predict the error rate of a tree.
#' @param tree.data: the tree that has been build.
#' @param dataset: dataset that has been used in the creation of the tree.
#' @param test.set: a list of indexes indicating the postions of the test dataset.
#' @param type.prediction: what kind of prediction has to be made in the predict().
#'
#' @return error rate of the prediction compared to the test set.
error.rate.prediction.trees <-
  function(tree.data,
           dataset,
           test.set,
           type.prediction = "class") {
    # Make a prediction using the tree data on the test data set.
    prediction.tree <-
      predict(tree.data, newdata = dataset[test.set, ], type = type.prediction)
    
    # Get the y data out of the dataset, create a list. Otherwise you get errors in the return statement.
    test.data <- dataset[test.set, ]$tissue
    
    # Error test of the prediction against the test set.
    return(1 - mean(prediction.tree == test.data))
  }

#' tree.pruner
#' Prune trees by using cross-validation and selection for the best set usign the min.set.selection().
#' @param tree.data: tree data that is going to be pruned.
#' @param seed: randomseed, will not be used if randomness is TRUE.
#' @param randomness: if a randomseed will be used or randomness will be allowed.
#' @param regression.tree: boolean that indicates if the tree.data is a boolean, default is on FALSE.
#' @return pruned tree
tree.pruner <-
  function(tree.data,
           seed = 1,
           randomness = FALSE,
           regression.tree = FALSE) {
    # Set a seed if randomness is not wanted.
    if (randomness == FALSE) {
      set.seed(seed)
    }
    
    # Perform cross-validation on the tree.
    if (regression.tree == FALSE) {
      # CV for the classification problem.
      cv.tissue <- cv.tree(tree.data, FUN = prune.misclass)
    }
    else{
      # CV for the regression problem.
      cv.tissue <- cv.tree(tree.data, FUN = prune.tree)
    }
    # Get the best.set using the min.set.selection.
    best.set <- min.set.selection(cv.tissue)
    
    # Prune the tree and return it.
    if (regression.tree == FALSE) {
      # Prunning method for the classification problem.
      return(prune.misclass(tree.data, best = best.set))
    }
    else{
      # Prunning method for the regression problem.
      return(prune.tree(tree.regression, best = best.set))
    }
  }

#' Cols
#' Creates a range of colors that can be used for plotting.
#' @param vec: is a vector that will be given colors ID's.
#'
#' @return Color ID's for the use in plotting
Cols <- function(vec) {
  cols <- rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

#' min.set.selection
#' Selection of the best size for trees and cross validation.
#' @param cvset: cross validated set.
#'
#' @return Minima set selection
min.set.selection <- function(cvset) {
  #rename
  min.dev <- which.min(cvset$dev)
  return (cvset$size[min.dev])
}

#' tissue.selection
#' Takes out two tissues
#' @param tissue1: first tissue name as a string.
#' @param tissue2: second tissue name as a string.
#' @param data: data as a data.frame.
#'
#' @return mydata: two tissue dataset as a data.frame
tissue.selection <-
  function(tissue1, tissue2, data = data.frame(expr4T.filtered)) {
    mydata <-
      data[which(data$tissue == tissue1 | data$tissue == tissue2), ]
    mydata <- droplevels(mydata)
    print(dim(mydata))
    print(table(mydata$tissue))
    return (mydata)
  }

#' is.finite.data.frame
#' Check if a dataframe does only contain finite data.
#' @param obj: dataframe that is being checked.
#'
#' @return True or False if a column contains finite data only.
is.finite.data.frame <- function(obj) {
  sapply(
    obj,
    FUN = function(x)
      all(is.finite(x))
  )
}

#' DetermineEffectPsRandomForest
#' Using different P's on the dataset, to see which one has the smallest mse.
#' @param dataset: dataset used for performing a randomForest upon.
#' @param train.set: list of indexes indicating the entries of the trainings set.
#' @param test.set: list of indexes indicating the entries of the test set.
#' @param column.index: column that is being used for the training and testing of data.
#' @return data.frame containing each of the P's MSE as a column.
DetermineEffectPsRandomForest <-
  function(dataset,
           train.set,
           test.set,
           column.index) {
    # Use different P's to see the effect on randomForest.
    p <- ncol(dataset) - 1
    p.2 <- p / 2
    p.3 <- sqrt(p)
    
    # Create different train and test sets for the data for the classification problem
    train.X <- dataset[train.set, -column.index]
    test.X <- dataset[test.set, -column.index]
    train.Y <- dataset[train.set, column.index]
    test.Y <- dataset[test.set, column.index]
    
    # Run the randomforest, each of them has a different P. But all of them have the same ntree and datasets.
    rf.tissue.p <- randomForest(
      train.X,
      train.Y,
      xtest = test.X,
      ytest = test.Y,
      mtry = p,
      ntree = 500
    )
    
    rf.tissue.p.2 <- randomForest(
      train.X,
      train.Y,
      xtest = test.X,
      ytest = test.Y,
      mtry = p.2,
      ntree = 500
    )
    
    rf.tissue.p.3 <- randomForest(
      train.X,
      train.Y,
      xtest = test.X,
      ytest = test.Y,
      mtry = p.3,
      ntree = 500
    )
    
    # Return a data.frame containg the mse of each of the randomForests
    return(
      data.frame(
        p1.mse = rf.tissue.p$test$mse,
        p2.mse = rf.tissue.p.2$test$mse,
        p3.mse = rf.tissue.p.3$test$mse
      )
    )
  }


#' steps for a complete glm prediction
#' @param tissues.pair: a vector of two tissues
#' @param data.set: used data
#' @param fit.type: glm fit as default, if change to lda problems in the alst part (change?)
#'
#' @return result: test error rate
glm.model <- function(fit.type = glm, data.set, tissues.pair) {
  set.seed(1)
  #pairs of tissues
  tissue1 <- tissues.pair[1]
  tissue2 <- tissues.pair[2]
  
  #specific data set
  mydat <- tissue.selection(tissue1, tissue2, data.set)
  
  train <- sample(1:nrow(mydat), round(nrow(mydat) * 0.35))
  test.set <- mydat[-train, ]
  
  #model
  tissues <- mydat$tissue
  new.data <- data.frame(mydat[, 1:10], tissues)
  fit <- fit.type(tissues ~ .,
                  data = new.data,
                  family = binomial,
                  subset = train)
  
  #predictions
  
  probs <-  predict(fit, test.set, type = "response")
  pred <- rep(tissue1, nrow(test.set))
  pred[probs > 0.5] = tissue2
  cat(table(pred, test.set$tissue))
  
  
  return (1 - mean(pred == test.set$tissue))
}

#############
# Main code #
#############

#### WEEK 1 ####
#### Pre-processing ####

# Load in the dataset.
# expr4T <- read.table("M:/Pattern Recognition/Project/expr4T.dat", sep="")
expr4T <-
  read.table("F:/Dropbox/Pattern Recognition/Week 1/Friday/expr4T.dat",
             sep = "")
###
dim(expr4T) #Checking dimentions

####REDUCTION OF THE DATA SET####

nn <- ncol(expr4T) - 1 #Last column tissue label
m <- sapply(expr4T[, 1:nn], mean)
s <- sapply(expr4T[, 1:nn], sd)
sm <- s / m

# Filter genes based on minimum mean and stdev/mean
minsm <- mean(sm)
minm <- mean(m)

m <- c(m, minm + 1)
sm <- c(sm, minsm + 1)

expr4T.filtered <- expr4T[, which(sm > minsm & m > minm)]
dim(expr4T.filtered)

#Remove variables that we are not going to use more
rm(m, minm, minsm, nn, s, sm)

# Check if these genes are in the dataset.
"ENSG00000271043.1_MTRNR2L2" %in% colnames(expr4T.filtered)
"ENSG00000229344.1_RP5.857K21.7" %in% colnames(expr4T.filtered)

####ANALYSIS OF THE GENE EXPRESSION LEVELS RELATIONSHIP####

#Creation of two tissue data set

tissue1 <- "brain_putamen"
tissue2 <- "brain_cortex"

mydat <- tissue.selection(tissue1, tissue2, data = expr4T.filtered)

#LINEAR REGReSSION MODELS

set.seed(1)
####Training set selection####
train.expr4T.data <-
  sample(1:nrow(expr4T.filtered), round(nrow(expr4T.filtered) * 0.35))

#Train set for two tissue data separated
train.mydat <- sample(1:nrow(mydat), round(nrow(mydat) * 0.35))

#### Performance of a linear model using all the data expected ENSG00000271043.1_MTRNR2L2 gene ####

lm.fit.1 <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = expr4T.filtered,
     subset = train.expr4T.data)

lm.fit.1.two.tissues <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = mydat[, 1:70],
     subset = train.mydat)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2
plot(lm.fit.1)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2 fpr two tissue data
plot(lm.fit.1.two.tissues)

#### Approach for avoiding patterns in linearity plot(NOT SUCCESS) Remove? It is useless so maybe we can just omit it. #####
reduced <-
  c(names(lm.fit.1$fitted.values[which(lm.fit.1$fitted.values > 500)]))
a <- expr4T.filtered[reduced, ]

lm.fit.reduced <-
  lm(ENSG00000225972.1_MTND1P23 ~ .,
     data = a,
     subset = train.expr4T.data)
plot(lm.fit.reduced)

rm(reduced, a, lm.fit.reduced)

#### Checking predictions ####

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = expr4T.filtered,
  training.data = train.expr4T.data,
  lm.training = lm.fit.1
)

# Check if there is collinearity, if true collinearity is present.
sqrt(vif(lm.fit.1)) > 2

# Data shows collinearity. Explains why the model is highly significant, while few individual predictors are significant.
# Therefore the results of the PredictivePerformanceLm are not trustworty

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = mydat,
  training.data = train.mydat,
  lm.training = lm.fit.1.two.tissues
)
# Negative R-square shows that there is dependance of y on x is non-linear. Therefore  negative R-square.

#### Second gene of interest ENSG00000229344.1_RP5.857K21.7 ####

# Perform a linear model using all the terms in te model
lm.fit.2 <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = expr4T.filtered,
     subset = train.expr4T.data)

lm.fit.2.two.tissues <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = mydat[, 1:70],
     subset = train.mydat)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2
plot(lm.fit.2)
# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2 in two tissue data set
plot(lm.fit.2.two.tissues)

PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = expr4T.filtered,
  training.data = train.expr4T.data,
  lm.training = lm.fit.2
)
# Seems the trainings set explains quite a lot of the test set. So it seems like a good model.

PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = mydat,
  training.data = train.mydat,
  lm.training = lm.fit.2.two.tissues
)

sqrt(vif(lm.fit.2)) > 2
#Again strong collinearity, results are not trustworthy.

####  OTHER TWO GENES ####
## ENSG00000198695.2_MT.ND6 was selected for being the one with the highest mean
## ENSG00000125144.9_MT1G lowest mean

# Perform a linear model using all the terms in te model
lm.fit.3 <-
  lm(ENSG00000198695.2_MT.ND6 ~ .,
     data = expr4T.filtered,
     subset = train.expr4T.data)

lm.fit.3.two.tissues <-
  lm(ENSG00000198695.2_MT.ND6 ~ .,
     data = mydat,
     #heRE WE USE THE WHOLE DATA SET, OTHERWISE THE GENE IS NOT FOUND
     subset = train.mydat)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2
plot(lm.fit.3)
# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2 in two tissue data set
plot(lm.fit.3.two.tissues)

PredictivePerformanceLm(
  y = "ENSG00000198695.2_MT.ND6",
  data.set = expr4T.filtered,
  training.data = train.expr4T.data,
  lm.training = lm.fit.3
)
# Seems the trainings set explains quite a lot of the test set. So it seems like a good model.

PredictivePerformanceLm(
  y = "ENSG00000198695.2_MT.ND6",
  data.set = mydat,
  training.data = train.mydat,
  lm.training = lm.fit.3.two.tissues
)

sqrt(vif(lm.fit.3)) > 2

#### SUMMARIES OF ALL LINEAR MODELS (Examples of coeff below in case we want to use it)####
summary(lm.fit.1) #14 SIGNIFICANT COEFF
#ENSG00000225972.1_MTND1P23                4.203e-02  1.210e-02   3.474 0.000595
#ENSG00000064787.8_BCAS1                   1.277e+00  8.954e-01   1.426 0.154862

summary(lm.fit.1.two.tissues) # 4
#ENSG00000225972.1_MTND1P23                  0.05533    0.03846   1.439   0.2237
#ENSG00000104419.10_NDRG1                   11.57953    3.43752   3.369   0.0281 *

summary(lm.fit.2) # 5
#ENSG00000125144.9_MT1G                   -5.575e+00  1.840e+00  -3.029 0.002681 **
#ENSG00000173267.9_SNCG                    6.329e+00  1.899e+00   3.332 0.000977 ***

summary(lm.fit.2.two.tissues) # 6
#ENSG00000234745.5_HLA.B                     58.3789    20.2977   2.876   0.0452 *
#ENSG00000237973.1_hsa.mir.6723               0.5232     0.1702   3.075   0.0371 *

summary(lm.fit.3) # 11
#tissuebrain_anteriorcortex                6.144e+01  1.379e+01   4.456 1.21e-05 ***
#ENSG00000125148.6_MT2A                    1.746e-01  4.950e-02   3.527 0.000491 ***

summary(lm.fit.3.two.tissues) # 0

###CONVERTION OF THE DATA (Log) ####
log.mydat <- log(mydat[-ncol(mydat)])
log.expr4T <- log(expr4T.filtered[-ncol(expr4T.filtered)])



# Check if there are any infs in the log-transformed datasets.
is.finite.data.frame(log.mydat)
is.finite.data.frame(log.expr4T)

# Works only on matrixes
log.mydat <- as.matrix(log.mydat)

# Replace any non finites with 0, this is needed to get the logs to work.
log.mydat[!is.finite(log.mydat)] <- 0
log.expr4T <- as.matrix(log.expr4T)
log.expr4T[!is.finite(log.expr4T)] <- 0

# Converted back to a data.frame for the lm's.
log.mydat <- as.data.frame(log.mydat)
log.expr4T <- as.data.frame(log.expr4T)

#### Perform a linear model on the log-transformed data. ####
lm.fit.1.log <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = log.expr4T,
     subset = train.expr4T.data)

plot(lm.fit.1.log)

#Predicted values
test.pred <-
  predict(lm.fit.1.log, newdata = log.expr4T[-train.expr4T.data,])
#Plot predicted against observed
plot(
  log.expr4T$ENSG00000271043.1_MTRNR2L2[-train.expr4T.data],
  test.pred,
  xlab = "Predicted values for ENSG00000271043.1_MTRNR2L2",
  ylab = "Observed values for ENSG00000271043.1_MTRNR2L2"
)
abline(0, 1, col = "red")

log.mydat <-
  log.mydat[, 1:70] #again we select a small data set, otherwise it doesn't work

#Two tissue model
lm.fit.1.two.tissues.log <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = log.mydat,
     #again we select a small data set, otherwise it doesn't work
     subset = train.mydat)

plot(lm.fit.1.two.tissues.log)

#Predicted values
test.pred <-
  predict(lm.fit.1.two.tissues.log, newdata = log.mydat[-train.mydat,])

#Plot predicted against observed
plot(
  log.mydat$ENSG00000271043.1_MTRNR2L2[-train.mydat],
  test.pred,
  xlab = "Predicted values for ENSG00000271043.1_MTRNR2L2",
  ylab = "Observed values for ENSG00000271043.1_MTRNR2L2"
)
abline(0, 1, col = "red")

#### Checking predictive performances ####

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = log.expr4T,
  training.data = train.expr4T.data,
  lm.training = lm.fit.1.log
)

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = log.mydat[, 1:50],
  training.data = train.mydat,
  lm.training = lm.fit.1.two.tissues.log
)
#We can see that the R squared for two tissue data is worst.

#### Gene ENSG00000229344.1_RP5.857K21.7 model ####
lm.fit.2.log <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = log.expr4T,
     subset = train.expr4T.data)
lm.fit.2.two.tissues.log <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = log.mydat,
     subset = train.mydat)

plot(lm.fit.2.log)
plot(lm.fit.2.two.tissues.log)

test.pred <-
  predict(lm.fit.2.log, newdata = log.expr4T[-train.expr4T.data,])

plot(
  log.expr4T$ENSG00000229344.1_RP5.857K21.7[-train.expr4T.data],
  test.pred,
  xlab = "Predicted values for ENSG00000229344.1_RP5.857K21.7",
  ylab = "Observed values for ENSG00000229344.1_RP5.857K21.7"
)
abline(0, 1, col = "red")

test.pred <-
  predict(lm.fit.2.two.tissues.log, newdata = log.mydat[-train.mydat,])

plot(
  log.mydat$ENSG00000229344.1_RP5.857K21.7[-train.mydat],
  test.pred,
  xlab = "Predicted values for ENSG00000229344.1_RP5.857K21.7",
  ylab = "Observed values for ENSG00000229344.1_RP5.857K21.7"
)
abline(0, 1, col = "red")

#### Checking predictive performances ####
PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = log.expr4T,
  training.data = train.expr4T.data,
  lm.training = lm.fit.2.log
)

PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = log.mydat,
  training.data = train.mydat,
  lm.training = lm.fit.2.two.tissues.log
)

#### COMPARATION OF PERFORMANCES BT TISSUES ####

tissues.pairs <- c(
  "brain_hippocampus" ,
  "brain_nucleusaccumbens",
  "brain_spinalcord" ,
  "brain_substantianigra",
  "brain_cerebellarhemisphere" ,
  "brain_cerebellum",
  "brain_cerebellum" ,
  "brain_amygdala",
  "brain_cortex" ,
  "brain_putamen"
)
set.seed(1)
#### FIRST PAIR: "brain_hippocampus" , "brain_nucleusaccumbens" ####
#pairs of tissues
tissue1 <- "brain_cerebellarhemisphere"
tissue2 <- "brain_cerebellum"

#specific data set
set.seed(1)
mydat <- tissue.selection(tissue1, tissue2, expr4T.filtered)

train <- sample(1:nrow(mydat), round(nrow(mydat) * 0.35))
test.set <- mydat[-train, ]

#GLM model####

glm.fit <- glm(tissue ~ .,
               data = mydat,
               family = binomial,
               subset = train)

#predictions
GlmPredictionErrorRate(glm.fit = glm.fit, tissue.1 = tissue1, tissue.2 = tissue2,test.data = test.set)

#LDA model####
lda.fit <- lda(tissue ~ .,
               data = mydat,
               subset = train)
QdaLdaPredictionError(model = lda.fit, test.data = test.set)

#QDA model ####
qda.fit <- qda(tissue ~ ., # DOES NOT WORK
               data = mydat,
               subset = train)
QdaLdaPredictionError(model = qda.fit, test.data = test.set)

#KNN model ####
knn.pred <-
  knn(mydat[train, ][, 1:10], test.set[, 1:10], mydat$tissue[train], k =
        1)
table(knn.pred, test.set$tissue)
mean(knn.pred == test.set$tissue)

##Adding a non-linear term to lda (best performance) ####
log.term <- log(mydat$ENSG00000131771.9_PPP1R1B)
log.added <- data.frame(mydat, log.term)
set.seed(1)
train <- sample(1:nrow(log.added), round(nrow(log.added) * 0.35))
test.set <- mydat[-train, ]

#LDA model####
lda.fit <- lda(tissue ~ .,
               data = log.added,
               subset = train)
QdaLdaPredictionError(model = lda.fit, test.data = test.set)

#GLM model####
glm.fit <- glm(tissue ~ .,
               data = log.added,
               family = binomial,
               subset = train)
GlmPredictionErrorRate(glm.fit = glm.fit, tissue.1 = tissue1, tissue.2 = tissue2,test.data = test.set)


#QDA model ####
qda.fit <- qda(tissue ~ .,
               data = log.added,
               subset = train)
QdaLdaPredictionError(model = qda.fit, test.data = test.set)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#################
#### Week 2 #####
#################

# 1. Apply at least two tree-based methods, and an SVM with two types of kernel
# (all presented this week) to some of the problems analysed last week. Do so
# for one regression and one classification problem. Compare the predictive
# performance of last weeks models to that of the new models using a
# cross-validation setup. Discuss the observed differences in performance of
# methods for the various classification/regression problems.


#### Tree-based methods  ####
tissue1 <- "brain_cerebellum"
tissue2 <- "brain_amygdala"

#specific data set
set.seed(1)
mydat <- tissue.selection(tissue1, tissue2, expr4T.filtered)


train.mydat <- sample(1:nrow(mydat), round(nrow(mydat) * 0.35))
test.set <- mydat[-train, ]

# Method 1 ####
# Create the tree for the classification problem. ####
tree.tissues <- tree(tissue ~ .,
                     data = mydat,
                     subset = train.mydat)
summary (tree.tissues)
plot(tree.tissues)
text(tree.tissues, pretty = 0)

# Prune the tree for the classification problem.
prune.tree.tissues <- tree.pruner(tree.tissues)
plot(prune.tree.tissues)
text(prune.tree.tissues, pretty = 0)

# Error rate of pruned tree
error.rate.prediction.trees(
  tree.data = prune.tree.tissues,
  dataset = mydat,
  test.set = -train.mydat
)
# Error rate of unpruned tree
error.rate.prediction.trees(
  tree.data = tree.tissues,
  dataset = mydat,
  test.set = -train.mydat
)

# Regression ####
tree.regression <- tree(ENSG00000271043.1_MTRNR2L2 ~ . - tissue,
                        data = mydat,
                        subset = train.mydat)

prune.tree.regression <-
  tree.pruner(tree.regression, regression.tree = TRUE)

error.rate.prediction.trees(
  tree.data = tree.regression,
  dataset = mydat,
  test.set = -train.mydat,
  type.prediction = "vector"
)
error.rate.prediction.trees(
  tree.data = prune.tree.regression,
  dataset = mydat,
  test.set = -train.mydat,
  type.prediction = "vector"
)

# Method 2 ####

p <- ncol(mydat) - 1
p.2 <- p / 2
p.3 <- sqrt(p)


p.error.classification <-
  DetermineEffectPsRandomForest(
    dataset = mydat,
    train.set = train.mydat,
    test.set = -train.mydat,
    column.index = 152
  )
p.error.regression <-
  DetermineEffectPsRandomForest(
    dataset = mydat,
    train.set = train.mydat,
    test.set = -train.mydat,
    column.index = grep("ENSG00000271043.1_MTRNR2L2", colnames(mydat))
  )

###Giving NA values, don't know why. The code should work since is working for you.
mean(p.error.classification$p1.mse)
mean(p.error.classification$p2.mse) # Lowest MSE so it is going to be used for classification.
mean(p.error.classification$p3.mse)

# Perform a randomForest on the tissues to see which predictor is important for the tissues.
rf.tissues <- randomForest(
  tissue ~ .,
  data = mydat,
  subset = train.mydat,
  mtry = p.2,
  ntree = 500,
  importance = TRUE
)
error.rate.prediction.trees(
  tree.data = rf.tissues,
  dataset = mydat,
  test.set = -train.mydat
)

rf.tissues

#### SVM with two types of kernel  ####

# Searching for the best cost ####
set.seed(42)
tune.out.linear <- tune(
  svm,
  tissue ~ .,
  data = mydat ,
  kernel = "linear",
  ranges = list(cost = seq(0.1, 1.0, .1))
)
tune.out.linear

# Perform a SVM using a linear kernel on the classification problem. ####
svm.linear <-
  svm(
    tissue ~ .,
    data = mydat,
    subset = train.mydat,
    kernel = "linear",
    cost = 0.1
  )
summary(svm.linear)

# Performing from c(1, 10, 20, 30, 40, 50, 60), lead to 40.
# c(40, 41, 42, 43, 44, 46, 48) leads to 43
# seq(42.7, 44, 0.1) leads to 43.2
tune.out.poly <- tune(
  svm,
  tissue ~ .,
  data = mydat ,
  kernel = "polynomial",
  ranges = list(cost = seq(1, 50, 1))
)
tune.out.poly

# Perform a SVM using a polynomial kernel on the classification problem.
svm.poly <-   svm(
  tissue ~ .,
  data = mydat,
  subset = train.mydat,
  kernel = "polynomial",
  cost = 23 #In the case of cerebellum and amygdala with 1 is enough
)
summary(svm.poly)

# Error Rates ####
# Error rate SVM linear classification problem.
error.rate.prediction.trees(svm.linear, dataset = mydat, test.set = -train.mydat)

# Error rate SVM polynomial classification problem.
error.rate.prediction.trees(svm.poly, dataset = mydat, test.set = -train.mydat)

# Regression ####
svm.linear.r <-
  svm(
    ENSG00000104888.5_SLC17A7 ~ . - tissue,
    data = mydat,
    subset = train.mydat,
    kernel = "linear",
    cost = 1
  )


# @TODO Error rate of 1, check it again.
error.rate.prediction.trees(svm.linear.r,
                            dataset = mydat,
                            test.set = -train.mydat)

# Perform a SVM using a polynomial kernel on the regression problem.
svm.poly.r <-   svm(
  ENSG00000104888.5_SLC17A7 ~ . - tissue,
  data = mydat,
  subset = train.mydat,
  kernel = "polynomial",
  cost = 1
)

error.rate.prediction.trees(svm.poly.r,
                            dataset = mydat,
                            test.set = -train.mydat)


# RandomForest find informative features  #####

importance(rf.tissues)
varImpPlot(rf.tissues, main = "Cerebellum vs amygdala")

# The results show that for all the trees considered in the rf.tissues, the
# genes ENSG00000221890.2_NPTXR & ENSG00000183379.4_SYNDIG1L are by the two most
# important variables. You base this upon their positions in the top of the
# MeanDecreaseAccuracy and MeanDecreaseGini. So using them to discriminate the
# two different tissues could prove to be usefull.

# Sort on the MeanDecreaseAccuracy to get the best and worst genes.
best.genes <-
  names(sort(importance(rf.tissues)[, 3], decreasing = T))[1:5]
worst.genes <-
  names(sort(importance(rf.tissues)[, 3], decreasing = F))[1:5]

# Calcualte the p, now based on the length.
p <- length(best.genes) - 1
p.2 <- p / 2

# GLM find informative features ###

#The glmnet() function has an alpha argument that determines what type
#of model is fit. If alpha=0 then a ridge regression model is fit, and if alpha=1
#then a lasso model is fit.

x <- model.matrix(tissue~.,mydat)
y <- mydat$tissue
grid =seq (from= 0.1, to =0, by = -0.001)

#Ridge regression ####
#Also, using lambda = grid we have a wide range of lambda to compare with (possibility of a graph!)
ridge.fit <- glmnet(x[train.mydat,],y[train.mydat], alpha=0, lambda = grid, family = "binomial")
plot(ridge.fit, main= "ridge classification", label = T) ##The plot shows that the model is better without the penalty!
ridge.pred <-predict(ridge.fit, newx = x[-train.mydat,])

length(coef(ridge.fit))


#Plot of predictions of two tissues (only if MULTINOMIAL!) ####
x <- model.matrix(tissue~.,expr4T.filtered)
y <- expr4T.filtered$tissue
grid =10^ seq (10,-2, length =100)

ridge.fit <- glmnet(x[train.expr4T.data,],y[train.expr4T.data], alpha=0, lambda = grid, family = "multinomial")
plot(ridge.fit)
ridge.pred <-predict(ridge.fit, type = "coefficients")

plot(ridge.pred$brain_amygdala, ridge.pred$brain_cerebellum)
quantile(ridge.pred)
ridge.pred

#Here we do again a ridge regression but with a specific gene ####
x <- model.matrix(ENSG00000229344.1_RP5.857K21.7~.-tissue, mydat)
y<- mydat$ENSG00000229344.1_RP5.857K21.7
grid =10^ seq (10,-2, length =100)
ridge.fit <- glmnet(x[train.mydat,],y[train.mydat], alpha=0, lambda = grid)
plot(ridge.fit, main = "Regression ridge", label = T)
ridge.pred <- predict(ridge.fit, newx = x[-train.mydat,], s = 0.01) ##Changes s, MSE is different
mean((ridge.pred - y[-train.mydat])^2) ##test MSE


##repeated model only with the three significant genes ####

x <- model.matrix(ENSG00000271043.1_MTRNR2L2~ENSG00000101200.5_AVP+ENSG00000183379.4_SYNDIG1L+ENSG00000101405.3_OXT+ENSG00000209082.1_MT.TL1+ENSG00000130643.4_CALY, mydat) #remove tissue to see differences!
y<- mydat$ENSG00000271043.1_MTRNR2L2
grid =10^ seq (10,-2, length =100)
ridge.fit <- glmnet(x[train.mydat,],y[train.mydat], alpha=1, lambda = grid)
plot(ridge.fit, main = "Regression ridge", label = T)
ridge.pred <- predict(ridge.fit, newx = x[-train.mydat,], s = 0.01) ##Changes s, MSE is different
mean((ridge.pred - y[-train.mydat])^2) ##test MSE



ridge.fit
mydat[,147:152]
table(mydat$tissue)
#The MSE increments with its value, bigger lambda, bigger MSE. If default values is used the MSE is very high (a lot)
##INCLUDE THIS IN THE REPORT!!! WITH A NICE GRAPH!!! WE choose s=0.01 because gives the lowest MSE



cv.out <- cv.glmnet(x[train.mydat,], y[train.mydat], alpha = 0, type.measure = "mse")
plot(cv.out)

ridge.pred = predict(ridge.fit, s = cv.out$lambda.min, newx = x[-train,])
mean((ridge.pred - y[-train])^2)

##We try lasso now ####
##Classification, not useful ####
x <- model.matrix(tissue~., mydat)
y <- mydat$tissue
grid =10^ seq (10,-2, length =100)


ridge.fit <- glmnet(x[train.mydat,],y[train.mydat], alpha=1, lambda = grid, family = "binomial")
plot(ridge.fit, main= "lasso", label = T) 
ridge.pred <-predict(ridge.fit, newx = x[-train.mydat,])

# Lasso regression ####
x <- model.matrix(ENSG00000229344.1_RP5.857K21.7~.-tissue, mydat)
y<- mydat$ENSG00000229344.1_RP5.857K21.7
grid =10^ seq (10,-2, length =100)
lasso.fit <- glmnet(x[train.mydat,],y[train.mydat], alpha=1, lambda = grid)
plot(lasso.fit, main = "Regression lasso", label = T)
lasso.pred <- predict(lasso.fit, newx = x[-train.mydat,], s = 0.01) ##Changes s, MSE is different. 0.01 seems to give the lowest MSE.
mean((lasso.pred - y[-train.mydat])^2) 

cv.out <- cv.glmnet(x[train.mydat,], y[train.mydat], alpha = 1, type.measure = "mse")
plot(cv.out, main = "Lasso regression")

lasso.pred = predict(lasso.fit, s = cv.out$lambda.min, newx = x[-train,])
mean((lasso.pred - y[-train])^2)


# Generate a random forest of the best genes for the classification problem. ####
rf.best.genes <- randomForest(
  tissue ~ ENSG00000258283.1_RP11.386G11.3 + ENSG00000100362.8_PVALB + ENSG00000198121.9_LPAR1 +ENSG00000139899.6_CBLN3 +ENSG00000165802.15_NSMF,
  data = mydat,
  subset = train.mydat,
  mtry = p.2,
  ntree = 500,
  importance = TRUE
)
error.rate.prediction.trees(
  tree.data = rf.best.genes,
  dataset = mydat,
  test.set = -train.mydat
)

# Generate a random forest of the worst genes for the classification problem. ####
rf.worst.genes <- randomForest(
  tissue ~ ENSG00000225972.1_MTND1P23 + ENSG00000225630.1_MTND2P28 +ENSG00000237973.1_hsa.mir.6723 + ENSG00000229344.1_RP5.857K21.7 + ENSG00000126709.10_IFI6 ,
  data = mydat,
  subset = train.mydat,
  mtry = p.2,
  ntree = 500,
  importance = TRUE
)
error.rate.prediction.trees(
  tree.data = rf.worst.genes,
  dataset = mydat,
  test.set = -train.mydat
)


###############
# Question 3  #

set.seed(1)

# Create a tree for the set of best genes.####
tree.best.genes <-
  tree(
    tissue ~ ENSG00000258283.1_RP11.386G11.3 + ENSG00000100362.8_PVALB + ENSG00000198121.9_LPAR1 +ENSG00000139899.6_CBLN3 +ENSG00000165802.15_NSMF,
    data = mydat,
    subset = train.mydat
  )
plot(tree.best.genes)
text(tree.best.genes, pretty = 0)

# Prune the best genes tree.
prune.tree.best.genes <- tree.pruner(tree.best.genes)
plot(prune.tree.best.genes)
text(prune.tree.best.genes, pretty = 0)

# Error rate of pruned tree.
error.rate.prediction.trees(
  tree.data = prune.tree.best.genes,
  dataset = mydat,
  test.set = -train.mydat
)
# Error rate of unpruned tree
error.rate.prediction.trees(
  tree.data = tree.best.genes,
  dataset = mydat,
  test.set = -train.mydat
)

# Error rate tree.tissue: 0.03649635

#The same error rate.

# SVM on the best features.
svm.linear.best <-
  svm(
    tissue ~ ENSG00000258283.1_RP11.386G11.3+ ENSG00000100362.8_PVALB + ENSG00000198121.9_LPAR1 + ENSG00000139899.6_CBLN3 + ENSG00000165802.15_NSMF,
    data = mydat,
    subset = train.mydat,
    kernel = "linear",
    cost = 0.1
  )

svm.poly.best <-
  svm(
    tissue ~ ENSG00000258283.1_RP11.386G11.3+ ENSG00000100362.8_PVALB + ENSG00000198121.9_LPAR1 + ENSG00000139899.6_CBLN3 + ENSG00000165802.15_NSMF,
    data = mydat,
    subset = train.mydat,
    kernel = "polynomial",
    cost = 1
  )

error.rate.prediction.trees(
  tree.data = svm.linear.best,
  dataset = mydat,
  test.set = -train.mydat
)

error.rate.prediction.trees(
  tree.data = svm.poly.best,
  dataset = mydat,
  test.set = -train.mydat
)


rf.best.genes <- randomForest(
  tissue ~ ENSG00000258283.1_RP11.386G11.3 + ENSG00000100362.8_PVALB + ENSG00000198121.9_LPAR1 +ENSG00000139899.6_CBLN3 +ENSG00000165802.15_NSMF,
  data = mydat,
  subset = train.mydat,
  mtry = (5 - 1) / 2,
  ntree = 500,
  importance = TRUE
)
rf.best.genes  # Higher OOB estimate of  error rate, also in the case of cerebellum and hemsphere example
rf.tissues


lm.fit.two.tissues.best.genes <-
  lm(
    ENSG00000271043.1_MTRNR2L2 ~ ENSG00000258283.1_RP11.386G11.3 + ENSG00000100362.8_PVALB + ENSG00000198121.9_LPAR1 +ENSG00000139899.6_CBLN3 +ENSG00000165802.15_NSMF,
    data = mydat,
    subset = train.mydat
  )

plot(lm.fit.two.tissues.best.genes)
preds <- predict(lm.fit.two.tissues.best.genes, newdata = mydat[- train.mydat,])

plot(mydat$ENSG00000271043.1_MTRNR2L2[- train.mydat], preds,
     xlab = "Predicted values for ENSG00000271043.1_MTRNR2L2", 
     ylab = "Observed values for ENSG00000271043.1_MTRNR2L2")
abline(0,1, col = "red")

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = mydat,
  training.data = train.mydat,
  lm.training = lm.fit.two.tissues.best.genes
)
# what is we transform the data with log again as last week?????? ####

log.mydat <- log(mydat[-ncol(mydat)])

is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) all(is.finite(x)))
}

# Check if there are any infs in the log-transformed datasets.
is.finite.data.frame(log.mydat)


# Works only on matrixes
log.mydat <- as.matrix(log.mydat)

# Replace any non finites with 0, this is needed to get the logs to work.
log.mydat[!is.finite(log.mydat)] <- 0


# Converted back to a data.frame for the lm's.
log.mydat <- as.data.frame(log.mydat)

lm.fit.two.tissues.best.genes.log <-
  lm(
    ENSG00000271043.1_MTRNR2L2  ~ ENSG00000258283.1_RP11.386G11.3 + ENSG00000100362.8_PVALB + ENSG00000198121.9_LPAR1 +ENSG00000139899.6_CBLN3 +ENSG00000165802.15_NSMF,
    data = log.mydat,
    subset = train.mydat
  )

plot(lm.fit.two.tissues.best.genes.log)
preds <- predict(lm.fit.two.tissues.best.genes.log, newdata = mydat[- train.mydat,])

plot(log.mydat$ENSG00000271043.1_MTRNR2L2 [- train.mydat], preds,
     xlab = "Predicted values for ENSG00000271043.1_MTRNR2L2", 
     ylab = "Observed values for ENSG00000271043.1_MTRNR2L2",
     main = "Log transformation best genes")
abline(0,1, col = "red")

PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = log.mydat,
  training.data = train.mydat,
  lm.training = lm.fit.two.tissues.best.genes.log
)


# R-square =  -0.03732637
# Fraction of variability explained by the model =  0.0712916


###############
# Question 4  #

set.seed(1)
# Create a tree for everything besides the best features
tree.all.besides.best <-
  tree(
    tissue ~ . - ENSG00000258283.1_RP11.386G11.3+-ENSG00000100362.8_PVALB+-ENSG00000198121.9_LPAR1+-ENSG00000139899.6_CBLN3+-ENSG00000165802.15_NSMF,
    data = mydat,
    subset = train.mydat
  )

# Prune the tree.
prune.all.besides.best <- tree.pruner(tree.all.besides.best)

# Error rate of pruned tree.
error.rate.prediction.trees(
  tree.data = prune.all.besides.best,
  dataset = mydat,
  test.set = -train.mydat
)
# Error rate of unpruned tree
error.rate.prediction.trees(
  tree.data = tree.all.besides.best,
  dataset = mydat,
  test.set = -train.mydat
)


svm.linear.all.besides.best <-
  svm(
    tissue ~ . - ENSG00000258283.1_RP11.386G11.3+-ENSG00000100362.8_PVALB+-ENSG00000198121.9_LPAR1+-ENSG00000139899.6_CBLN3+-ENSG00000165802.15_NSMF,
    data = mydat,
    subset = train.mydat,
    kernel = "linear",
    cost = 0.1
  )

svm.poly.all.besides.best <-
  svm(
    tissue ~ . - ENSG00000258283.1_RP11.386G11.3+-ENSG00000100362.8_PVALB+-ENSG00000198121.9_LPAR1+-ENSG00000139899.6_CBLN3+-ENSG00000165802.15_NSMF,
    data = mydat,
    subset = train.mydat,
    kernel = "polynomial",
    cost = 1
  )

# Error rate SVM linear classification problem for worst genes.
error.rate.prediction.trees(svm.linear.all.besides.best,
                            dataset = mydat,
                            test.set = -train.mydat)

# Error rate SVM polynomial classification problem for worst genes.
error.rate.prediction.trees(svm.poly.all.besides.best,
                            dataset = mydat,
                            test.set = -train.mydat)

rf.all.besides.best <- randomForest(
  tissue ~ . - ENSG00000258283.1_RP11.386G11.3+-ENSG00000100362.8_PVALB+-ENSG00000198121.9_LPAR1+-ENSG00000139899.6_CBLN3+-ENSG00000165802.15_NSMF,
  data = mydat,
  subset = train.mydat,
  mtry = (5 - 1) / 2,
  ntree = 500,
  importance = TRUE
)
error.rate.prediction.trees(rf.all.besides.best,
                            dataset = mydat,
                            test.set = -train.mydat)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


#################
#### Week 3 #####
#################

# 1. Apply PCA to the expression data from the brain subregions. Perform PCA
# both with and without scaling to unit standard deviation. Address the
# following questions:

# Perform a PCA on the data, x must be numeric so remove the tissues from the PCA.
pr.out.scale <- prcomp(mydat[-152], scale = TRUE)
pr.out <- prcomp(mydat[-152], scale = FALSE)

# a. Analyze how much of the variation is explained by the first few PCs.

# Calculate the variance explained by each principal component, do this for the scale and non-scale.
pr.var.scale <- pr.out.scale$sdev ^ 2
pr.var <- pr.out$sdev ^ 2

# Calculate the proportion of the variance that is being explained by each principal component.
pve.scale <- pr.var.scale / sum(pr.var.scale)
pve <- pr.var / sum(pr.var)

# Plot the PVE for each component and the cumulative PVE as well.
par(mfrow = c(2, 2))

# Plot the scaled data.
plot(
  pve.scale,
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Proportion of Variance Explained Scaled PCA",
  ylim = c(0, 1),
  type = 'b'
)
plot(
  cumsum(pve.scale),
  xlab = "Principal Component",
  ylab = "Cumulative Proportion of Variance Explained",
  main = "Cumulative Proportion of Variance Explained Scaled PCA",
  ylim = c(0, 1) ,
  type = 'b'
)
# Plot the non-scaled data.
plot(
  pve,
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Proportion of Variance Explained Unscaled PCA",
  ylim = c(0, 1),
  type = 'b'
)
plot(
  cumsum(pve),
  xlab = "Principal Component",
  ylab = "Cumulative Proportion of Variance Explained",
  main = "Cumulative Proportion of Variance Explained Unscaled PCA",
  ylim = c(0, 1) ,
  type = 'b'
)

# There is quite a difference between the scale and the non-scale, the first PC
# in the non-scale explains 81.18%. While the scale only explains 27.68%.

# b. How are the different tissues placed in the space defined by
# these PCs? Comment on how this relates to observations when you classified
# these tissues compared to each other (week 1,2).

# See the number of unique tissues
unique(expr4T.filtered$tissue)

# Each entry corrosponds to the same index as the tissue, so first tissue has the first color.
sapply(Cols(unique(expr4T.filtered$tissue)), color.id)

#plotting samples after transformation:

par(mfrow = c(1, 2))

plot(pr.out$x[, 1:2],
     main = "PCA Unscaled Tissue Placement",
     col = Cols(expr4T.filtered$tissue),
     pch = 19)

# pr.out$rotation



plot(pr.out.scale$x[, 1:2],
     main = "PCA Scaled Tissue Placement",
     col = Cols(expr4T.filtered$tissue),
     pch = 19)

# pr.out.scale$rotation

# cbind(pr.out.scale$x[, 1:2], expr4T.filtered$tissue)
# cbind(pr.out$x[, 1:2], expr4T.filtered$tissue)



# Which genes contribute substantially to first or second PC:

pr.out.scale$rotation[which(abs(pr.out.scale$rotation[, 1]) > 0.2 |
                              abs(pr.out.scale$rotation[, 2]) > 0.2), 1:2]

# Which genes contribute substantially to first or second PC:

pr.out$rotation[which(abs(pr.out$rotation[, 1]) > 0.2 |
                        abs(pr.out$rotation[, 2]) > 0.2), 1:2]

# c. Use principal component loading vectors to identify genes which might be relevant for separation of
# different tissues. Compare these with your observation on informative features in week 2.

# @ TODO: Compare with week 2

# Look at which genes contribute substantially to each principal component loading vectors.
pr.out$rotation[which(abs(pr.out$rotation[, 1]) > 0.2), ]
pr.out.scale$rotation[which(abs(pr.out.scale$rotation[, 1]) > 0.2), ]



# 2. Apply hierarchical clustering, testing two different distance methods and
# two different linkages. Also apply K-means clustering. Finally, use the PCA
# results from step 1 (above), pick the first few principal components, and
# apply clustering to the principal component score vectors. Address the
# following questions:

# a. Compare the different clustering results. Do different methods give very
# different or very similar results?


# Use all the data for clustering, besides the tissue data.
clustering.data <- expr4T.filtered[-152]

# Set the number of centers equal to the number of unique tissues.
ncenters <- length(unique(expr4T.filtered$tissue))

# Perform a kmeans clustering on the dataset, use different settings to see the effect of them.
km.out.scale.13 <-
  kmeans(scale(clustering.data), centers = ncenters)
km.out.13 <- kmeans(clustering.data, centers = ncenters)
km.out.nstart.50 <-
  kmeans(clustering.data, centers = ncenters, nstart = 50)
km.out.nstart.50.scale <-
  kmeans(scale(clustering.data),
         centers = ncenters,
         nstart = 50)
km.out.13.pca <- kmeans(pr.out$x[, 1:2], centers = ncenters)
km.out.nstart.50.pca <-
  kmeans(pr.out$x[, 1:2], centers = ncenters, nstart = 20)

table(km.out.scale.13$cluster)
table(km.out.13$cluster)
# They give very different results for each cluster.

table(km.out.nstart.50$cluster)
table(km.out.nstart.50.scale$cluster)
# They also give very different results for each cluster.

table(km.out.13.pca$cluster)
table(km.out.nstart.50.pca$cluster)
# They also give very different results for each cluster.

# Perform hierarchical clustering, using three methods on the pca data.
hc.complete.pca <- hclust(dist(pr.out$x[, 1:2]), method = "complete")
hc.average.pca <- hclust(dist(pr.out$x[, 1:2]), method = "average")
hc.single.pca <- hclust(dist(pr.out$x[, 1:2]), method = "single")

# Perform hierarchical clustering, using three methods on the clustering data.
hc.complete <- hclust(dist(clustering.data), method = "complete")
hc.average <- hclust(dist(clustering.data), method = "average")
hc.single <- hclust(dist(clustering.data), method = "single")

# See the groups present within the dendrogram.
cutree(hc.complete , ncenters)
cutree(hc.average , ncenters)
cutree(hc.single , ncenters)

# b. Compare the clustering results with the results obtained when classifying
# tissues in weeks 1 and 2.

# K-means generate class labels, each of these class labels corrosponds to a entry in the tissues.
# So unique(expr4T.filtered$tissue)[1] corrosponds to k-means class label 1.
tissues.clustering.assigned <-
  unique(expr4T.filtered$tissue)[as.vector(km.out.13$cluster)]

# Put them into a data.frame, the tissues.clustering.assigned as a new column.
dat <- data.frame(expr4T.filtered, K = tissues.clustering.assigned)

# See the results of the K-means compared to the real tissues.
table(dat$K)
table(dat$tissue)

# Plot to see how much difference there is between the real tissues and the assigned k-means class labels.
plot(
  as.vector(table(dat$tissue)) - as.vector(table(dat$K)),
  main = "Differences between the class labels",
  ylab = "Difference true class labels and clustering labels",
  xlab = "Class label",
  type = "b"
)

# @ TODO: compare the results with the week 1 and 2.

# c. Can you say something on how many clusters there are in this dataset?

#http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters

# Calculate the within groups sum of squares
wss <-
  (nrow(clustering.data) - 1) * sum(apply(clustering.data, 2, var))

# Loop trough all centers that have been used to create the clustering.
for (i in 1:ncenters) {
  wss[i] <- sum(kmeans(clustering.data,
                       centers = i)$withinss)
}

# Plot the SSE plot to see if clusters are present within the data.
plot(
  1:ncenters,
  wss,
  type = "b",
  main = "Compare SSE for a number of clusters",
  xlab = "Number of Clusters",
  ylab = "Within groups sum of squares"
)

# There seems to be clusters in the data, this can be seen in the decrease of the within groups sum of squares.
# The number of clusters seem to be 6, since there is no longer a sharp decrease of SSE. AFter ti the reduction slows down quite a bit.

# 3. Use one of the approaches from step 1 or 2 above to pre-process the gene
# expression data. After this preprocessing, re-train at least two of the
# classification approaches applied in week 1 or 2. Compare the prediction
# performance with what was obtained previously.

# Use clustering to find a intreseting group of tissues, you can train a model of that and use test data to validate it.

interesting.genes <- c()
# Loop trough all the centers in the tree.
for (i in 1:ncenters) {
  # If there are just 2 group members within a single cut, than save it.
  if (length(which(cutree(hc.average , ncenters) == i)) <= 2) {
    # Take the name of the  interesting gene.
    interesting.genes <-
      c(interesting.genes, names(which(cutree(
        hc.average , ncenters
      ) == i)))
  }
}

# Perform a randomforst on the pre-processed data of the clustering, predictors were taken from the interesting.genes vector.
rf.tissues.clustering.preprocess <- randomForest(
  tissue ~ ENSG00000225972.1_MTND1P23 + ENSG00000225630.1_MTND2P28 + ENSG00000237973.1_hsa.mir.6723 + ENSG00000229344.1_RP5.857K21.7 + ENSG00000154146.8_NRGN + ENSG00000131095.7_GFAP + ENSG00000266844.1_RP11.862L9.3 + ENSG00000197971.10_MBP + ENSG00000123560.9_PLP1 + ENSG00000226958.1_CTD.2328D6.1 + ENSG00000198695.2_MT.ND6 + ENSG00000101405.3_OXT + ENSG00000101200.5_AVP,
  data = expr4T.filtered,
  subset = train.expr4T.data,
  mtry = (length(interesting.genes) - 1) / 2,
  # p.2
  ntree = 500,
  importance = TRUE
)
# OOB estimate of  error rate: 43.17%

error.rate.prediction.trees(
  rf.tissues.clustering.preprocess,
  dataset = expr4T.filtered,
  test.set = -train.expr4T.data
)


# rf.clus.pred <-
#   predict (rf.tissues.clustering.preprocess, newdata = expr4T.filtered[-train.expr4T.data,])
# 1 - mean(rf.clus.pred == expr4T.filtered[-train.expr4T.data,]$tissue)

# Perform a randomforst on the 'normal'data, predictors were taken from the all the genes.
rf.tissues.all.genes <- randomForest(
  tissue ~ .,
  data = expr4T.filtered,
  subset = train.expr4T.data,
  mtry = (ncol(expr4T.filtered) - 1) / 2,
  # p.2
  ntree = 500,
  importance = TRUE
)
# OOB estimate of  error rate: 16.67%

error.rate.prediction.trees(rf.tissues.all.genes,
                            dataset = expr4T.filtered,
                            test.set = -train.expr4T.data)


# rf.pred <-
#   predict (rf.tissues.all.genes, newdata = expr4T.filtered[-train.expr4T.data,])
# 1 - mean(rf.pred == expr4T.filtered[-train.expr4T.data,]$tissue)

# Performance of the randomforest

# Create a tree that is based on the predictors of the pre-processed data.
tree.expr.clus <-
  tree(
    tissue ~ ENSG00000225972.1_MTND1P23 + ENSG00000225630.1_MTND2P28 + ENSG00000237973.1_hsa.mir.6723 + ENSG00000229344.1_RP5.857K21.7 + ENSG00000154146.8_NRGN + ENSG00000131095.7_GFAP + ENSG00000266844.1_RP11.862L9.3 + ENSG00000197971.10_MBP + ENSG00000123560.9_PLP1 + ENSG00000226958.1_CTD.2328D6.1 + ENSG00000198695.2_MT.ND6 + ENSG00000101405.3_OXT + ENSG00000101200.5_AVP,
    data = expr4T.filtered,
    subset = train.expr4T.data
  )

prune.expr.clus <- tree.pruner(tree.expr.clus)


tree.all.genes <- tree(tissue ~ .,
                       data = expr4T.filtered,
                       subset = train.expr4T.data)


prune.expr.all.genes <- tree.pruner(tree.all.genes)

# Error rate unpruned tree of the clustering data.
error.rate.prediction.trees(
  tree.data = tree.expr.clus,
  dataset = expr4T.filtered,
  test.set = -train.expr4T.data
)

# Error rate pruned tree of the clustering data.
error.rate.prediction.trees(
  tree.data = prune.expr.clus,
  dataset = expr4T.filtered,
  test.set = -train.expr4T.data
)

# Error rate of unpruned all gene data.
error.rate.prediction.trees(
  tree.data = tree.all.genes,
  dataset = expr4T.filtered,
  test.set = -train.expr4T.data
)

# Error rate of pruned all gene data.
error.rate.prediction.trees(
  tree.data = prune.expr.all.genes,
  dataset = expr4T.filtered,
  test.set = -train.expr4T.data
)
