rm(list = ls())
library(tidyverse)
library(caret)
library(ROCR)
library(InformationValue)
library(pscl)
library(e1071)

## yeast dataset
# A imbalanced version of the Yeast data set, 
# where the positive examples belong to class ME3 and the negative examples belong to the rest.
yeast_dat = read.table("yeast1.dat",sep = ",",skip = 15) %>% as_tibble()

# convert target variable to an integer
yeast_dat = yeast_dat %>%
            mutate(target = ifelse(V9 ==" positive",1,0)) %>%
            mutate(V9 = NULL)

# check there are no strings or characters in any column
str(yeast_dat)

# create test and train data
index = unlist(createDataPartition(yeast_dat$target,p =  0.7))
train = yeast_dat[index,]
test = yeast_dat[-index,]


# run logisitc regression for classification
reg1 = glm(target ~ .,data = train,family='binomial')

# summary of the model
summary(reg1)

## checking accuracy of glm
##########################################################
# get probabilty scores for validation set
preds = predict(reg1, test %>% select(-target), steps = 1)

ROCRpred = prediction(preds, test$target)
ROCRperf = performance(ROCRpred, 'tpr','fpr')
plot(ROCRperf)

# get AUROC
auc = performance(ROCRpred, measure = "auc")
auc = auc@y.values[[1]]
auc

# probablity density plot
dat = data.frame(dens = preds, lines = test$target)
densityplot(~dens,data=dat,groups = lines,
            plot.points = FALSE, ref = TRUE, 
            auto.key = list(space = "right"))

# get confusion matrix
pred = ifelse (preds > 0, 1, 0)
caret::confusionMatrix(as.factor(pred), as.factor(test$target), positive = "1")
############################################################################################





# cascaded classifier iterations (used glm) - 2 level cascade
##########################################################
###### implementing cascade (level 1)
## fit training set
train$predictions = predict(reg1, train %>% select(-target))

#################################################
## get density plots to determine threshold
dat = data.frame(dens = train$predictions, lines = train$target)
densityplot(~dens,data=dat,groups = lines,
            plot.points = FALSE, ref = TRUE,
            auto.key = list(space = "right"))

# opt =  optimalCutoff(train$target, train$predictions, optimiseFor = "missclasserror", returnDiagnostics = FALSE)


## classify based on threshold
train$pred_yeast = ifelse(train$predictions > 0, 1, 0)
#################################################

cfm_level_1 = caret::confusionMatrix(as.factor(train$pred_yeast),
                                     as.factor(train$target), positive = "1")
cfm_level_1

#################################################
## removing correctly predicted 0's 2nd level cascade
## for train
train_lvl2 = train %>%
  filter((pred_yeast == 1) | (target == 1))

train_lvl2 = train_lvl2 %>%
  select(-c(pred_yeast, predictions))



################################################
## for test
test$predictions = predict(reg1, test %>% select(-target))


dat = data.frame(dens = test$predictions
                 , lines = test$target)
densityplot(~dens,data=dat,groups = lines,
            plot.points = FALSE, ref = TRUE,
            auto.key = list(space = "right"))

opt = optimalCutoff(test[,c("target")], test$predictions, optimiseFor = "misclasserror", returnDiagnostics = FALSE)


## classify based on threshold (considered less false negatives - yeasts classified as non yeasts)
test$pred_yeast = ifelse(test$predictions > opt, 1, 0)
#################################################

cfm_level_1_test = caret::confusionMatrix(as.factor(test$pred_yeast),
                                          as.factor(test$target), positive = "1")
cfm_level_1_test

#################################################
test_lvl2 = test %>%
  filter((pred_yeast == 1))

test_lvl2 = test_lvl2 %>%
  select(-c(pred_yeast, predictions))

##############################################

test_check = test

assign_labels = function(a, p){
  # print(a)
  if((a == 1) & (p == 1)){
    return("TP")}
  else if((a == 1) & (p == 0)){
    return("FN")}
  else if((a == 0) & (p == 1)){
    return("FP")}
  else if((a == 0) & (p == 0))
    return("TN")}

test_check = test_check %>%
  rowwise()%>%
  mutate(labels = assign_labels(target, pred_yeast)) %>%
  filter((labels == "TN") | (labels == "FN"))


###############################################

dim(train)
dim(train_lvl2)




# level 2 - cascade
##########################################################

## fitting
##########################################################

reg2 = glm(target ~ .,data = train_lvl2,family='binomial')


## high difference in validation and training accuracy hints
## high variance


train_lvl2$predictions = predict(reg2, train_lvl2 %>% select(-c(target)))

#################################################
## get density plots to determine threshold
dat = data.frame(dens = train_lvl2$predictions
                 , lines = train_lvl2$target)
densityplot(~dens,data=dat,groups = lines,
            plot.points = FALSE, ref = TRUE,
            auto.key = list(space = "right"))

opt = optimalCutoff(train_lvl2[,c("target")], train_lvl2$predictions, optimiseFor = "misclasserror", returnDiagnostics = FALSE)


## classify based on threshold
train_lvl2$pred_yeast = ifelse(train_lvl2$predictions > opt, 1, 0)
#################################################

cfm_level_2 = caret::confusionMatrix(as.factor(train_lvl2$pred_yeast),
                                     as.factor(train_lvl2$target), positive = "1")
cfm_level_2

#################################################
## subsetting predicted "1s" for 2nd level cascade
## for train
train_lvl3 = train_lvl2 %>%
  filter((pred_yeast == 1) | (target == 1))

train_lvl3 = train_lvl3 %>%
  select(-c(pred_yeast, predictions))



################################################
## for test
test_lvl2$predictions = predict(reg2, test_lvl2 %>% select(-target))


dat = data.frame(dens = test_lvl2$predictions
                 , lines = test_lvl2$target)
densityplot(~dens,data=dat,groups = lines,
            plot.points = FALSE, ref = TRUE,
            auto.key = list(space = "right"))
## no visible separation of probability densities
## this should be the last level of cascade

opt = optimalCutoff(test_lvl2[,c("target")], test_lvl2$predictions, optimiseFor = "misclasserror", returnDiagnostics = FALSE)

## classify based on threshold
test_lvl2$pred_yeast = ifelse(test_lvl2$predictions > opt, 1, 0)
#################################################

cfm_level_2_test = caret::confusionMatrix(as.factor(test_lvl2$pred_yeast),
                                          as.factor(test_lvl2$target), positive = "1")
cfm_level_2_test

#################################################
test_lvl3 = test_lvl2 %>%
  filter((pred_yeast == 1))


test_lvl3 = test_lvl3 %>%
  select(-c(pred_yeast, predictions))


#################################################
test_check_2 = test_lvl2


test_check_2 = test_check_2 %>%
  rowwise()%>%
  mutate(labels = assign_labels(target, pred_yeast))



test_check_2 = test_check %>%
  rbind(test_check_2)

### comparing cascaded confusion matrix with first iteration confusion matrix
#################################################
cfm_level_1_test


caret::confusionMatrix(as.factor(test_check_2$pred_yeast),
                       as.factor(test_check_2$target), positive = "1")



