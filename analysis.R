# Install versioned packages used for analysis (if needed): ----
# require(remotes)
# install_version('caret', version = '6.0-85', repos = 'http://cran.us.r-project.org')
# install_version('quantregForest', version = '1.3-7', repos = 'http://cran.us.r-project.org')
# install_version('tidyr', version = '1.0.2', repos = 'http://cran.us.r-project.org')
# install_version('ggthemes', version = '4.2.0', repos = 'http://cran.us.r-project.org')
# install_version('ggridges', version = '0.5.2', repos = 'https://cran.us.r-project.org')

# Load Libraries ----
# needed for analysis
library(caret)
library(quantregForest)
# needed for plotting
library(tidyr)
library(ggthemes)
library(ggridges)

# Load Data
data <- read.csv('GRF_IMU_data.csv', header = T)
# Factor Sub ID and Sex 
data$Sub <- as.factor(data$Sub)
data$IsFemale <- as.factor(data$IsFemale)

# Split Data ----
set.seed(541) # set seed for reproducibility
#determine which subjects will be in train/test. 
#Want similar distribution of M/F in test/train AND no subjects in both test/train.
#split by sex first
male.data <- data[data$IsFemale == 0,]
female.data <- data[data$IsFemale == 1,]
folds = 5
#make 80/20 split of each sex according to subject
male.split <- groupKFold(male.data$Sub, k = folds)
female.split <- groupKFold(female.data$Sub, k = folds)
#rejoin male and female data for train/test
train.data <- rbind(male.data[male.split$Fold1,], female.data[female.split$Fold1,]) #will use to do kfold cv with each model
test.data <- rbind(male.data[-male.split$Fold1,], female.data[-female.split$Fold1,]) #test different models on totally new data

#define folds for kfold cross validation (cv) ----
#determine randomly which subjects will be in train/test.
Male <- unique(train.data$Sub[train.data$IsFemale == 0])
Female <- unique(train.data$Sub[train.data$IsFemale == 1])
Male.Folds <- createFolds(1:length(Male), k =folds, returnTrain = TRUE)
Female.Folds <- createFolds(1:length(Female), k =folds, returnTrain = TRUE)

Male.Index <- vector("list", folds)
for(i in 1:length(Male.Folds)){
  x <- which(train.data$Sub %in% Male[Male.Folds[[i]]])
  Male.Index[[i]] <- x
}

Female.Index <- vector("list", folds)
for(i in 1:length(Female.Folds)){
  x <- which(train.data$Sub %in% Female[Female.Folds[[i]]])
  Female.Index[[i]] <- x
}
#join indexes for training data
train.data.Index <- mapply(c,Male.Index, Female.Index, SIMPLIFY = FALSE)
#name folds for caret
names(train.data.Index) <- sapply(1:folds, function(x) paste(c('fold',x), collapse = ''))

#define caret training control parameter
train.control <- trainControl(savePredictions = 'all', index = train.data.Index, method = 'cv', number = folds)

# BUILD MODELS ----
#where to store results
peak.accuracy <- data.frame(matrix(nrow = 1, ncol = 4))
colnames(peak.accuracy) <- c('model', 'RMSE','Rsquared','MAE')
impulse.accuracy <- peak.accuracy
tc.accuracy <- peak.accuracy
# Equations
formula.peak <- GRFPeak ~ Speed+IMUPeak+Mass+IMUsf
formula.impulse <- GRFImpulse ~ Speed+IMUImpulse+Mass+IMUsf
formula.tc <- GRFtc ~ Speed+IMUtc+Mass+IMUsf

# Linear Regression ----
#train model
set.seed(541) 
lm.peak<- train(formula.peak, data = train.data, trControl = train.control, method = 'lm')
lm.impulse<- train(formula.impulse, data = train.data, trControl = train.control, method = 'lm')
lm.tc<- train(formula.tc, data = train.data, trControl = train.control, method = 'lm')

#test model & store results
test.data$lm.peak <- predict(lm.peak, newdata = test.data)
peak.accuracy <- rbind(peak.accuracy,c('lm', postResample(pred = test.data$GRFPeak, obs = test.data$lm.peak)))

test.data$lm.impulse <- predict(lm.impulse, newdata = test.data)
impulse.accuracy <- rbind(impulse.accuracy, c('lm',postResample(pred = test.data$GRFImpulse, obs = test.data$lm.impulse)))

test.data$lm.tc <- predict(lm.tc, newdata = test.data)
tc.accuracy <- rbind(tc.accuracy, c('lm', postResample(pred = test.data$GRFtc, obs = test.data$lm.tc)))

# Quantile Regression/Random Forest ----
set.seed(541)
qrf.peak<- train(formula.peak, data = train.data, trControl = train.control, method = 'qrf')
qrf.impulse<- train(formula.impulse, data = train.data, trControl = train.control, method = 'qrf')
qrf.tc<- train(formula.tc, data = train.data, trControl = train.control, method = 'qrf')

#test model & store results
test.data$qrf.peak <- predict(qrf.peak, newdata = test.data)
peak.accuracy <- rbind(peak.accuracy,c('qrf', postResample(pred = test.data$GRFPeak, obs = test.data$qrf.peak)))

test.data$qrf.impulse <- predict(qrf.impulse, newdata = test.data)
impulse.accuracy <- rbind(impulse.accuracy, c('qrf',postResample(pred = test.data$GRFImpulse, obs = test.data$qrf.impulse)))

test.data$qrf.tc <- predict(qrf.tc, newdata = test.data)
tc.accuracy <- rbind(tc.accuracy, c('qrf', postResample(pred = test.data$GRFtc, obs = test.data$qrf.tc)))

# Store Model Accuracy for Metrics ----
peak.accuracy <- peak.accuracy[-1,]
impulse.accuracy <- impulse.accuracy[-1,]
tc.accuracy <- tc.accuracy[-1,]

peak.accuracy$metric <- 'Peak vGRF'
impulse.accuracy$metric <- 'Vertical Impulse'
tc.accuracy$metric <- 'Contact Time'
 
peak.accuracy$formula <- Reduce(paste, deparse(formula.peak))
impulse.accuracy$formula <- Reduce(paste, deparse(formula.impulse))
tc.accuracy$formula <- Reduce(paste, deparse(formula.tc))
 
model_accuracy <- rbind(peak.accuracy,impulse.accuracy, tc.accuracy)
View(model_accuracy)

# Calculate Average (+/- SD) MAPE and RMSE across test subjects and compare models ----
mape <- function (act, pred)
{
  return(abs(act - pred)/abs(act)*100)
}
# Quantile random forest (qrf)
#Peak vGRF
mean(mape(test.data$GRFPeak, test.data$qrf.peak))
sd(mape(test.data$GRFPeak, test.data$qrf.peak))

#Vertical Impulse
mean(mape(test.data$GRFImpulse, test.data$qrf.impulse))
sd(mape(test.data$GRFImpulse, test.data$qrf.impulse))

#Contact Time
mean(mape(test.data$GRFtc, test.data$qrf.tc))
sd(mape(test.data$GRFtc, test.data$qrf.tc))

# Linear model (lm)
#Peak vGRF
mean(mape(test.data$GRFPeak, test.data$lm.peak))
sd(mape(test.data$GRFPeak, test.data$lm.peak))

#Vertical Impulse
mean(mape(test.data$GRFImpulse, test.data$lm.impulse))
sd(mape(test.data$GRFImpulse, test.data$lm.impulse))

#Contact Time
mean(mape(test.data$GRFtc, test.data$lm.tc))
sd(mape(test.data$GRFtc, test.data$lm.tc))

# Calculate Paired T-Tests Bewtween QRF & LR ----
#QRF
mape.qrf.peak <- (abs(test.data$GRFPeak - test.data$qrf.peak)/abs(test.data$GRFPeak))*100
mape.qrf.impulse <- (abs(test.data$GRFImpulse - test.data$qrf.impulse)/abs(test.data$GRFImpulse))*100
mape.qrf.tc <- (abs(test.data$GRFtc - test.data$qrf.tc)/abs(test.data$GRFtc))*100

#LR
mape.lm.peak <- (abs(test.data$GRFPeak - test.data$lm.peak)/abs(test.data$GRFPeak))*100
mape.lm.impulse <- (abs(test.data$GRFImpulse - test.data$lm.impulse)/abs(test.data$GRFImpulse))*100
mape.lm.tc <- (abs(test.data$GRFtc - test.data$lm.tc)/abs(test.data$GRFtc))*100

# T-Tests
t.test(mape.lm.peak, mape.qrf.peak, paired = TRUE, alternative = 'two.sided')
t.test(mape.lm.impulse, mape.qrf.impulse, paired = TRUE, alternative = 'two.sided')
t.test(mape.lm.tc, mape.qrf.tc, paired = TRUE, alternative = 'two.sided')

# GENERATE FIGURE 1 ----
# QRF Plots ----
test <- test.data
test$Sex1 <- as.numeric(levels(test.data$IsFemale))[test.data$IsFemale] #remove factored sex in test dataset bc caret used dummy variables?

# Reformat models to work with predict.randomForest(predict.all = TRUE)
# peak 
peak_qrf_model <- qrf.peak$finalModel
class(peak_qrf_model) <- 'randomForest'
pred.peak <- predict(peak_qrf_model, test, predict.all = T)
# impulse 
impulse_qrf_model <- qrf.impulse$finalModel
class(impulse_qrf_model) <- 'randomForest'
pred.impulse <- predict(impulse_qrf_model, test, predict.all = T)
# tc
tc_qrf_model <- qrf.tc$finalModel
class(tc_qrf_model) <- 'randomForest'
pred.tc <- predict(tc_qrf_model, test, predict.all = T)

# plot density of trees in forest
plot_treelines <- function(rf_model, tree_pred, obs, pred, sub, ntrees, titlename, show.lm = F){
  if(titlename == 'Peak'){
    lims <- c(2.4,3.6) 
    ytit <- 'Observed [BW]'
    xtit <- 'Predicted [BW]'
    break_nums <- seq(2.4, 3.6, 0.4)
  } else if(titlename == 'tc'){
    lims <- c(0.15,0.23) 
    ytit <- 'Observed [s]'
    xtit <- 'Predicted [s]'
    break_nums <- seq(0.16, 0.22, 0.02)
  } else if(titlename == 'Impulse'){
    lims <- c(0.27,0.38)
    ytit <- 'Observed [BW-s]'
    xtit <- 'Predicted [BW-s]'
    break_nums <- seq(0.28, 0.38, 0.03)
  }
  tree_pred <- as.data.frame(tree_pred$individual)
  colnames(tree_pred) <- 1:ntrees
  tree_pred$obs <- obs
  tree_pred$sub <- sub
  
  tree_pred_long <- gather(tree_pred, tree, pred, 1:ntrees)
  agg <- aggregate(pred ~ obs + sub, data = tree_pred_long, FUN = mean)
  colnames(agg) <- c('obs', 'sub', 'mean_pred')
  
  pred_obs <- data.frame(pred = pred, obs = obs)
  pred_obs <- merge(pred_obs, agg, by = 'obs')
  p <- ggplot(tree_pred_long) + 
    geom_abline(slope = 1, intercept = 0, lty = 2)+
    geom_density_ridges(aes(x = pred, y = obs, group = obs, fill = as.factor(sub)),
                        alpha = 0.2, color = 'black', rel_min_height = 0.02, size = 0.2)+
    geom_point(data = pred_obs, aes(x = pred, y = obs, fill = sub), pch = 21, size = 2)
  
  if(show.lm){p <-p + geom_smooth(data = pred_obs, aes(x = pred, y = obs), 
                                  method = 'lm', se = F, color = 'black')}
  
  p <- p + 
    theme_classic()+
    ggtitle(paste(titlename, '- Quantile Regression Forest'))+
    coord_fixed(xlim = lims, ylim = lims)+
    scale_y_continuous(ytit, breaks = break_nums)+
    scale_x_continuous(xtit, breaks = break_nums)+
    scale_fill_tableau(palette = 'Classic 10')+
    theme(
      text = element_text(color = 'black', size = 11, face = 'plain'),
      axis.text = element_text(color = 'black', size = 10, face = 'plain'),
      axis.ticks = element_line(color = 'black'),
      legend.position = 'none',
      
    )
  
  print(p)
}
# Make predicted vs observed plots (with individual trees)
plot_treelines(peak_qrf_model, pred.peak, test.data$GRFPeak, test.data$qrf.peak, test.data$Sub,
               500, 'Peak', show.lm = F)  # forest_trees_tc.eps
plot_treelines(impulse_qrf_model, pred.impulse, test.data$GRFImpulse, test.data$qrf.impulse, test.data$Sub, 
               500, 'Impulse', show.lm = F)  # forest_trees_impulse.eps
plot_treelines(tc_qrf_model, pred.tc, test.data$GRFtc, test.data$qrf.tc, test.data$Sub, 
               500, 'tc', show.lm = F)  # forest_trees_tc.eps

# LR Plots ----
plot_MR_pred_obs <- function(df, x, y, titlename, show.lm){
  if(titlename == 'Peak'){
    lims <- c(2.4,3.6) 
    ytit <- 'Observed [BW]'
    xtit <- 'Predicted [BW]'
    break_nums <- seq(2.4, 3.6, 0.4)
  } else if(titlename == 'tc'){
    lims <- c(0.15,0.23) 
    ytit <- 'Observed [s]'
    xtit <- 'Predicted [s]'
    break_nums <- seq(0.16, 0.22, 0.02)
  } else if(titlename == 'Impulse'){
    lims <- c(0.27,0.38)
    ytit <- 'Observed [BW-s]'
    xtit <- 'Predicted [BW-s]'
    break_nums <- seq(0.28, 0.38, 0.03)
  }
  p <- ggplot(df, aes(x = x, y = y))+
    geom_abline(slope = 1, intercept = 0, lty = 2)+
    geom_point(aes(fill = Sub), color = 'black', pch = 23, size = 2)
  
  if(show.lm){p <-p + geom_smooth(data = df, aes(x = x, y = y), method = 'lm', se = F, color = 'black')}
  
  p <- p + 
    theme_classic()+
    ggtitle(paste(titlename, '- Quantile Regression Forest'))+
    coord_fixed(xlim = lims, ylim = lims)+
    scale_y_continuous(ytit, breaks = break_nums)+
    scale_x_continuous(xtit, breaks = break_nums)+
    scale_fill_tableau(palette = 'Classic 10')+
    theme(
      text = element_text(color = 'black', size = 11, face = 'plain'),
      axis.text = element_text(color = 'black', size = 10, face = 'plain'),
      axis.ticks = element_line(color = 'black'),
      legend.position = 'none',
      
    )
  
  print(p)
}
# make LR predicted vs observed plots
plot_MR_pred_obs(test.data, test.data$lm.peak, test.data$GRFPeak, 'Peak', show.lm = F)  # MR_peak.eps
plot_MR_pred_obs(test.data, test.data$lm.impulse, test.data$GRFImpulse, 'Impulse', show.lm = F)  # MR_impulse.eps
plot_MR_pred_obs(test.data, test.data$lm.tc, test.data$GRFtc, 'tc', show.lm = F)  # MR_tc.eps
