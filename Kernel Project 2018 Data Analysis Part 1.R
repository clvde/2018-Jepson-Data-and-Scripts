### Jepson Prairie Kernel Project May 2018
### Summary Statistical Analysis of Neighbour Plants

## General Experimental Design
# Follows closely the fitness dispersal project. 
# For each species, 10 sets of 2 maternal plants chosen randomly along a 2D transect grid
# Transects were run first along the longest axis of the pool or patch, then along the perpendicular axis. 
# Random 2D point chosen, nearest plant that hadn't dropped seeds yet was chosen.
# Three species total = 20*3 = 60 maternal plants
# surrounding each maternal plant, 5 neighbours chosen. Total of 300 neighbour plants in theory (actually fewer)
# For each neighbour:
# (1) Reproductive biomass, 
# (2) Inflorescense diameter (For single-flower plants)
# (3) 

## Packages
library(tidyverse)

## Set Working Directory
setwd("/Users/Courtney/Documents/2018-Jepson-Data-and-Scripts")

## Read and clean neighbour data
neighb_data <- read.csv("./Field Kernel Data 2018 - Neighbour Data.csv", header = TRUE, na.strings="NA")

# Fix data classes
neighb_data$Number <- as.factor(neighb_data$Number)
neighb_data$Mom <- as.factor(neighb_data$Mom)
neighb_data$Neighbour <- as.factor(neighb_data$Neighbour)
neighb_data$Tube <- as.factor(neighb_data$Tube)
neighb_data$Date_Collected <- as.Date(neighb_data$Date_Collected, "%m/%d/%Y")  
neighb_data$Date_Measured <- as.Date(neighb_data$Date_Measured, "%m/%d/%Y") 
neighb_data$Viable_Ray_Seed_Count <- as.numeric(neighb_data$Viable_Ray_Seed_Count)
neighb_data$Undev_Ray_Seed_Count <- as.numeric(neighb_data$Undev_Ray_Seed_Count)

# Delete empty columns until data entry is complete.
neighb_data$Inflor_Height <- NULL
neighb_data$Branch_Lengths <- NULL
neighb_data$Vegetative_Biomass <- NULL



## Read and clean maternal plant data
mom_data <- read.csv("./Field Kernel Data 2018 - Maternal Plant Data.csv", header = TRUE, na.strings="NA")
mom_data_filt <- filter(mom_data, !is.na(Measure_Date))

## Question 1: How many seeds per inflorescence by size and species?

# Filter out NA rows for Diameter measurements

neighb_data_1flr <- filter(neighb_data, !is.na(Diameter))

# Species-specific data frames
neighb_data_1flr_cal <- filter(neighb_data_1flr, Species == "cal")
neighb_data_1flr_fre <- filter(neighb_data_1flr, Species == "fre")
neighb_data_1flr_gla <- filter(neighb_data_1flr, Species == "gla")

# Plot of viable disc seeds by floral diameter
ggplot(neighb_data_1flr_cal, aes(Diameter, Viable_Disc_Seed_Count)) + 
  geom_point() +
  geom_smooth (method = "lm", formula = y ~ exp(x)) +
  geom_smooth (method = "lm", formula = y ~ x + I(x^2)) +
  geom_smooth (method = "lm", formula = y ~ x) +
  ggtitle("L. californica")

ggplot(neighb_data_1flr_fre, aes(Diameter, Viable_Disc_Seed_Count)) + 
  geom_point() +
  #geom_smooth (method = "lm", formula = y ~ exp(x)) +
  geom_smooth (method = "lm", formula = y ~ x + I(x^2)) +
  #geom_smooth (method = "lm", formula = y ~ x) +
  ggtitle("L. fremontii")

ggplot(neighb_data_1flr_gla, aes(Diameter, Viable_Disc_Seed_Count)) + 
  geom_point() +
  geom_smooth (method = "lm", formula = y ~ exp(x)) +
  geom_smooth (method = "lm", formula = y ~ x + I(x^2)) +
  geom_smooth (method = "lm", formula = y ~ x) +
  ggtitle("L. glaberrima")

## PREDICTION ERROR RATE
## Measuring prediction error rate between the two models - drop 1 cross validation for cal data

# Number of observations total
n_obs <- dim(neighb_data_1flr_cal)[1]

# Number of observations to drop
ndrop <- 1

# Initialize matrices to store predicted data
drop_pred.lm <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.lm.poly <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.lm.exp <- matrix(nrow = n_obs, ncol = ndrop)

# Initialize matrix to store observed data
drop_obs <- matrix(nrow = n_obs, ncol = ndrop)

# Iterate through all data points dropping ndrop at a time, then fit model and predict dropped pts
for (i in 1:n_obs){
  if (ndrop == 1){ # drop one data point at a time
    rand_drop <- i
  
    # training dataset
    train_set <- neighb_data_1flr_cal[-rand_drop,]
    
    # testing dataset
    drop_set <- neighb_data_1flr_cal[rand_drop,]
    
    # fit a linear model, a quadratic model, and an exponential model
    cv.lm <- lm(Viable_Disc_Seed_Count ~ Diameter, data=train_set)
    cv.lm.poly <- lm(Viable_Disc_Seed_Count ~ Diameter+I(Diameter^2), data=train_set)
    cv.lm.exp <- lm(Viable_Disc_Seed_Count ~ exp(Diameter), data=train_set)
    
    # predict the dropped values and store them
    dropdf <- data.frame(Diameter=drop_set$Diameter)
    drop_pred.lm[i,] <- predict(cv.lm, dropdf)
    drop_pred.lm.poly[i,] <- predict(cv.lm.poly, dropdf)
    drop_pred.lm.exp[i,] <- predict(cv.lm.exp, dropdf)
    drop_obs[i,] <- drop_set$Viable_Disc_Seed_Count
    
  } else { # more than one point dropped at a time
    rand_drop <- round(runif(ndrop, min=1, max=n_obs))
    
    # all remaining data points are training data set
    train_set <- neighb_data_1flr_cal[-rand_drop,]
    drop_set <- neighb_data_1flr_cal[rand_drop,]
    
    # fit a linear model
    cv.lm <- lm(Viable_Disc_Seed_Count~Diameter,data=train_set)
    cv.lm.poly <- lm(Viable_Disc_Seed_Count~Diameter + I(Diameter^2), data=train_set)
    cv.lm.exp <- lm(Viable_Disc_Seed_Count~exp(Diameter), data=train_set)
    
    # predict the dropped values and store them
    dropdf <- data.frame(Diameter=drop_set$Diameter)
    drop_pred.lm[i,] <- predict(cv.lm, dropdf)
    drop_pred.lm.poly[i,] <- predict(cv.lm.poly, dropdf)
    drop_pred.lm.exp[i,] <- predict(cv.lm.exp, dropdf)
    drop_obs[i,] <- drop_set$Viable_Disc_Seed_Count
  }
}

# Calculate the MSE (CV score) for each model
cal_lm_err <- sum((drop_obs-drop_pred.lm)^2)/n_obs
cal_polylm_err <- sum((drop_obs-drop_pred.lm.poly)^2)/n_obs
cal_explm_err <- sum((drop_obs-drop_pred.lm.exp)^2)/n_obs

# Lowest score means lowest prediction error. Linear model works best here.
cal_lm_err
cal_polylm_err
cal_explm_err


## Drop 1 cross validation for fre data
n_obs <- dim(neighb_data_1flr_fre)[1]
ndrop <- 1

drop_pred.lm <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.lm.poly <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.lm.exp <- matrix(nrow = n_obs, ncol = ndrop)

drop_obs <- matrix(nrow = n_obs, ncol = ndrop)

for (i in 1:n_obs){
  if (ndrop == 1){
    # randomly choose one data point to drop
    rand_drop <- i
    
    # all remaining data points are training data set
    train_set <- neighb_data_1flr_fre[-rand_drop,]
    
    drop_set <- neighb_data_1flr_fre[rand_drop,]
    
    # fit a linear model
    cv.lm <- lm(Viable_Disc_Seed_Count ~ Diameter, data=train_set)
    cv.lm.poly <- lm(Viable_Disc_Seed_Count ~ Diameter+I(Diameter^2), data=train_set)
    cv.lm.exp <- lm(Viable_Disc_Seed_Count ~ exp(Diameter), data=train_set)
    
    # predict the dropped values and store them
    dropdf <- data.frame(Diameter=drop_set$Diameter)
    drop_pred.lm[i,] <- predict(cv.lm, dropdf)
    drop_pred.lm.poly[i,] <- predict(cv.lm.poly, dropdf)
    drop_pred.lm.exp[i,] <- predict(cv.lm.exp, dropdf)
    drop_obs[i,] <- drop_set$Viable_Disc_Seed_Count
    
  } else {
    # randomly choose one data point to drop
    rand_drop <- round(runif(ndrop, min=1, max=n_obs))
    
    # all remaining data points are training data set
    train_set <- neighb_data_1flr_fre[-rand_drop,]
    drop_set <- neighb_data_1flr_fre[rand_drop,]
    
    # fit a linear model
    cv.lm <- lm(Viable_Disc_Seed_Count~Diameter,data=train_set)
    cv.lm.poly <- lm(Viable_Disc_Seed_Count~Diameter + I(Diameter^2), data=train_set)
    cv.lm.exp <- lm(Viable_Disc_Seed_Count~exp(Diameter), data=train_set)
    
    # predict the dropped values and store them
    dropdf <- data.frame(Diameter=drop_set$Diameter)
    drop_pred.lm[i,] <- predict(cv.lm, dropdf)
    drop_pred.lm.poly[i,] <- predict(cv.lm.poly, dropdf)
    drop_pred.lm.exp[i,] <- predict(cv.lm.exp, dropdf)
    drop_obs[i,] <- drop_set$Viable_Disc_Seed_Count
  }
}

fre_lm_err <- sum((drop_obs-drop_pred.lm)^2)/n_obs
fre_polylm_err <- sum((drop_obs-drop_pred.lm.poly)^2)/n_obs
fre_explm_err <- sum((drop_obs-drop_pred.lm.exp)^2)/n_obs

# Quadratic model is best
fre_lm_err
fre_polylm_err
fre_explm_err

## Drop 1 cross validation for gla data
n_obs <- dim(neighb_data_1flr_gla)[1]
ndrop <- 1

drop_pred.lm <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.lm.poly <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.lm.exp <- matrix(nrow = n_obs, ncol = ndrop)

drop_obs <- matrix(nrow = n_obs, ncol = ndrop)

for (i in 1:n_obs){
  if (ndrop == 1){
    # randomly choose one data point to drop
    rand_drop <- i
    
    # all remaining data points are training data set
    train_set <- neighb_data_1flr_gla[-rand_drop,]
    drop_set <- neighb_data_1flr_gla[rand_drop,]
    
    # fit a linear model
    cv.lm <- lm(Viable_Disc_Seed_Count ~ Diameter, data=train_set)
    cv.lm.poly <- lm(Viable_Disc_Seed_Count ~ Diameter+I(Diameter^2), data=train_set)
    cv.lm.exp <- lm(Viable_Disc_Seed_Count ~ exp(Diameter), data=train_set)
    
    # predict the dropped values and store them
    dropdf <- data.frame(Diameter=drop_set$Diameter)
    drop_pred.lm[i,] <- predict(cv.lm, dropdf)
    drop_pred.lm.poly[i,] <- predict(cv.lm.poly, dropdf)
    drop_pred.lm.exp[i,] <- predict(cv.lm.exp, dropdf)
    drop_obs[i,] <- drop_set$Viable_Disc_Seed_Count
    
  } else {
    # randomly choose one data point to drop
    rand_drop <- round(runif(ndrop, min=1, max=n_obs))
    
    # all remaining data points are training data set
    train_set <- neighb_data_1flr_gla[-rand_drop,]
    drop_set <- neighb_data_1flr_gla[rand_drop,]
    
    # fit a linear model
    cv.lm <- lm(Viable_Disc_Seed_Count~Diameter,data=train_set)
    cv.lm.poly <- lm(Viable_Disc_Seed_Count~Diameter + I(Diameter^2), data=train_set)
    cv.lm.exp <- lm(Viable_Disc_Seed_Count~exp(Diameter), data=train_set)
    
    # predict the dropped values and store them
    dropdf <- data.frame(Diameter=drop_set$Diameter)
    drop_pred.lm[i,] <- predict(cv.lm, dropdf)
    drop_pred.lm.poly[i,] <- predict(cv.lm.poly, dropdf)
    drop_pred.lm.exp[i,] <- predict(cv.lm.exp, dropdf)
    drop_obs[i,] <- drop_set$Viable_Disc_Seed_Count
  }
}

gla_lm_err <- sum((drop_obs-drop_pred.lm)^2)/n_obs
gla_polylm_err <- sum((drop_obs-drop_pred.lm.poly)^2)/n_obs
gla_explm_err <- sum((drop_obs-drop_pred.lm.exp)^2)/n_obs

# Linear model is best
gla_lm_err
gla_polylm_err
gla_explm_err

## PREDICTING MATERNAL PLANT DATA
## Use the winning model for each species to predict number of seeds expected for each maternal plant

# Get diameters for maternal plants to predict number of seeds.
mom_cal <- filter(mom_data_filt, Species == "cal")
mom_fre <- filter(mom_data_filt, Species == "fre")
mom_gla <- filter(mom_data_filt, Species == "gla")

# Fit linear models for each species
cal.lm = lm(Viable_Disc_Seed_Count ~ Diameter, neighb_data_1flr_cal)
fre.lm = lm(Viable_Disc_Seed_Count ~ Diameter, neighb_data_1flr_fre)
gla.lm = lm(Viable_Disc_Seed_Count ~ Diameter, neighb_data_1flr_gla)
cal.lm.poly = lm(Viable_Disc_Seed_Count~Diameter + I(Diameter^2), neighb_data_1flr_cal)
fre.lm.poly = lm(Viable_Disc_Seed_Count~Diameter + I(Diameter^2), neighb_data_1flr_fre)
gla.lm.poly = lm(Viable_Disc_Seed_Count~Diameter + I(Diameter^2), neighb_data_1flr_gla)

# californica - linear model

# predict values
mom_diams_cal <- data.frame(Diameter = mom_cal$Flr_Size_mm)
mom_disc_pred_cal <- as.data.frame(predict(cal.lm, newdata = mom_diams_cal, interval="prediction"))

# add predicted values and standard errors to the data frame
mom_cal_wpred <- mom_cal
mom_cal_wpred$Predicted_Disc_Seed_Count <- mom_disc_pred_cal$fit
mom_cal_wpred$Predicted_Disc_Seed_PI_LWR <- mom_disc_pred_cal$lwr
mom_cal_wpred$Predicted_Disc_Seed_PI_UPR <- mom_disc_pred_cal$upr
mom_cal_wpred

# fremontii - quadratic model

mom_diams_fre <- data.frame(Diameter = mom_fre$Flr_Size_mm)
mom_disc_pred_fre <- as.data.frame(predict(fre.lm.poly, newdata = mom_diams_fre, interval="prediction"))

mom_fre_wpred <- mom_fre
mom_fre_wpred$Predicted_Disc_Seed_Count <- mom_disc_pred_fre$fit
mom_fre_wpred$Predicted_Disc_Seed_PI_LWR <- mom_disc_pred_fre$lwr
mom_fre_wpred$Predicted_Disc_Seed_PI_UPR <- mom_disc_pred_fre$upr
mom_fre_wpred

# glaberrima - linear model

mom_diams_gla <- data.frame(Diameter = mom_gla$Flr_Size_mm)
mom_disc_pred_gla <-  as.data.frame(predict(gla.lm, newdata = mom_diams_gla, interval="prediction"))

mom_gla_wpred <- mom_gla
mom_gla_wpred$Predicted_Disc_Seed_Count <- mom_disc_pred_gla$fit
mom_gla_wpred$Predicted_Disc_Seed_PI_LWR <- mom_disc_pred_fre$lwr
mom_gla_wpred$Predicted_Disc_Seed_PI_UPR <- mom_disc_pred_fre$upr

mom_gla_wpred

# Create modified dataframe with all three species

mom_wpred <- arrange(rbind(mom_cal_wpred,mom_fre_wpred,mom_gla_wpred), Hub)
write.csv(mom_wpred, file = "./maternal_data_with_predictions.csv")

# Create mod dataframe with all three species but BY WHOLE PLANT - ONE ROW PER PLANT INSTEAD OF PER FLOWER

# Reduce so one row per maternal plant
mom_cal_filt <- filter(mom_cal, Flr_ID == 1)
mom_cal_filt$Flr_Size_mm <- NULL

# Sum all inflorescences to find is a total predicted number of seeds for the entire plant
sum_predcount <- as.data.frame(summarize(group_by(mom_cal_wpred,Full_ID), Predicted_Disc_Seed_Count=sum(Predicted_Disc_Seed_Count)))

# Merge the two so we have predicted and found seed counts for each plant
mom_cal_filt2 <- merge(mom_cal_filt, sum_predcount,by="Full_ID")

# Calculate the 'percentage found' for each plant
mom_cal_filt2$Proportion_Found <- mom_cal_filt2$N_Found_Seeds/mom_cal_filt2$Predicted_Disc_Seed_Count

mom_cal_filt2
