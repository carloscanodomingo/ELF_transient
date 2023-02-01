library(easyGgplot2)
library(latex2exp)
library(R.matlab)
library(ggplot2)
library(qqplotr)
library(qualityTools)
library(goft)
library(fitdistrplus)


matlabFile <- readMat("r_am_sr.mat")




varNames    <- names(matlabFile$data[,,1])
datList     <- matlabFile$data
datList     <- lapply(datList, unlist, use.names=FALSE)
data        <- as.data.frame(datList)
names(data) <- varNames
# Histogram from a single numeric vector 
# ggplot2.histogram(data=numVector)
# Basic histogram plot from the vector "weight"
weight = c(data$FWHMT);

data$Type[data$peak.value < 2e6] <- "Normal";
data$Type[data$peak.value >= 2e6] <- "Qburst";

fitl<-fitdist(data$peak.value[data$Type != "Qburst"]/1e6,"lnorm");
summary(fitl)
mu = fitl$estimate[1];
sigma = fitl$estimate[2];
mean_fit = exp(mu + 1/2* sigma^2);mean_fit
var_fit =  exp(2*mu + sigma^2)*(exp(sigma^2) - 1); std_fit = sqrt(var_fit);std_fit

fitl<-fitdist(data$T1[data$Type != "Qburst"],"lnorm");
summary(fitl)
mu = fitl$estimate[1];
sigma = fitl$estimate[2];
mean_fit = exp(mu + 1/2* sigma^2);mean_fit
var_fit =  exp(2*mu + sigma^2)*(exp(sigma^2) - 1); std_fit = sqrt(var_fit);std_fit

fitl<-fitdist(data$PSD.in.band[data$Type != "Qburst" & data$PSD.in.band != 0]/1e10,"lnorm");
summary(fitl)
mu = fitl$estimate[1];
sigma = fitl$estimate[2];
mean_fit = exp(mu + 1/2* sigma^2);mean_fit
var_fit =  exp(2*mu + sigma^2)*(exp(sigma^2) - 1); std_fit = sqrt(var_fit)  ;std_fit


ks.test(data$peak.value[data$Type != "Qburst"]/1e6,"plnorm", meanlog=coef(fitg)[1], sdlog =coef(fitg)[2] )
