library(ensembleBMA)
rm(list=ls())

?ensembleBMAnormal

data(srft)
memberLabels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")
srftData <- ensembleData(forecasts = srft[,memberLabels],
                         dates = srft$date, observations = srft$obs,
                         latitude = srft$lat, longitude = srft$lon,
                         forecastHour = 48, initializationTime = "00")

?julTOymdh
JulianDates=ymdhTOjul(srftData$date)
head(JulianDates)
head(as.character(srftData$date))
head(as.character(julTOymdh(JulianDates-1461)))

exGroups <- c( CMCG=1, ETA=2, GASP=3, GFS=2, JMA=4, NGPS=5, TCWB=6, UKMO=7)
exGroups <- c( 1:3,2,4:7)
exGroups <- rep(1,8)
srftDataX <- ensembleData(forecasts = srft[,memberLabels],
                          dates = srft$date, observations = srft$obs,
                          latitude = srft$lat, longitude = srft$lon,
                          forecastHour = 48, initializationTime = "00",
                          exchangeable = exGroups)

srftFit <- ensembleBMA( srftDataX, dates = "2004013100",
                        model = "normal", trainingDays = 25)
plot( srftFit, srftDataX, dates = "2004013100")

train <- trainingData( srftData, date = "2004013100",
                       trainingDays = 25)
srftTrainFit <- fitBMA( train,  model = "normal")

srftForc <- quantileForecast( srftFit, srftData,
                              quantiles = c( .1, .5, .9))

use <- as.character(srftData$dates) == "2004013100"
lat <- srftData$latitude[use]; lon <- srftData$longitude[use]
lonRange <- range(lon); latRange <- range(lat)
range(srftForc[,"0.5"])
color <- "brown"; mapColor <- "black"

library(fields)
library(maps)
plotProbcast( srftForc[,"0.5"], lon, lat, interpolate = TRUE, col = color,
              type = "contour", levels = seq(from=264, to=284, by=2))
title("Median Forecast")
points(lon, lat, pch = 16, cex = 0.5, col = color)  # observation locations

data(ensBMAtest)

ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

obs <- paste("T2","obs", sep = ".")
ens <- paste("T2", ensMemNames, sep = ".")


tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                              dates = ensBMAtest[,"vdate"],
                              observations = ensBMAtest[,obs],
                              station = ensBMAtest[,"station"],
                              forecastHour = 48,
                              initializationTime = "00")

## Not run:  #  R check
tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
MAE( tempTestFit, tempTestData)
CRPS( tempTestFit, tempTestData)


##########

setwd("~/Documents/courbariaux/WORKLUC/TempEnsemblesNormal")
rm(list=ls())
library(ensembleBMA)
load("Ain@Vouglans.Rdata")
#load("Buech@Chambons.Rdata")
#load("Drac@Sautet.Rdata")
library(tidyverse)
library(tidyselect)
library(lubridate)
library(R2jags)
library(verification)
library(GGally)
dtout %>% filter(Echeance==4, Date <"2005-03-31") %>% 
  dplyr::select(-Year, -Month, -Calendaire, -Echeance) %>% 
  mutate(Date=(gsub("-",'',as.character(Date))))->d
memberLabels <- paste0("Run",1:50)
dData <- ensembleData(forecasts = d[,memberLabels],
                         dates = d$Date, observations = d$Obs,
                         forecastHour = 96, initializationTime = "00")

dFit <- ensembleBMA( dData, dates = c("20050205","20050206"),
                        model = "normal", trainingDays = 25)
plot( dFit, dData, dates = c("20050205","20050206"))
modelParameters(dFit)->md
dData[34,51:52]
w<-md$weights[,2]
coef<-md$biasCoefs[,,2]
std<-md$sd[2]
runs=dData[34,1:50]
centres= coef[1,]+coef[2,]*runs
repet=1000
qui=sample(x = 1:50,size=repet,replace = T,prob=w)
samp=rnorm(n = repet,mean = as.numeric(centres[qui]),sd=rep(std,repet))
hist(samp,nc=50,freq=F,add=T,fill=T, col='lightgrey')
exch=rep(1,50)
colnames(exch)<-as.character(memberLabels)
dDataX <- ensembleData(forecasts = d[,memberLabels],
                      dates = d$Date, observations = d$Obs,
                      forecastHour = 96, initializationTime = "00",
                      exchangeable = c(Run1=1,rep(1,49)))
dFitX <- ensembleBMA( dDataX, dates = c("20050205","20050206"),
                     model = "normal", trainingDays = 25,
                     exchangeable = c(Run1=1,rep(1,49)))
plot( dFitX, dDataX, dates = c("20050205","20050206"))