#
# Title: ABMI Amphibian habitat models (South)
# Created: August 8th, 2018
# Last Updated: January 24th, 2022
# Author: Brandon Allen
# Objective: Species habitat models for the South analysis region
# Keywords: Soil models
# Notes: 
# 

###############
# Soil models # 
###############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, lookups, and data set
library(arm)
library(binom)  # For exact binomial confidence intervals
library(mapproj)  # For projected maps
library(mgcv)  # For binomial GAM
library(MuMIn)
library(pROC)
library(RcmdrMisc)

source("src/habitat-models-functions_2020-10-28.R")

# Landscape, climate, space
load("data/processed/site/amphibian-abundance-data-south-150m_2022-01-20.RData")
south.in <- soil.cur
rm(soil.cur)

# Impute the mean value for weather stations distance and temperature values with missing data
south.in$Weather_Station_Distance[is.na(south.in$Weather_Station_Distance)] <- mean(south.in$Weather_Station_Distance, na.rm = TRUE) 
south.in$Hourly_Temp[is.na(south.in$Hourly_Temp)] <- mean(south.in$Hourly_Temp, na.rm = TRUE) 

# Standardize the detection terms
south.in$ToY <- as.numeric(scale(south.in$ToY))
south.in$ToY2 <- as.numeric(scale(south.in$ToY2))
south.in$Hourly_Temp <- as.numeric(scale(south.in$Hourly_Temp))

#################
# All Recording #
#################

# Prediction matrix No water
pm <- read.csv("data/base/prediction-matrix/soil-prediction-matrix_2021.csv")
rownames(pm) <- pm$VegType
pm <- pm[, -1]

# List containing models that will be looped through by the function No water
habitat.models <- list(as.formula(paste("pcount ~ Blowout + ClaySub + Loamy + RapidDrain + SandyLoam + ThinBreak + Other + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                       as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                       as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                       as.formula(paste("pcount ~ Blowout + ClaySub + Loamy + RapidDrain + SandyLoam + ThinBreak  + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Cult + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + Alien + SoftLin")),
                       as.formula(paste("pcount ~ Blowout + ClaySub + Loamy + RapidDrain + SandyLoam + ThinBreak  + Other + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid +  Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid +  Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbInd + Wellsites + Rural + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Blowout + ClaySub + Loamy +RapidDrain + SandyLoam + ThinBreak  + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Cult + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Productive + Nonproductive + Other + Alien + SoftLin + paspen")))

# Water excluded
landscape.names <- c("Loamy", "SandyLoam", "ClaySub", "RapidDrain",
                     "Blowout", "ThinBreak", "Other", "EnSeismic",
                     "EnSoftLin", "TrSoftLin", "HardLin", "UrbInd", "Rural",
                     "Wellsites", "Crop", "TameP", "RoughP")

climate.space.names <- c("Intercept", "AHM", "FFP", "MAP", "MAT", 
                         "MCMT", "MWMT", "PET", "Lat", "Long", 
                         "Lat2", "Long2", "LatLong", "MAT2", "MWMT2", 
                         "MAPPET", "MAPFFP", "MATAHM", "Waterkm", "Waterkm2")

species.names <- c("BCFR", "CATO", "WETO", "WOFR")

coef.template <- list(matrix(ncol = length(landscape.names), nrow = 100, dimnames = list(c(1:100), c(landscape.names))), 
                     matrix(ncol = length(landscape.names), nrow = 100, dimnames = list(c(1:100), c(landscape.names))), 
                     matrix(ncol = 3, nrow = 100, dimnames = list(c(1:100), c("ToY", "ToY2", "Hourly_Temp"))),
                     matrix(ncol = length(climate.space.names), nrow = 100, dimnames = list(c(1:100), c(climate.space.names))), 
                     matrix(ncol = 1, nrow = 100, dimnames = list(c(1:100), c("paspen"))), 
                     matrix(ncol = 1, nrow = 100, dimnames = list(c(1:100), c("paspen"))),
                     matrix(ncol = 4, nrow = 100, dimnames = list(c(1:100), c("auc_LC", "auc_both", "dect", "survey"))))

names(coef.template) <- c("landscape.coef", "landscape.se", "detection.coef", "climate.coef", "paspen.coef", "paspen.se", "fit")

####################
# Bootstrap Blocks # Uses Peters block. As we don't have many years of data, we don't block by year
####################

site.block <- data.frame(LongBlock = cut(south.in$Long, c(-121, -116, -112,-109)),
                         LatBlock = cut(south.in$Lat, c(48, 51, 54, 57, 61)))
site.block["Block"] <- interaction(droplevels(site.block$LongBlock), droplevels(site.block$LatBlock), sep="::", drop=TRUE)

# Reclassify
reclass.site <- data.frame(Orig = unique(site.block$Block),
                           Update = letters[1:length(unique(site.block$Block))])
south.in["Block"] <- reclass.site$Update[match(site.block$Block, reclass.site$Orig)]

rm(site.block, reclass.site)

####################
# Species modeling # 
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################
# Presence/Absence #
####################

species.boot.results <- NULL

for (spp in 1:length(species.names)) {
    
    species.coef <- coef.template
    
    for (boot.set in 1:100) {
        
        species.coef <- pa_models_3.0(data.analysis = south.in, results.store = species.coef,
                                      landscape.models = habitat.models, prediction.matrix = pm,
                                      region = "South", species.ID = species.names[spp],
                                      age.class = FALSE, weight.method = "IVW", detection = TRUE,
                                      coef.adjust = FALSE, boot = boot.set)
        
        print(boot.set)
        
    }
    
    species.boot.results[[species.names[spp]]] <- species.coef
    
}

# Save the bootstrap coefficients
save(species.boot.results, file = paste("results/coef/soil-coefficients-amphibians-detection-IVW_", Sys.Date(), ".Rdata", sep = ""))
