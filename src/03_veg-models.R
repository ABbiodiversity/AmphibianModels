#
# Title: ABMI Amphibian habitat models (North)
# Created: August 8th, 2018
# Last Updated: January 24th, 2022
# Author: Brandon Allen
# Objective: Species habitat models for the North analysis region
# Keywords: Veg models
# Notes: 
# 

##############
# Veg models # 
##############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
load("data/processed/site/amphibian-abundance-data-north-150m_2022-01-20.RData")
north.in <- veg.cur
rm(veg.cur)

# Impute the mean value for weather stations distance and temperature values with missing data
north.in$Weather_Station_Distance[is.na(north.in$Weather_Station_Distance)] <- mean(north.in$Weather_Station_Distance, na.rm = TRUE) 
north.in$Hourly_Temp[is.na(north.in$Hourly_Temp)] <- mean(north.in$Hourly_Temp, na.rm = TRUE) 

# Standardize the detection terms
north.in$ToY <- as.numeric(scale(north.in$ToY))
north.in$ToY2 <- as.numeric(scale(north.in$ToY2))
north.in$Hourly_Temp <- as.numeric(scale(north.in$Hourly_Temp))

#################
# All Recording #
#################

# Prediction matrix without water
pm <- read.csv("data/base/prediction-matrix/veg-prediction-matrix-CC_2021.csv")
rownames(pm) <- pm$VegType
pm <- pm[, -1]

# List containing models that will be looped through by the function

# Version 2022 Included in model approach
habitat.models <- list(as.formula(paste("pcount ~ WhiteSpruce + Pine + Deciduous + Mixedwood + BlackSpruce + TreedFen + TreedSwamp + GraminoidFen + ShrubbyFen + ShrubbyBog + ShrubbySwamp + Marsh + Grass + Shrub + CCWhiteSpruceR + CCWhiteSpruce1 + CCWhiteSpruce234 + CCPineR + CCPine1234 + CCDecidMixedR + CCDecidMixed1 + CCDecidMixed234 + HardLin + EnSeismic + EnSoftLin + TrSoftLin + UrbInd + Wellsites + Crop + RoughP + TameP + Rural")),
                       as.formula(paste("pcount ~ WhiteSpruce + Pine + Deciduous + Mixedwood + TreedBogFen + NonTreedBogFen + TreedSwamp + ShrubbySwamp + Grass + Shrub + CCWhiteSpruceR + CCWhiteSpruce1 + CCWhiteSpruce234 + CCPineR + CCPine1234 + CCDecidMixedR + CCDecidMixed1 + CCDecidMixed234 + HardLin + EnSeismic + EnSoftLin + TrSoftLin + UrbInd + Wellsites + Crop + RoughP + TameP + Rural")),
                       as.formula(paste("pcount ~ WhiteSpruce + Pine + Deciduous + Mixedwood + TreedBogFen + NonTreedBogFen + Swamp + GrassShrub + CCWhiteSprucePineR + CCWhiteSprucePine1 + CCWhiteSprucePine234 + CCDecidMixedR + CCDecidMixed1 + CCDecidMixed234 + HardLin + SoftLin + UrbInd + Wellsites + Crop + RoughP + TameP + Rural")),
                       as.formula(paste("pcount ~ WhiteSpruce + Pine + Deciduous + Mixedwood + Fen + Bog + MarshSwamp + GrassShrub + CCWhiteSprucePineR1234 + CCDecidMixedR1234 + HardLin + SoftLin + UrbIndWellsites + Crop + Pasture + Rural")),
                       as.formula(paste("pcount ~ WhiteSpruce + Pine + Deciduous + Mixedwood + Lowland + GrassShrub + CCWhiteSprucePineR1234 + CCDecidMixedR1234 + SoftLin + Alien")),
                       as.formula(paste("pcount ~ Upland + BlackSpruce + TreedFen + TreedSwamp + GraminoidFen + ShrubbyFen + ShrubbyBog + ShrubbySwamp + Marsh + Grass + Shrub + CCR1234 + HardLin + EnSeismic + EnSoftLin + TrSoftLin + UrbInd + Wellsites + Crop + RoughP + TameP + Rural")),
                       as.formula(paste("pcount ~ Upland + TreedBogFen + NonTreedBogFen + TreedSwamp + ShrubbySwamp + Grass + Shrub + CCR1234 + HardLin + EnSeismic + EnSoftLin + TrSoftLin + UrbInd + Wellsites + Crop + RoughP + TameP + Rural")),
                       as.formula(paste("pcount ~ Upland + TreedBogFen + NonTreedBogFen + Swamp + GrassShrub + CCR1234 + HardLin + SoftLin + UrbInd + Wellsites + Crop + RoughP + TameP + Rural")),
                       as.formula(paste("pcount ~ Upland + Fen + Bog + MarshSwamp + GrassShrub + CCR1234 + HardLin + SoftLin + UrbIndWellsites + Crop + Pasture + Rural")),
                       as.formula(paste("pcount ~ Upland + Lowland + GrassShrub + CCR1234 + SoftLin + Alien")))

# Add a combined clear clear group
north.in$CCR1234 <- rowSums(north.in[, c("CCWhiteSprucePineR1234", "CCDecidMixedR1234")])

landscape.names <- colnames(north.in)[c(41:129)]

climate.space.names <- c("Intercept", "AHM", "FFP", "MAP", "MAT",
                         "MCMT", "MWMT", "PET", "Lat", "Long",
                         "Lat2", "Long2", "LatLong", "MAT2", "MWMT2",
                         "MAPPET", "MAPFFP", "MATAHM", "Waterkm", "Waterkm2")

species.names <- c("BCFR", "CATO", "WETO", "WOFR")

coef.template <- list(matrix(ncol = length(landscape.names), nrow = 100, dimnames = list(c(1:100), c(landscape.names))),
                      matrix(ncol = length(landscape.names), nrow = 100, dimnames = list(c(1:100), c(landscape.names))),
                      matrix(ncol = 3, nrow = 100, dimnames = list(c(1:100), c("ToY", "ToY2", "Hourly_Temp"))),
                      matrix(ncol = length(climate.space.names), nrow = 100, dimnames = list(c(1:100), c(climate.space.names))),
                      matrix(ncol = 4, nrow = 100, dimnames = list(c(1:100), c("auc_LC", "auc_both", "dect", "survey"))))

names(coef.template) <- c("landscape.coef", "landscape.se", "detection.coef", "climate.coef", "fit")

####################
# Bootstrap Blocks # Uses Peters block. As we don't have many years of data, we don't block by year
####################

site.block <- data.frame(LongBlock = cut(north.in$Long, c(-121, -116, -112,-109)),
                         LatBlock = cut(north.in$Lat, c(48, 51, 54, 57, 61)))
site.block["Block"] <- interaction(droplevels(site.block$LongBlock), droplevels(site.block$LatBlock), sep="::", drop=TRUE)

# Reclassify
reclass.site <- data.frame(Orig = unique(site.block$Block),
                           Update = letters[1:length(unique(site.block$Block))])
north.in["Block"] <- reclass.site$Update[match(site.block$Block, reclass.site$Orig)]

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
        
        species.coef <- pa_models_3.0(data.analysis = north.in, results.store = species.coef,
                                      landscape.models = habitat.models, prediction.matrix = pm,
                                      region = "North", species.ID = species.names[spp],
                                      age.class = TRUE, weight.method = "IVW", detection = TRUE,
                                      coef.adjust = TRUE, boot = boot.set)
        
        print(boot.set)
        
    }
    
    species.boot.results[[species.names[spp]]] <- species.coef
    
}

# Save the bootstrap coefficients
save(species.boot.results, file = paste0("results/coef/veg-coefficients-amphibians-detection-IVW_", Sys.Date(), ".Rdata"))
