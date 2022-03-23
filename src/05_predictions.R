#
# Title: Prediction calculation script
# Created: February 12th, 2020
# Last Updated: February 1st, 2022
# Author: Brandon Allen
# Objectives: Calculate provincial predictions for each validated model
# Keywords: 
# Notes: 1) The landcover coefficients need to be backtransformed binomial()$linkfun in order to proper predictions

#####################
# Transition Matrix #
#####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
rm(list=ls())
gc()

# load libraries and functions
library(mefa4)
library(openxlsx)

# Define Current year
cur.year <- 2018

# Load the Kgrid and transition matrix
load("data/base/kgrid/kgrid_table_km.RData")
load(paste0("data/base/kgrid/veghf_w2w_ref_", cur.year, "_transitions_wide_water.RData"))

# if cur.year = 2010, update names

if (cur.year == 2010) {
        
        # There is no sector effect information for the soil lookup in 2010, match results from 2018 to correct
        chSoil2018 <- read.csv("data/lookup/transition/soil-transition_ref-2018.csv")
        chSoil <- chSoil10
        chSoil$sector <- chSoil2018$sector[match(chSoil$label, chSoil2018$label)]
        chSoil["LenT->CultivationRoughPasture", "sector"] <- "Agriculture"
        chSoil$sector <- as.factor(chSoil$sector)
        rm(chSoil2018)
        
        # Everything else is okay
        chVeg <- chVeg10
        trSoil <- trSoil10
        trVeg <- trVeg10
        
        rm(chSoil10, chVeg10, trSoil10, trVeg10)
        
}

# Align the kgrid with the transition matrix
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]

#########################################
# Climate and Weighting Standardization #
#########################################

# Organize the space and climate information for the north and south
make_clim <- function(x, birds=FALSE, truncate_latitude=FALSE) {
        if (birds) {
                z <- with(x, cbind(
                        pWater_KM=pWater,
                        pWater2_KM=pWater^2,
                        xPET=(PET - 0) / 800,
                        xMAT=(MAT - 0) / 6,
                        xAHM=(AHM - 0) / 50,
                        xFFP=(FFP - 0) / 130,
                        xMAP=(MAP - 0) / 2300,
                        xMWMT=(MWMT - 0) / 20,
                        xMCMT=(MCMT - 0) / 25,
                        xY=(POINT_Y - 54.1) / 2.43,
                        xX=(POINT_X - (-113.3)) / 2.12))
                z <- cbind(z,
                           xY2=z[,"xY"]^2,
                           xX2=z[,"xX"]^2,
                           `xFFP:xMAP`=z[,"xFFP"]*z[,"xMAP"],
                           `xMAP:xPET`=z[,"xMAP"]*z[,"xPET"],
                           `xAHM:xMAT`=z[,"xAHM"]*z[,"xMAT"],
                           `xX:xY`=z[,"xX"]*z[,"xY"])
        } else {
                LAT <- pmin(x$POINT_Y, 56.5)
                z <- with(x, cbind(
                        Intercept=1,
                        Lat=LAT,
                        Long=x$POINT_X,
                        AHM=x$AHM,
                        PET=x$PET,
                        FFP=x$FFP,
                        MAP=x$MAP,
                        MAT=x$MAT,
                        MCMT=x$MCMT,
                        MWMT=x$MWMT,
                        Lat2=LAT^2,
                        Long2=x$POINT_X^2,
                        LatLong=x$POINT_X*LAT,
                        MAPPET=x$MAP*x$PET,
                        MATAHM=x$MAT*x$AHM,
                        MAPFFP=x$MAP*x$FFP,
                        MAT2=x$MAT^2,
                        MWMT2=x$MWMT^2,
                        Waterkm=x$pWater,
                        Waterkm2=x$pWater^2))
        }
        rownames(z) <- rownames(x)
        z
}

## making space-climate model matrix
Xclim_nonb <- make_clim(kgrid, birds=FALSE)

kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

UseN <- rownames(kgrid)[!kgrid$useS]
UseS <- rownames(kgrid)[!kgrid$useN]

Xclim_nonb_S <- Xclim_nonb[UseS,]
pA <- kgrid[UseS, "pAspen"]
Xclim_nonb_N <- Xclim_nonb[UseN,]

# Remove unecessary information
rm(Xclim_nonb, make_clim, UseS, UseN)
gc()

FORE <- c("CCDecid1", "CCDecid2",
          "CCDecid3", "CCDecid4", "CCDecidR", "CCMixedwood1", "CCMixedwood2",
          "CCMixedwood3", "CCMixedwood4", "CCMixedwoodR", "CCPine1", "CCPine2",
          "CCPine3", "CCPine4", "CCPineR", "CCSpruce1", "CCSpruce2", "CCSpruce3",
          "CCSpruce4", "CCSpruceR",
          "Decid1", "Decid2", "Decid3", "Decid4", "Decid5", "Decid6", "Decid7", "Decid8",
          "DecidR", "Mixedwood1", "Mixedwood2",
          "Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7",
          "Mixedwood8", "MixedwoodR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6",
          "Pine7", "Pine8", "PineR",
          "Spruce1", "Spruce2", "Spruce3", "Spruce4", "Spruce5", "Spruce6",
          "Spruce7", "Spruce8", "SpruceR", "TreedBog1",
          "TreedBog2", "TreedBog3", "TreedBog4", "TreedBog5", "TreedBog6",
          "TreedBog7", "TreedBog8", "TreedBogR", "TreedFen1", "TreedFen2",
          "TreedFen3", "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7",
          "TreedFen8", "TreedFenR", "TreedSwamp1", "TreedSwamp2", "TreedSwamp3",
          "TreedSwamp4", "TreedSwamp5", "TreedSwamp6", "TreedSwamp7", "TreedSwamp8",
          "TreedSwampR")

UNKN <- chSoil$rf=="UNK" | chSoil$cr=="UNK"
kgrid$pSoilUnk <- rowSums(trSoil[,UNKN,drop=FALSE]) / rowSums(trSoil)
kgrid$pFor <- rowSums(trVeg[,chVeg$cr %in% FORE]) / rowSums(trVeg)

kgrid$wS <- (1-kgrid$pAspen)*(1-kgrid$pSoilUnk)*(1-kgrid$pFor)
kgrid$wS[kgrid$useN] <- 0
kgrid$wS[kgrid$useS] <- 1
kgrid$wN <- kgrid$pAspen / (kgrid$pAspen + (1-kgrid$pFor))
kgrid$wN[kgrid$useN] <- 1
kgrid$wN[kgrid$useS] <- 0
wsm <- kgrid$wS + kgrid$wN
kgrid$wS <- kgrid$wS / wsm
kgrid$wN <- kgrid$wN / wsm

kgrid <- kgrid[, c("Row_Col", "Row", "Col", "NSRNAME", "NRNAME", "LUF_NAME", "useN", "useS", "pAspen",
                   "pSoilUnk", "pFor", "wS", "wN")]

# Clear memory
rm(FORE, pA, UNKN, wsm)
gc()

# Create lookup table for matching coefficients
lt <- list(
        south=structure(list(Label = c("ClaySub", "Other", "Other", "Other",
                                       "Other", "Other", "RapidDrain", "Loamy", "Other", "Other", "Other",
                                       "Other", "Other", "Other", "ClaySub", "SandyLoam", "RapidDrain",
                                       "RapidDrain", "RapidDrain", "RapidDrain", "RapidDrain", "ThinBreak",
                                       "Blowout", "Other", "Other", "Blowout", "UNK", "Water", "Urban",
                                       "Urban", "Rural", "Industrial", "Industrial", "Rural", "Mine",
                                       "Mine", "Wellsites", "EnSoftLin", "EnSoftLin", "EnSeismic", "EnSeismic",
                                       "HardLin", "HardLin", "TrSoftLin", "TrSoftLin", "TrSoftLin",
                                       "Crop", "RoughP", "RoughP", "TameP", "Industrial", "Water", "Water",
                                       "Water", "Water", "UNK"), Sector = c("Native", "Native", "Native",
                                                                            "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                            "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                            "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                            "Native", "Native", "Native", "Native", "RuralUrban", "RuralUrban",
                                                                            "RuralUrban", "RuralUrban", "Energy", "Misc", "Energy", "Misc",
                                                                            "Energy", "Energy", "Energy", "Energy", "Energy", "Transportation",
                                                                            "Transportation", "Transportation", "Transportation", "Transportation",
                                                                            "Agriculture", "Agriculture", "Agriculture", "Agriculture", "Agriculture",
                                                                            "Misc", "Misc", "Misc", "Misc", "Forestry")), row.names = c("Cy",
                                                                                                                                        "Len", "LenS", "LenSP", "LenT", "LenW", "Li", "Lo", "Ltc", "LtcC",
                                                                                                                                        "LtcD", "LtcH", "LtcS", "LtcR", "Sb", "Sy", "BdL", "CS", "Gr",
                                                                                                                                        "Sa", "SwG", "TB", "BlO", "LenA", "Ov", "SL", "UNK", "Water",
                                                                                                                                        "UrbanIndustrial", "UrbanResidence", "RuralResidentialIndustrial",
                                                                                                                                        "IndustrialSiteRural", "WindGenerationFacility", "OtherDisturbedVegetation",
                                                                                                                                        "MineSite", "PeatMine", "WellSite", "Pipeline", "TransmissionLine",
                                                                                                                                        "SeismicLineNarrow", "SeismicLineWide", "RoadHardSurface", "RailHardSurface",
                                                                                                                                        "RoadTrailVegetated", "RoadVegetatedVerge", "RailVegetatedVerge",
                                                                                                                                        "CultivationCrop", "CultivationAbandoned", "CultivationRoughPasture",
                                                                                                                                        "CultivationTamePasture", "HighDensityLivestockOperation", "BorrowpitsDugoutsSumps",
                                                                                                                                        "MunicipalWaterSewage", "Reservoirs", "Canals", "CutBlocks"), class = "data.frame"),
        north=structure(list(Label = c("DeciduousR", "Deciduous1", "Deciduous2",
                                       "Deciduous3", "Deciduous4", "Deciduous5", "Deciduous6", "Deciduous7",
                                       "Deciduous8", "MixedwoodR", "Mixedwood1", "Mixedwood2", "Mixedwood3",
                                       "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7", "Mixedwood8",
                                       "PineR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6",
                                       "Pine7", "Pine8", "WhiteSpruceR", "WhiteSpruce1", "WhiteSpruce2",
                                       "WhiteSpruce3", "WhiteSpruce4", "WhiteSpruce5", "WhiteSpruce6",
                                       "WhiteSpruce7", "WhiteSpruce8", "TreedBogR", "TreedBog1", "TreedBog2",
                                       "TreedBog3", "TreedBog4", "TreedBog5", "TreedBog6", "TreedBog7",
                                       "TreedBog8", "TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3",
                                       "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8",
                                       "TreedSwamp", "TreedSwamp", "TreedSwamp", "TreedSwamp", "TreedSwamp",
                                       "TreedSwamp", "TreedSwamp", "TreedSwamp", "TreedSwamp", "GrassHerb",
                                       "Shrub", "GraminoidFen", "Marsh", "ShrubbyBog", "ShrubbyFen",
                                       "ShrubbySwamp", "Water", "Urban", "Urban", "Rural", "Industrial",
                                       "Industrial", "Rural", "Mine", "Mine", "Wellsites", "EnSoftLin",
                                       "EnSoftLin", "EnSeismic", "EnSeismic", "HardLin", "HardLin",
                                       "TrSoftLin", "TrSoftLin", "TrSoftLin", "Crop", "RoughP", "RoughP",
                                       "TameP", "Industrial", "CCDeciduousR", "CCDeciduous1", "CCDeciduous2",
                                       "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR", "CCMixedwood1",
                                       "CCMixedwood2", "CCMixedwood3", "CCMixedwood4", "CCPineR", "CCPine1",
                                       "CCPine2", "CCPine3", "CCPine4", "CCWhiteSpruceR", "CCWhiteSpruce1",
                                       "CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4", "Bare",
                                       "Water", "Water", "Water", "Water", "SnowIce"), Sector = c("Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "Native", "Native", "Native", "Native", "Native", "Native", "Native",
                                                                                                  "RuralUrban", "RuralUrban", "RuralUrban", "RuralUrban", "Energy",
                                                                                                  "Misc", "Energy", "Misc", "Energy", "Energy", "Energy", "Energy",
                                                                                                  "Energy", "Transportation", "Transportation", "Transportation",
                                                                                                  "Transportation", "Transportation", "Agriculture", "Agriculture",
                                                                                                  "Agriculture", "Agriculture", "Agriculture", "Forestry", "Forestry",
                                                                                                  "Forestry", "Forestry", "Forestry", "Forestry", "Forestry", "Forestry",
                                                                                                  "Forestry", "Forestry", "Forestry", "Forestry", "Forestry", "Forestry",
                                                                                                  "Forestry", "Forestry", "Forestry", "Forestry", "Forestry", "Forestry",
                                                                                                  "Native", "Misc", "Misc", "Misc", "Misc", "Native")), row.names = c("DecidR",
                                                                                                                                                                      "Decid1", "Decid2", "Decid3", "Decid4", "Decid5", "Decid6", "Decid7",
                                                                                                                                                                      "Decid8", "MixedwoodR", "Mixedwood1", "Mixedwood2", "Mixedwood3",
                                                                                                                                                                      "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7", "Mixedwood8",
                                                                                                                                                                      "PineR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6",
                                                                                                                                                                      "Pine7", "Pine8", "SpruceR", "Spruce1", "Spruce2", "Spruce3",
                                                                                                                                                                      "Spruce4", "Spruce5", "Spruce6", "Spruce7", "Spruce8", "TreedBogR",
                                                                                                                                                                      "TreedBog1", "TreedBog2", "TreedBog3", "TreedBog4", "TreedBog5",
                                                                                                                                                                      "TreedBog6", "TreedBog7", "TreedBog8", "TreedFenR", "TreedFen1",
                                                                                                                                                                      "TreedFen2", "TreedFen3", "TreedFen4", "TreedFen5", "TreedFen6",
                                                                                                                                                                      "TreedFen7", "TreedFen8", "TreedSwampR", "TreedSwamp1", "TreedSwamp2",
                                                                                                                                                                      "TreedSwamp3", "TreedSwamp4", "TreedSwamp5", "TreedSwamp6", "TreedSwamp7",
                                                                                                                                                                      "TreedSwamp8", "GrassHerb", "Shrub", "GraminoidFen", "Marsh",
                                                                                                                                                                      "ShrubbyBog", "ShrubbyFen", "ShrubbySwamp", "Water", "UrbanIndustrial",
                                                                                                                                                                      "UrbanResidence", "RuralResidentialIndustrial", "IndustrialSiteRural",
                                                                                                                                                                      "WindGenerationFacility", "OtherDisturbedVegetation", "MineSite",
                                                                                                                                                                      "PeatMine", "WellSite", "Pipeline", "TransmissionLine", "SeismicLineNarrow",
                                                                                                                                                                      "SeismicLineWide", "RoadHardSurface", "RailHardSurface", "RoadTrailVegetated",
                                                                                                                                                                      "RoadVegetatedVerge", "RailVegetatedVerge", "CultivationCrop",
                                                                                                                                                                      "CultivationAbandoned", "CultivationRoughPasture", "CultivationTamePasture",
                                                                                                                                                                      "HighDensityLivestockOperation", "CCDecidR", "CCDecid1", "CCDecid2",
                                                                                                                                                                      "CCDecid3", "CCDecid4", "CCMixedwoodR", "CCMixedwood1", "CCMixedwood2",
                                                                                                                                                                      "CCMixedwood3", "CCMixedwood4", "CCPineR", "CCPine1", "CCPine2",
                                                                                                                                                                      "CCPine3", "CCPine4", "CCSpruceR", "CCSpruce1", "CCSpruce2",
                                                                                                                                                                      "CCSpruce3", "CCSpruce4", "Bare", "BorrowpitsDugoutsSumps", "Canals",
                                                                                                                                                                      "MunicipalWaterSewage", "Reservoirs", "SnowIce"), class = "data.frame"))
# Addressing issue with which sector OtherDisturbedVegetation is assigned to.
lt$south["OtherDisturbedVegetation", "Sector"] <- "RuralUrban"
lt$north["OtherDisturbedVegetation", "Sector"] <- "RuralUrban"

###############################
# Processes transition matrix #
###############################

# Define sectors (could be attribute or sector effects here)
chVeg$sector_use <- chVeg$sector
chSoil$sector_use <- chSoil$sector

#
# Southern matrix processing
#

compare_sets(chSoil$cr, rownames(lt$south))
compare_sets(chSoil$rf, rownames(lt$south))
chSoil$cr2 <- lt$south$Label[match(chSoil$cr, rownames(lt$south))]
chSoil$rf2 <- lt$south$Label[match(chSoil$rf, rownames(lt$south))]

# either soil is unknown or Sector is Forestry --> exclude these
s <- chSoil$cr2=="UNK" | chSoil$rf2=="UNK" | rownames(chSoil) == "UNK"
trSoil <- trSoil[,!s]
chSoil <- chSoil[!s,]

# reference=water is not part of the landbase, so does not count for averaging
s <- chSoil$rf2=="Water" | chSoil$cr2=="Water"

# but some of this is not Water->Water: we cannot attribute and it is implausible
sum(trSoil[,s])/sum(trSoil) - sum(trSoil[,chSoil$rf2=="Water" & chSoil$cr2=="Water"])/sum(trSoil)

# so we drop this ~3% together with open water
trSoil <- trSoil[,!s]
chSoil <- chSoil[!s,]

# Adjusting row sums after the removal of water and unknown transitions have occurred
# now we make sure rows sum to 1 (or 0)
rs <- rowSums(trSoil)
trSoil <- trSoil / ifelse(rs > 0, rs, 1) # Make sure row is at least divided by 1 to avoid 0 division issues

# take subset that contains only the south study region
trSoil <- trSoil[!kgrid$useN,]

rm(rs)
gc()

# 
# Northern matrix processing
#

compare_sets(chVeg$cr, rownames(lt$north))
setdiff(chVeg$cr, rownames(lt$north))
setdiff(rownames(lt$north), chVeg$cr)
compare_sets(chVeg$rf, rownames(lt$north))
chVeg$cr2 <- lt$north$Label[match(chVeg$cr, rownames(lt$north))]
chVeg$rf2 <- lt$north$Label[match(chVeg$rf, rownames(lt$north))]

# reference=water is not part of the landbase, so does not count for averaging
s <- chVeg$rf2=="Water" | chVeg$cr2=="Water"
chVeg[s,]

## but some of this is not Water->Water: we cannot attribute and it is implausible
sum(trVeg[,s])/sum(trVeg) - sum(trVeg[,chVeg$rf2=="Water" & chVeg$cr2=="Water"])/sum(trVeg)

## so we drop this ~3% together with open water
trVeg <- trVeg[,!s]
chVeg <- chVeg[!s,]

## now we make sure rows sum to 1 (or 0)
rn <- rowSums(trVeg)
trVeg <- trVeg / ifelse(rn > 0, rn, 1)

## take subset that contains only the north study region
trVeg <- trVeg[!kgrid$useS,]

# For each landcover type, create the necessary aggregated matrices

# Making matrices for each type
vegref <- list()
vegcur <- list()
soilref <- list()
soilcur <- list()

# Remove unessecary levels
chVeg <- droplevels(chVeg)
chSoil <- droplevels(chSoil)

# Veg
for (i in levels(chVeg$sector_use)) {
        
        j <- chVeg$sector_use == i
        
        vegref[[i]] <- groupSums(trVeg[,j], 2, chVeg$rf2[j])
        vegcur[[i]] <- groupSums(trVeg[,j], 2, chVeg$cr2[j])
        
}

# Soil
for (i in levels(chSoil$sector_use)) {
        
        j <- chSoil$sector_use == i
        
        soilref[[i]] <- groupSums(trSoil[,j], 2, chSoil$rf2[j])
        soilcur[[i]] <- groupSums(trSoil[,j], 2, chSoil$cr2[j])
        
        
}

rm(trSoil, trVeg)

# This is a template for storing predictions before summing
north.template <- matrix(data = NA, nrow = nrow(vegref[[1]]),
                         ncol = 14, 
                         dimnames = list(rownames(vegref[[1]]),
                                         c("Cur_Agriculture", "Ref_Agriculture",
                                           "Cur_Energy", "Ref_Energy",
                                           "Cur_Forestry", "Ref_Forestry",
                                           "Cur_Misc", "Ref_Misc",
                                           "Cur_RuralUrban", "Ref_RuralUrban",
                                           "Cur_Transportation", "Ref_Transportation",
                                           "Cur_Native", "Ref_Native")))

# This is a template for storing predictions before summing
south.template <- matrix(data = NA, nrow = nrow(soilref[[1]]),
                         ncol = 14, 
                         dimnames = list(rownames(soilref[[1]]),
                                         c("Cur_Agriculture", "Ref_Agriculture",
                                           "Cur_Energy", "Ref_Energy",
                                           "Cur_Forestry", "Ref_Forestry",
                                           "Cur_Misc", "Ref_Misc",
                                           "Cur_RuralUrban", "Ref_RuralUrban",
                                           "Cur_Transportation", "Ref_Transportation",
                                           "Cur_Native", "Ref_Native")))

# Template for the areas
area.north <- matrix(data = NA, nrow = nrow(vegref[[1]]),
                     ncol = 7, 
                     dimnames = list(rownames(vegref[[1]]),
                                     c("Area_Agriculture",
                                       "Area_Energy",
                                       "Area_Forestry", 
                                       "Area_Misc",
                                       "Area_RuralUrban",
                                       "Area_Transportation", 
                                       "Area_Native")))

area.south <- matrix(data = NA, nrow = nrow(soilref[[1]]),
                     ncol = 7, 
                     dimnames = list(rownames(soilref[[1]]),
                                     c("Area_Agriculture",
                                       "Area_Energy",
                                       "Area_Forestry", 
                                       "Area_Misc",
                                       "Area_RuralUrban",
                                       "Area_Transportation", 
                                       "Area_Native")))

# Calculate the areas of each sector by cell
for (i in levels(chVeg$sector)) {
        
        area.north[, paste0("Area_", i)] <- rowSums(vegcur[[i]])
        
}

for (i in levels(chSoil$sector)) {
        
        area.south[, paste0("Area_", i)] <- rowSums(soilcur[[i]])
        
}

# Subset south to the Grassland
area.south <- area.south[rownames(kgrid)[kgrid$NRNAME == "Grassland"], ]

# Merge the areas and save
area.sector <- rbind(area.north, area.south)
save(area.sector, file = paste0("data/processed/predictions/areas/sector-effect-areas_2021-01-12.Rdata"))

rm(area.south, area.sector, area.north)

#######################
# Species Predictions #
#######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load species coefficient tables and merge
load("results/coef/veg-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")
north.coef <- species.boot.results
load("results/coef/soil-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")
south.coef <- species.boot.results
rm(species.boot.results)

SPPn <- c("BCFR", "CATO", "WETO", "WOFR")
SPPs <- c("BCFR", "CATO", "WETO", "WOFR")
SPP <- sort(union(SPPn, SPPs))

for (spp in SPP) {
        
        print(paste(spp))
        flush.console()
        
        #########################
        # Coefficient alignment #
        #########################
        
        type <- "C" # combo species (N+S)
        M <- list(N=spp %in% SPPn, S=spp %in% SPPs)
        if (M$N & !M$S)
                type <- "N"
        if (!M$N & M$S)
                type <- "S"

        cfn <- if (type == "S")
                NULL else north.coef[[spp]]
        cfs <- if (type == "N")
                NULL else south.coef[[spp]]
        XclimS <- Xclim_nonb_S
        XclimN <- Xclim_nonb_N
        
        link.fun <- binomial()$linkfun
        ivnlink.fun <- binomial()$linkinv
        
        ##############
        # Prediction #
        ##############
        
        if (type != "N") {
                
                # Create the blank table for storing sector effects
                south.sector <- south.template
                
                gc()
                
                ## space-climate coefs
                bscl <- cfs$climate.coef[1, colnames(XclimS)]
                bscl[is.na(bscl)] <- 0 # this happens for habitat elements
                bspa <- cfs$paspen.coef[1, "paspen"]
                
                ## additive components for south
                muscl <- drop(XclimS %*% bscl)
                
                muspa <- kgrid[!kgrid$useN, "pAspen"] * bspa
                
                for (i in levels(chSoil$sector)) {
                        
                        ## south calculations for the i'th run
                        bscr <- link.fun(cfs$landscape.coef[1, colnames(soilcur[[i]])]) # current land cover
                        bsrf <- link.fun(cfs$landscape.coef[1, colnames(soilref[[i]])]) # reference land cover
                        
                        NScr <- matrix(muscl + muspa, nrow = nrow(soilcur[[i]]), ncol = ncol(soilcur[[i]]),
                                       dimnames = list(rownames(soilcur[[i]])))
                        NScr <- t(t(NScr) + bscr)
                        south.sector[, paste0("Cur_", i)] <- rowSums(soilcur[[i]] * ivnlink.fun(NScr))
                        
                        NSrf <- matrix(muscl + muspa, nrow = nrow(soilref[[i]]), ncol = ncol(soilref[[i]]),
                                       dimnames = list(rownames(soilcur[[i]])))
                        NSrf <- t(t(NSrf) + bsrf)
                        south.sector[, paste0("Ref_", i)]  <- rowSums(soilref[[i]]  * ivnlink.fun(NSrf))
                        
                        rm(NScr, NSrf, bscr, bsrf)
                        gc()
                        
                }
                
                # Adjust Forestry to 0
                south.sector[, c("Cur_Forestry", "Ref_Forestry")] <- 0
                
        } else {
                
                south.sector <- NULL
                
        }
        if (type != "S") {
                
                gc()
                
                # Create the blank table for storing sector effects
                north.sector <- north.template
                
                ## space-climate coefs
                bncl <- cfn$climate.coef[1, colnames(XclimN)]
                bncl[is.na(bncl)] <- 0 # this happens for habitat elements
                
                ## additive components for north
                muncl <- drop(XclimN %*% bncl)
                
                for (i in levels(chVeg$sector)) {
                        
                        ## north calculations for the i'th run
                        bncr <- link.fun(cfn$landscape.coef[1, colnames(vegcur[[i]])]) # current land cover
                        bnrf <- link.fun(cfn$landscape.coef[1, colnames(vegref[[i]])]) # reference land cover
                        
                        NNcr <- matrix(muncl, nrow = nrow(vegcur[[i]]), ncol = ncol(vegcur[[i]]))
                        NNcr <- t(t(NNcr) + bncr)
                        north.sector[, paste0("Cur_", i)] <- rowSums(vegcur[[i]] * ivnlink.fun(NNcr))
                        
                        NNrf <- matrix(muncl, nrow = nrow(vegref[[i]]), ncol = ncol(vegref[[i]]))
                        NNrf <- t(t(NNrf) + bnrf)
                        north.sector[, paste0("Ref_", i)] <- rowSums(vegref[[i]] * ivnlink.fun(NNrf))
                        
                        rm(NNcr, NNrf, bncr, bnrf)
                        gc()
                        
                }
                
        } else {
                
                north.sector <- NULL
                
        }
        
       
        save(north.sector, south.sector, 
             file=file.path(paste0("data/processed/predictions/sectoreffects/", spp, ".RData")))
        
        rm(south.sector, north.sector)
        
}

rm(list=ls())
gc()


