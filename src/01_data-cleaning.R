#
# Title: Cleaning of ABMI Amphibian Abundance data
# Created: October 28th, 2020
# Last Updated: January 24th, 2022
# Author: Brandon Allen
# Objective: Clean and filter the WildTrax data for Amphibian analyses.
# Keywords: Notes, Occurrence summaries, Calling Intensity, Landcover summaries, Climate summaries
#

######### 
# Notes # 
#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 1) Each recording represents a single observaton.
# 2) For definitions of the labels used to distinguish the types of ARU data, see the BU_acoustic_recording_analysis_protocol-v10
# 3) There are new labels for Wildtrax, we are going to use the 1m and 3m species recordings, but only the detections in the first minute.
# 4) New lookup for the vegetation and soil classes (v2020) are implemented
# 5) This analysis includes only data included in WildTrax. This means some information from other BU projects may not be included. 
# 6) Data was downloaded from all projects with amphibian detections on December 2, 2021
#
########################
# Occurrence summaries # 
########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##############################
# Occurrence standardization # 
##############################

# Clear memory
rm(list=ls())
gc()

# Load libraries
library(ggplot2)
library(mefa4)
library(opticut)

# Using the original data sets, create one merged document
column.lookup <- read.csv("data/lookup/column-id_lookup_v2021.csv")

occ.files <- list.files(path = "data/base/species/2021/", full.names = TRUE)
occ.files <- occ.files[!(occ.files %in% c("data/base/species/2021/english_column_definitions.csv",
                                          "data/base/species/2021/APPENDED_REPORT_2021.csv"))]

amphibian.data <- NA
for (x in 1:length(occ.files)) {
  
  # Load file
  temp.data <- read.csv(occ.files[x])
  
  # Standardize columns
  temp.col <- column.lookup$ColID[!(column.lookup$ColID %in% colnames(temp.data))]
  new.data <- data.frame(matrix(data = NA, nrow = nrow(temp.data), ncol = length(temp.col), dimnames = list(NULL, temp.col)))
  temp.data <- cbind(temp.data, new.data)
  temp.data <- temp.data[, column.lookup$ColID]
  
  amphibian.data <- rbind(amphibian.data, temp.data)
  print(x)
  
  rm(temp.col, new.data,temp.data)
  
}

write.csv(amphibian.data, file = "data/base/species/2021/APPENDED_REPORT_2021.csv", row.names = FALSE)
rm(x, amphibian.data, occ.files, column.lookup)

#####################
# Project Selection # 
#####################

# Load amphibian data 
amphibian.data <- read.csv("data/base/species/2021/APPENDED_REPORT_2021.csv")

# Keep projects based on the lookup table
project.lookup <- read.csv("data/lookup/project-lookup_v2021.csv")
project.lookup <- project.lookup$Project[project.lookup$Inclusion]
amphibian.data <- amphibian.data[amphibian.data$project %in% project.lookup, ]

#########################
# Survey Standarization # 
#########################

# Load amphibian lookup tables
tax <- read.csv("data/lookup/amphibian_lookup_v2020.csv")

# Standardize the recording intervals as characters
amphibian.data$min0_start <- as.character(amphibian.data$min0_start)
amphibian.data$min1_start <- as.character(amphibian.data$min1_start)
amphibian.data$min2_start <- as.character(amphibian.data$min2_start)

# Remove "Heavy" rain, wind, and noise
amphibian.data <- amphibian.data[!(amphibian.data[, "rain"] %in% "Heavy"), ]
amphibian.data <- amphibian.data[!(amphibian.data[, "wind"] %in% "Heavy"), ]
amphibian.data <- amphibian.data[!(amphibian.data[, "industry_noise"] %in% "Heavy"), ]
amphibian.data <- amphibian.data[!(amphibian.data[, "noise"] %in% "Heavy"), ]

# Keep only Confident and Confirmed Calls
amphibian.data <- amphibian.data[amphibian.data$confidence %in% c("Confident", "Confirmed"), ]

# Take the 1 and 3 minute recording intervals
amphibian.data <- amphibian.data[amphibian.data$method %in% c("1m 1SPM", "1m 2SPM", "3m 1SPM", "3m 2SPM"),]
amphibian.data$maxdur <- as.integer(sapply(strsplit(as.character(amphibian.data$method), "m"), "[[", 1))

# Split the characters of the start minutes to capture the first detection time
f <- function(v) {
  v <- strsplit(as.character(v), ",")
  v <- sapply(v, function(z) if (length(z) < 1) NA else z[1])
  v <- as.numeric(v)
  v
}
amphibian.data$min0_start_num <- f(amphibian.data$min0_start)
amphibian.data$min1_start_num <- f(amphibian.data$min1_start)
amphibian.data$min2_start_num <- f(amphibian.data$min2_start)

amphibian.data$Start <- strptime(paste(amphibian.data$recording_date, amphibian.data$recording_time),  "%Y-%m-%d %H:%M:%S")

# Compare the species codes from Wildtrax to the lookup table
compare_sets(tax$code, amphibian.data$species_code)
intersect(tax$code, amphibian.data$species_code)
setdiff(tax$code, amphibian.data$species_code)
setdiff(amphibian.data$species_code, tax$code)

# Standardize the common and scientific names for each recording
tx <- nonDuplicated(amphibian.data[,c("species_code", "scientific_name","species_english_name")], species_code, TRUE)
tx$m1 <- as.character(tax$code[match(tx$species_code, tax$code)])
tx$m1[is.na(tx$m1)] <- "Other"

# Match the species lookup table to create the presence/absence information for each station
amphibian.data$Spp <- as.factor(tx$m1)[match(amphibian.data$species_code, tx$species_code)]

###############################
# Survey time standardization # 
###############################

# Create information for the time of year, time of year cuts, time of day, and subset base on the recording time (11:30 - 2:30)
amphibian.data$ToY <- amphibian.data$Start$yday
amphibian.data$ToY2 <- amphibian.data$ToY * amphibian.data$ToY
amphibian.data$Year <- amphibian.data$Start$year + 1900
amphibian.data$site_stn <- paste0(amphibian.data$location, "_", amphibian.data$Year)
amphibian.data$ToYc <- as.integer(cut(amphibian.data$ToY, c(0, 105, 120, 140, 150, 160, 170, 180, 365)))
amphibian.data$replicate <- as.character(amphibian.data$Start)
amphibian.data$replicate <- gsub("[[:punct:]]", "", amphibian.data$replicate)
amphibian.data$replicate <- gsub("[[:space:]]", "", amphibian.data$replicate)
amphibian.data$visit <- paste0("ABMISM::", amphibian.data$site_stn, "::", amphibian.data$replicate)

amphibian.data$ToD <- amphibian.data$Start$hour + amphibian.data$Start$min / 60
amphibian.data$ToDx <- round(amphibian.data$ToD, 0)
amphibian.data$ToDc <- as.factor(ifelse(amphibian.data$ToDx %in% c(0, 1, 2, 24), "Midnight", "Morning"))

# Remove all recordings not completed at "Midnight"
amphibian.data <- amphibian.data[amphibian.data$ToDc == "Midnight", ]

#################################
# Vocalization standardization  # 
#################################

# Adjust abundance scores (CI 1 = 1, CI 2 = 2, CI 3 = 3, N/A = 0, ONE = 1, TMTC = 3, "TMTT" = 3)
abund.lookup <- data.frame(Orig = c("1", "2", "3", "4", "5", "7",
                                    "CI 1", "CI 2", "CI 3", "N/A", "ONE", "TMTC", "TMTT"),
                           Update = c(1, 2, 3, 3, 3, 3,
                                      1, 2, 3, 0, 1, 3, 3))
amphibian.data$abundance <- abund.lookup$Update[match(amphibian.data$abundance, abund.lookup$Orig)]

# If all species that are not in tax, abundance code of 0
amphibian.data$abundance[!(amphibian.data$species_code %in% tax$code)] <- 0

# Standardize the output
amphibian.occ <- amphibian.data[, c("project", "organization", "location", "Year", "site_stn", "latitude", "longitude", "hourly_weather_station_distance", "hourly_temp", "hourly_precipitation_mm", "ToY", "ToY2")]
amphibian.occ <- amphibian.occ[!duplicated(amphibian.occ), ]
colnames(amphibian.occ) <- c("Project", "Organization", "Location", "Year", "Site_Location", "Lat", "Long", "Weather_Station_Distance",
                                  "Hourly_Temp", "Hourly_Precipitation", "ToY", "ToY2")

# Add factor to account for stations that have been visited multiple times
visit.template <- data.frame(Site = names(table(amphibian.occ$Location)), 
                             visits = as.numeric(table(amphibian.occ$Location)))

amphibian.occ["visit"] <- visit.template$visits[match(amphibian.occ$Location, visit.template$Site)]
amphibian.occ["Site_Loc_ToY"] <- paste0(amphibian.occ$Site_Location, "_", amphibian.occ$ToY)
amphibian.occ <- amphibian.occ[!duplicated(amphibian.occ$Site_Loc_ToY), ]

# Append the abundance and the time to first detection information (SpeciesCode; SpeciesCode_Dect)

for (spp in tax$code) {
  
  spp.long <- NULL
  
  for (site in unique(amphibian.occ$Site_Location)) {

    site.subset <- amphibian.data[amphibian.data$site_stn == site, ]
    
    temp.subset <- NULL
    
    # If there are multiple recordings on a single day, take only the earliest recordings
    for(unique.days in unique(site.subset$recording_date)) {
      
      day.subset <- site.subset[site.subset$recording_date == unique.days, ]
      day.subset <- day.subset[day.subset$ToD == min(day.subset$ToD), ]
      temp.subset <- rbind(day.subset, temp.subset)
      
    }
    
    temp.subset$abundance[temp.subset$Spp != spp] <- 0
    temp.subset$min0_start_num[temp.subset$Spp != spp] <- 0
    temp.subset <- temp.subset[, c("site_stn", "ToY", "abundance", "min0_start_num")]
    temp.subset <- temp.subset[!duplicated(temp.subset), ]
    temp.subset <- temp.subset[order(temp.subset$abundance, decreasing = TRUE), ]
    temp.subset <- temp.subset[!duplicated(temp.subset[, 1:2]), ]
    temp.subset$Site_Loc_ToY <- paste0(temp.subset$site_stn, "_", temp.subset$ToY)
    
    spp.long <- rbind(spp.long, temp.subset)
    
    rm(temp.subset)
    
  }
  
  spp.long <- spp.long[, c("Site_Loc_ToY", "abundance", "min0_start_num")]
  colnames(spp.long)[2:3] <- c(spp, paste0(spp, "_", "min0_start"))
  amphibian.occ <- merge.data.frame(amphibian.occ, spp.long, by = "Site_Loc_ToY")
  
}

# Convert weather station distance to numeric
amphibian.occ$Weather_Station_Distance <- as.numeric(gsub("km", "", amphibian.occ$Weather_Station_Distance))

#
# Removal of incorrect species identifications
#

# WETO
amphibian.occ[amphibian.occ$Site_Loc_ToY %in% c("1578-NE_2017_123",
                                                "34-NE_2015_156",
                                                "1362-NE_2016_179",
                                                "8-SW_2018_157"), "WETO"] <- 0

# PLSP
amphibian.occ[amphibian.occ$Site_Loc_ToY %in% c("183-NW_2017_195"), "PLSP"] <- 0

# NLFR
amphibian.occ[amphibian.occ$Site_Loc_ToY %in% c("592-SW_2018_190",
                                                "W479-CL_2018_169"), "NLFR"] <- 0

# CSFR
amphibian.occ[amphibian.occ$Site_Loc_ToY %in% c("OG-EI-750-41-1_2017_148",
                                                "1311-NE_2018_164"), "CSFR"] <- 0

# CATO
amphibian.occ[amphibian.occ$Site_Loc_ToY %in% c("1644-NW_2016_146",
                                                "1645-NE_2016_144",
                                                "497-NW_2016_133",
                                                "555-NW_2017_169",
                                                "588-SW_2017_139",
                                                "653-NE_2017_186"), "CATO"] <- 0

# WOFR
# No corrections

# BCFR
# No corrections

# GPTO
# No corrections

# Save results
save(amphibian.occ, file = paste0("data/processed/amphibians-abundance-2015-2021_2022-01-05.Rdata"))

rm(list=ls())
gc()

#######################
# Landcover summaries # Information is spread out across multiple files. 
#######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Reload the cleaned occurrence data
load("data/processed/amphibians-abundance-2015-2021_2022-01-05.Rdata")

############### 
# 2015 - 2016 # 150m buffer, current and reference
###############

load("data/base/landcover/veg-hf_CameraARU_v6verified.Rdata")

# Load lookup table for landcover name changes
site.lookup <- read.csv("data/lookup/naming/landcover-lookup_2021-12-06.csv")

# Remove flagged sites from Eric (Remove from analysis or No information)
site.lookup <- site.lookup[!(site.lookup$Notes %in% c("No information", "Remove from analysis")), ]

# Filter to relevant subset
site.list <- site.lookup[site.lookup$Veg_File == "veg-hf_CameraARU_v6verified.Rdata", ]

# Subset the relevant sites
landcover.2015.2016 <- dd_150m 
landcover.2015.2016$veg_current <- landcover.2015.2016$veg_current[rownames(landcover.2015.2016$veg_current) %in% site.list$Veg_Name, ]
landcover.2015.2016$veg_reference <- landcover.2015.2016$veg_reference[rownames(landcover.2015.2016$veg_reference) %in% site.list$Veg_Name, ]
landcover.2015.2016$soil_current <- landcover.2015.2016$soil_current[rownames(landcover.2015.2016$soil_current) %in% site.list$Veg_Name, ]
landcover.2015.2016$soil_reference <- landcover.2015.2016$soil_reference[rownames(landcover.2015.2016$soil_reference) %in% site.list$Veg_Name, ]

# Get the 564m water
water.2015.2016 <- dd_564m
water.2015.2016$veg_current <- water.2015.2016$veg_current[rownames(water.2015.2016$veg_current) %in% site.list$Veg_Name, ]
water.2015.2016$soil_current <- water.2015.2016$soil_current[rownames(water.2015.2016$soil_current) %in% site.list$Veg_Name, ]

rm(site.list, dd_150m, dd_564m, dd_point)

############### 
# 2017 - 2019 # 150m buffer, current and reference
###############

load("data/base/landcover/veg-hf_ARU-2017-2019_Veg61-vHF.Rdata")
site.list <- site.lookup[site.lookup$Veg_File == "veg-hf_ARU-2017-2019_Veg61-vHF.Rdata", ]

# Subset the relevant sites
landcover.2017.2019 <- d_wide_150m # 150m
landcover.2017.2019$veg_current <- landcover.2017.2019$veg_current[rownames(landcover.2017.2019$veg_current) %in% site.list$Veg_Name, ]
landcover.2017.2019$veg_reference <- landcover.2017.2019$veg_reference[rownames(landcover.2017.2019$veg_reference) %in% site.list$Veg_Name, ]
landcover.2017.2019$soil_current <- landcover.2017.2019$soil_current[rownames(landcover.2017.2019$soil_current) %in% site.list$Veg_Name, ]
landcover.2017.2019$soil_reference <- landcover.2017.2019$soil_reference[rownames(landcover.2017.2019$soil_reference) %in% site.list$Veg_Name, ]

# Get the 564m water
water.2017.2019 <- d_wide_1km
water.2017.2019$veg_current <- water.2017.2019$veg_current[rownames(water.2017.2019$veg_current) %in% site.list$Veg_Name, ]
water.2017.2019$soil_current <- water.2017.2019$soil_current[rownames(water.2017.2019$soil_current) %in% site.list$Veg_Name, ]

rm(clim, d_long_pt, d_wide_150m, d_wide_1km, site.list)

############### Also includes new BU projects
# 2013 - 2021 # 150m buffer, current and reference
###############

load("data/base/landcover/amphibian-landcover-summaries_2021.RData")
site.list <- site.lookup[site.lookup$Veg_File == "amphibian-landcover-summaries_2021.RData", ]

# Subset the relevant sites
landcover.2013.2021 <- d_wide_150 # 150m
landcover.2013.2021$veg_current <- landcover.2013.2021$veg_current[rownames(landcover.2013.2021$veg_current) %in% site.list$Veg_Name, ]
landcover.2013.2021$veg_reference <- landcover.2013.2021$veg_reference[rownames(landcover.2013.2021$veg_reference) %in% site.list$Veg_Name, ]
landcover.2013.2021$soil_current <- landcover.2013.2021$soil_current[rownames(landcover.2013.2021$soil_current) %in% site.list$Veg_Name, ]
landcover.2013.2021$soil_reference <- landcover.2013.2021$soil_reference[rownames(landcover.2013.2021$soil_reference) %in% site.list$Veg_Name, ]

# Get the 564m water
water.2013.2021 <- d_wide_564
water.2013.2021$veg_current <- water.2013.2021$veg_current[rownames(water.2013.2021$veg_current) %in% site.list$Veg_Name, ]
water.2013.2021$soil_current <- water.2013.2021$soil_current[rownames(water.2013.2021$soil_current) %in% site.list$Veg_Name, ]

rm(d_long_150, d_long_564, d_wide_150, d_wide_564,  site.list)

###########
# Combine # 
###########

# Merge into single landcover file, standardize column names
veg.cur <- as.matrix(rbind(landcover.2013.2021$veg_current[, colnames(landcover.2013.2021$veg_current)], 
                           landcover.2015.2016$veg_current[, colnames(landcover.2013.2021$veg_current)], 
                           landcover.2017.2019$veg_current[, colnames(landcover.2013.2021$veg_current)]))

soil.cur <- as.matrix(rbind(landcover.2013.2021$soil_current[, colnames(landcover.2013.2021$soil_current)], 
                           landcover.2015.2016$soil_current[, colnames(landcover.2013.2021$soil_current)], 
                           landcover.2017.2019$soil_current[, colnames(landcover.2013.2021$soil_current)]))

veg.water <- as.matrix(rbind(water.2013.2021$veg_current[, colnames(water.2013.2021$veg_current)], 
                             water.2015.2016$veg_current[, colnames(water.2013.2021$veg_current)], 
                             water.2017.2019$veg_current[, colnames(water.2013.2021$veg_current)]))

soil.water <- as.matrix(rbind(water.2013.2021$soil_current[, colnames(water.2013.2021$soil_current)], 
                              water.2015.2016$soil_current[, colnames(water.2013.2021$soil_current)], 
                              water.2017.2019$soil_current[, colnames(water.2013.2021$soil_current)]))

###############################
# Clean to current categories #  
###############################

source("src/data-cleaning-functions_2020-10-28.R")
veg.lookup <- read.csv("data/lookup/lookup-veg-hf-age-v2020.csv")

veg.cur <- veg.cur[complete.cases(veg.cur), ] # Remove sites where vegetation information is unavailable
veg.cur <- landscape_hf_summary(data.in = veg.cur, landscape.lookup = veg.lookup, class.in = "ID", class.out = "UseInAnalysis_Simplified")
veg.cur <- as.data.frame(veg.cur / rowSums(veg.cur))

# Using the vegetation prediction matrix, create the necessary combined categories
pred.matrix <- read.csv("data/base/prediction-matrix/veg-prediction-matrix-CC_2021.csv")

for (col.id in colnames(pred.matrix)[-1]) {
        
        # Identify coefficient set
        coef.set <- pred.matrix$VegType[as.logical(pred.matrix[, col.id])]
        
        if (length(coef.set) == 1) {
                
                if (coef.set %in% c("WhiteSpruce", "Pine", "Deciduous", "Mixedwood", "BlackSpruce")) {
                        
                        coef.set <- paste0(rep(coef.set, 8), c("R", 1:8))
                        
                }
        }
        
        if (length(coef.set) == 1) {
                
                veg.cur[, col.id] <- veg.cur[, colnames(veg.cur) %in% coef.set]
                
        } else {
                
                veg.cur[, col.id] <- rowSums(veg.cur[, colnames(veg.cur) %in% coef.set])
                
        }

}

########
# Soil #
########

soil.lookup <- read.csv("data/lookup/lookup-soil-hf-v2020.csv")

# Using the vegetation prediction matrix, create the necessary combined categories
pred.matrix <- read.csv("data/base/prediction-matrix/soil-prediction-matrix_2021.csv")

# Remove the time of year and hourly temp information

soil.cur <- soil.cur[complete.cases(soil.cur),] # Remove sites without soil
soil.cur <- landscape_hf_summary(data.in = soil.cur, landscape.lookup = soil.lookup, class.in = "ID", class.out = "UseInAnalysis_Simplified")
soil.cur <- as.data.frame(soil.cur / rowSums(soil.cur))

for (col.id in colnames(pred.matrix)[-1]) {
        
        # Identify coefficient set
        coef.set <- pred.matrix$VegType[as.logical(pred.matrix[, col.id])]
        
        if (length(coef.set) == 1) {
                
                if (coef.set %in% c("WhiteSpruce", "Pine", "Deciduous", "Mixedwood", "BlackSpruce")) {
                        
                        coef.set <- paste0(rep(coef.set, 8), c("R", 1:8))
                        
                }
        }
        
        if (length(coef.set) == 1) {
                
                soil.cur[, col.id] <- soil.cur[, colnames(soil.cur) %in% coef.set]
                
        } else {
                
                soil.cur[, col.id] <- rowSums(soil.cur[, colnames(soil.cur) %in% coef.set])
                
        }
        
}

# Convert 564m water to proportions
veg.water <- veg.water / rowSums(veg.water)
soil.water <- soil.water / rowSums(soil.water)

veg.water <- data.frame(Site_Location = rownames(veg.water),
                        Waterkm = as.numeric(veg.water[, "Water"]))

soil.water <- data.frame(Site_Location = rownames(soil.water),
                        Waterkm = as.numeric(soil.water[, "Water"]))

##################### 
# Climate summaries # Information is spread out across multiple files. 
#####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# There are stations where no recordings occurred, but space climate information is available. These sites can be ignroed. 
# There are also sites where we have recordings, but there is no space climate information. These sites need to be identified.

############### 
# 2015 - 2016 # AHM, FFP, MAP, MAT, MCMT, MWMT, PET, pAspen
###############

load("data/base/landcover/ab-birds-dd-2018-11-29.RData")
dd <- dd[dd$PCODE == "ABMISM", ]
dd <- dd[dd$YEAR %in% c(2015, 2016), ]
dd["Site_Location"] <- paste0(gsub("_", "-", dd$SS), "_", dd$YEAR)

# Filter to relevant subset
site.list <- site.lookup[site.lookup$Veg_File == "veg-hf_CameraARU_v6verified.Rdata", ]

climate.2015.2016 <- dd[dd$Site_Location %in% site.list$Occ_Name, c("Site_Location", "AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "PET", "pAspen")]
climate.2015.2016 <- climate.2015.2016[!duplicated(climate.2015.2016), ]

# There are sites with missing climate data. Use the information from the broad ABMI site list
missing.sites <- site.list[!(site.list$Occ_Name %in% climate.2015.2016$Site_Location), "Occ_Name"]
missing.sites <- c(climate.2015.2016[is.na(climate.2015.2016$AHM), "Site_Location"], missing.sites)
climate.2015.2016 <- climate.2015.2016[!is.na(climate.2015.2016$AHM), ]

# Create blank data frame
missing.sites <- data.frame(Site_Location = missing.sites,
                            AHM = NA,
                            FFP = NA,
                            MAP = NA,
                            MAT = NA, 
                            MCMT = NA,
                            MWMT = NA,
                            PET = NA,
                            pAspen = NA)

# Load substitute values and lookup
missing.lookup <- read.csv("data/lookup/naming/climate-lookup_2020-11-18.csv")
climate.info <- read.csv("data/lookup/site-climate-summary_v2020.csv")

# Append values
missing.lookup <- cbind(missing.lookup, climate.info[match(missing.lookup$ABMI_Site, climate.info$SITE_ID), c("AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "PET", "pAspen_mean")])

missing.sites[, 2:9] <- missing.lookup[match(missing.lookup$Site_Location, missing.lookup$Site_Location), c("AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "PET", "pAspen_mean")]

# Merge
climate.2015.2016 <- rbind(climate.2015.2016, missing.sites)

# Correct naming
climate.2015.2016$Site_Location <- site.lookup$Veg_Name[match(climate.2015.2016$Site_Location, site.lookup$Occ_Name)]

rm(missing.lookup, missing.sites, climate.info, dd)

############### 
# 2017 - 2019 # AHM, FFP, MAP, MAT, MCMT, MWMT, PET, pAspen
###############

load("data/base/landcover/veg-hf_ARU-2017-2019_Veg61-vHF.Rdata")
climate.2017.2019 <- clim[rownames(clim) %in% rownames(landcover.2017.2019$veg_current), ]
climate.2017.2019["Site_Location"] <- rownames(climate.2017.2019)
climate.2017.2019 <- climate.2017.2019[, c("Site_Location", "AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "PET", "pAspen")]

# There is one quadrant with a missing paspen value. Replace with the mean for that site.
climate.2017.2019["1546-SW_2017", "pAspen"] <- 0.8869087

############### Also includes new BU projects
# 2013 - 2021 # AHM, FFP, MAP, MAT, MCMT, MWMT, PET, pAspen
############### NOTE THERE IS ONE EXTRA ROW FOR SOME REASON

load("data/base/landcover/amphibian-climate-summaries_2021.Rdata")
climate.2013.2021 <- climate[climate$UID %in% rownames(landcover.2013.2021$veg_current), ]
climate.2013.2021["Site_Location"] <- climate.2013.2021$UID
climate.2013.2021 <- climate.2013.2021[, c("Site_Location", "AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "Eref", "paspen")]
colnames(climate.2013.2021)[colnames(climate.2013.2021) %in% c("Eref", "paspen")] <- c("PET", "pAspen")

# Remove the one duplicated station
climate.2013.2021 <- climate.2013.2021[!duplicated(climate.2013.2021$Site_Location),]

###########
# Combine #
###########

climate.summaries <- rbind(climate.2015.2016, climate.2017.2019, climate.2013.2021)
climate.summaries$Site_Location <- site.lookup$Occ_Name[match(climate.summaries$Site_Location, site.lookup$Veg_Name)]
rownames(climate.summaries) <- climate.summaries$Site_Location

# Rename pAspen to paspen
colnames(climate.summaries)[colnames(climate.summaries) == "pAspen"] <- "paspen"

rm(clim, climate.2015.2016, climate.2017.2019, climate.2013.2021, d_long_pt, d_wide_150m, d_wide_1km, climate)

#####################################
# Match with Climate and Occurrence # We currently lose the 2021 occurrence information and the SAAM sites (waiting for Erin)
#####################################

rownames(veg.cur) <- site.lookup$Occ_Name[match(rownames(veg.cur), site.lookup$Veg_Name)]
rownames(soil.cur) <- site.lookup$Occ_Name[match(rownames(soil.cur), site.lookup$Veg_Name)]
veg.water$Site_Location <- site.lookup$Occ_Name[match(veg.water$Site_Location, site.lookup$Veg_Name)]
soil.water$Site_Location <- site.lookup$Occ_Name[match(soil.water$Site_Location, site.lookup$Veg_Name)]

veg.cur["Site_Location"] <- rownames(veg.cur)
site.lookup$Site_Location <- site.lookup$Occ_Name
climate.summaries <- merge.data.frame(climate.summaries, site.lookup[, c("Site_Location", "Natural_Region", "Natural_Subregion")], "Site_Location")
veg.cur <- merge.data.frame(climate.summaries, veg.cur, "Site_Location")
veg.cur <- merge.data.frame(amphibian.occ, veg.cur, "Site_Location")

# Add the 564m water
veg.cur <- merge.data.frame(veg.cur, veg.water, "Site_Location")

# veg.ref["Site_Location"] <- rownames(veg.ref)
# veg.ref <- merge.data.frame(comb.lookup, veg.ref, "Site_Location")
# veg.ref <- merge.data.frame(amphibian.occ, veg.ref, "Site_Location")
# veg.ref <- merge.data.frame(climate.summaries, veg.ref, "Site_Location")

soil.cur["Site_Location"] <- rownames(soil.cur)
soil.cur <- merge.data.frame(climate.summaries, soil.cur, "Site_Location")
soil.cur <- merge.data.frame(amphibian.occ, soil.cur, "Site_Location")

# Add the 564m water
soil.cur <- merge.data.frame(soil.cur, soil.water, "Site_Location")

# soil.ref["Site_Location"] <- rownames(soil.ref)
# soil.ref <- merge.data.frame(comb.lookup, soil.ref, "Site_Location")
# soil.ref <- merge.data.frame(amphibian.occ, soil.ref, "Site_Location")
# soil.ref <- merge.data.frame(climate.summaries, soil.ref, "Site_Location")

# Add additional climate variables
veg.cur["Lat2"] <- veg.cur[, "Lat"] * veg.cur[, "Lat"]
veg.cur["Long2"] <- veg.cur[, "Long"] * veg.cur[, "Long"]
veg.cur["LatLong"] <- veg.cur[, "Lat"] * veg.cur[, "Long"]
veg.cur["MWMT2"] <- veg.cur[, "MWMT"] * veg.cur[, "MWMT"]
veg.cur["MAT2"] <- veg.cur[, "MAT"] * veg.cur[, "MAT"]
veg.cur["MAPPET"] <- veg.cur[, "MAP"] * veg.cur[, "PET"]
veg.cur["MAPFFP"] <- veg.cur[, "MAP"] * veg.cur[, "FFP"]
veg.cur["MATAHM"] <- veg.cur[, "MAT"] * veg.cur[, "AHM"]
veg.cur["Waterkm2"] <- veg.cur[, "Waterkm"] * veg.cur[, "Waterkm"]

# Subset vegetation and save
# Modify the latitude of each site so those below 52.8 Lat are treated as further north
veg.cur["TrueLat"] <- veg.cur$Lat
veg.cur$Lat[veg.cur$Lat <= 52.8] <- 52.8
veg.cur <- veg.cur[!(veg.cur$Natural_Region %in% c("Grassland")), ]
save(veg.cur, file = paste0("data/processed/site/amphibian-abundance-data-north-150m_", Sys.Date(), ".RData"))

# Add additional climate variables
soil.cur["Lat2"] <- soil.cur[, "Lat"] * soil.cur[, "Lat"]
soil.cur["Long2"] <- soil.cur[, "Long"] * soil.cur[, "Long"]
soil.cur["LatLong"] <- soil.cur[, "Lat"] * soil.cur[, "Long"]
soil.cur["MWMT2"] <- soil.cur[, "MWMT"] * soil.cur[, "MWMT"]
soil.cur["MAT2"] <- soil.cur[, "MAT"] * soil.cur[, "MAT"]
soil.cur["MAPPET"] <- soil.cur[, "MAP"] * soil.cur[, "PET"]
soil.cur["MAPFFP"] <- soil.cur[, "MAP"] * soil.cur[, "FFP"]
soil.cur["MATAHM"] <- soil.cur[, "MAT"] * soil.cur[, "AHM"]
soil.cur["Waterkm2"] <- soil.cur[, "Waterkm"] * soil.cur[, "Waterkm"]

# Subset soil and save
# Add a TrueLat column so the habitat model scripts and the plotting scripts can use the same column naming structure
soil.cur["TrueLat"] <- soil.cur$Lat
soil.cur <- soil.cur[(soil.cur$Natural_Region %in% c("Grassland", "Parkland") | soil.cur$Natural_Subregion %in% "Dry Mixedwood") & soil.cur$Lat <= 56.4, ]
save(soil.cur, file = paste0("data/processed/site/amphibian-abundance-data-south-150m_", Sys.Date(), ".RData"))

rm(list=ls())
gc()

