#
# Title: Intactness and sector effect summaries
# Created: January 22, 2022
# Last Updated: February 1st, 2022
# Author: Brandon Allen
# Objectives: Calculate sector effect and intactness summaries for the two modeling regions
# Keywords: Summaries
# Notes: 

#############
# Summaries # 
#############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
rm(list=ls())
gc()

# Load functions
source("src/intactness-calculation_functions.R")

# Load the kgrid for Reporting Regions
load("data/base/kgrid/kgrid_table_km.Rdata")

kgrid <- kgrid[, c(1:8)]
kgrid$LinkID <- as.character(kgrid$Row_Col)

# Load the kgrid and define the reporting regions
kgrid$Regions <- ifelse(kgrid$NR %in% c("Boreal", "Foothills", "Canadian Shield", "Rocky Mountain") |
                                kgrid$NSR == "Dry Mixedwood", "Forested", "Prairie")
kgrid$Regions[kgrid$NSR == "Dry Mixedwood" & kgrid$Lat <= 56.7] <- "Prairie"
kgrid$Regions[kgrid$NSR == "Montane" & kgrid$Long > -114] <- "Prairie"

# Define the overlap region
load("data/base/kgrid/overlap-region.Rdata")

# Create species list
se.2018.path <- c(list.files("data/processed/predictions/sectoreffects/", full.names = TRUE))

species.names  <- c(list.files("data/processed/predictions/sectoreffects/", full.names = FALSE))

species.names <- gsub(".RData", "", species.names)

# Create species lookup table
species.status <- data.frame(SpeciesID = species.names,
                             Taxon = "Amphibians",
                             ScientificName = c("Pseudacris maculata", "Anaxyrus hemiophrys",
                                                "Anaxyrus boreas", "Lithobates sylvaticus"),
                             CommonName = c("Boreal Chorus Frog", "Canadian Toad",
                                            "Western Toad", "Wood Frog"),
                             Genus = "Species",
                             Model_north = TRUE,
                             Model_south = TRUE,
                             ExpertInclusion = TRUE)

# Regional summary
regional.results <- regional_reports(species.list = se.2018.path,
                                     species.lookup = species.status,
                                     summary.region = kgrid,
                                     overlap.region = overlap.region)

save(regional.results, file = "results/tables/sectoreffects/regional-sector-effects_HFI2018_2023-01-06.Rdata")


