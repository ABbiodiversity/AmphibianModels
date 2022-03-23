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

# Create thresholds
occ.threshold <- 0.01

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

# Load the areas with sector effect information
load("data/processed/predictions/areas/sector-effect-areas_2021-01-12.Rdata")
area.2018 <- as.data.frame(as.matrix(area.sector))
area.2018$LinkID <- rownames(area.2018)

rm(area.sector)
gc()

# Regional summary
regional.results <- regional_reports(species.list = se.2018.path,
                                     species.lookup = species.status,
                                     summary.region = kgrid)

save(regional.results, file = "results/tables/sectoreffects/regional-sector-effects_HFI2018_2022-02-01.Rdata")


