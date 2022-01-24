#
# Title: Intactness and sector effect summaries
# Created: January 22, 2022
# Last Updated: January 24th, 2022
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
source("src/intactness-calculation_functions_2019-10-28.R")

# Load the kgrid for Reporting Regions
load("data/base/kgrid/kgrid_table_km.Rdata")

kgrid <- kgrid[, c(1:8)]
kgrid$LinkID <- as.character(kgrid$Row_Col)

# Define the sector effect boundaries
# North
kgrid$North_Sector <- ifelse(kgrid$NR %in% c("Boreal", "Foothills", "Canadian Shield", "Rocky Mountain") |
                                     kgrid$NSR == "Dry Mixedwood", TRUE, FALSE)
kgrid$North_Sector[kgrid$NSR == "Dry Mixedwood" & kgrid$POINT_Y <= 56.7] <- FALSE
kgrid$North_Sector[kgrid$NSR == "Montane" & kgrid$POINT_X > -114] <- FALSE

# South
kgrid$South_Sector <- ifelse(kgrid$NR %in% c("Grassland", "Parkland") |
                                     kgrid$NSR == "Dry Mixedwood", TRUE, FALSE)
kgrid$South_Sector[kgrid$NSR == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- FALSE
kgrid$South_Sector[kgrid$NSR == "Montane" & kgrid$POINT_X > -114] <- TRUE

# # Rocky
# kgrid$Rocky_Sector <- ifelse(kgrid$NR %in% c("Rocky Mountain"), TRUE, FALSE)
# kgrid$Rocky_Sector[kgrid$NSR == "Montane" & kgrid$POINT_X > -114] <- FALSE

# Create thresholds
occ.threshold <- 0.01

# Create species list (With Mammals)
se.2018.path <- c(list.files("data/processed/predictions/comb/", full.names = TRUE))

species.names  <- c(list.files("data/processed/predictions/comb/", full.names = FALSE))

species.names <- gsub("-comb_2022-01-20.RData", "", species.names)

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

# Identify the list of reporting regions
reporting.regions <- data.frame(Folder = c("SectorEffects", "SectorEffects"),
                                Name = c("North_Sector", "South_Sector"),
                                Regions = c("North", "South"))

# Intactness summary

for (region.area in reporting.regions$Name) {
        
        # Subset the region of interest
        link.id <- kgrid[kgrid[, region.area], "LinkID"]
        folder.id <- reporting.regions[reporting.regions$Name == region.area, "Folder"]
        
        intactness.results <- NULL
        intactness.results <- si_se_summary_table(species.list.2018 = se.2018.path,
                                                  se.areas.2018 = area.2018,
                                                  region.id = link.id,
                                                  species.lookup = species.names, 
                                                  status.lookup = species.status,
                                                  threshold = occ.threshold, 
                                                  results.out = intactness.results)
        
        write.csv(intactness.results, file = paste0("results/tables/", folder.id, "/", region.area, "-intactness-se-results_", Sys.Date(), ".csv"),
                  row.names = FALSE)
        
        print(region.area)
        
}
