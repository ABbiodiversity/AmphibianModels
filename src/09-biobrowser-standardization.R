#
# Title: Biobrowser Standardization
# Created: February 28th, 2022
# Last Updated: February 28th, 2022
# Author: Brandon Allen
# Objective: Standardization of the amphibian coefficients and sector effects for the biobrowser
# Keywords: Initialization
# Notes: 
#

##################
# Initialization #
##################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
rm(list=ls())
gc()

# Load north coefficients
load("results/coef/veg-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")

# Define the array
north <- array(dim = c(4, 120, 100))
dimnames(north)[1] <- list(c("BCFR", "CATO", "WETO", "WOFR"))
dimnames(north)[2] <- list(c(colnames(species.boot.results[["BCFR"]]$landscape.coef), 
                             colnames(species.boot.results[["BCFR"]]$climate.coef)))

north.auc <- c(0,0,0,0)
names(north.auc) <- names(species.boot.results)

for(spp in names(species.boot.results)) {
            
            north[spp,,] <- rbind(qlogis(t(species.boot.results[[spp]]$landscape.coef)),
                                  t(species.boot.results[[spp]]$climate.coef))
            
            north.auc[spp] <- mean(species.boot.results[[spp]]$fit[, "auc_both"])
            
}

# Load south coefficients
load("results/coef/soil-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")

# Define the array
south <- array(dim = c(4, 41, 100))
dimnames(south)[1] <- list(c("BCFR", "CATO", "WETO", "WOFR"))
dimnames(south)[2] <- list(c(colnames(species.boot.results[["BCFR"]]$landscape.coef), 
                             "pAspen",
                             colnames(species.boot.results[["BCFR"]]$climate.coef)))

south.auc <- c(0,0,0,0)
names(south.auc) <- names(species.boot.results)

for(spp in names(species.boot.results)) {
            
            south[spp,,] <- rbind(qlogis(t(species.boot.results[[spp]]$landscape.coef)),
                                  t(species.boot.results[[spp]]$paspen.coef),
                                  t(species.boot.results[[spp]]$climate.coef))
            
            south.auc[spp] <- mean(species.boot.results[[spp]]$fit[, "auc_both"])
            
}

# Load base occurrence data
load("data/processed/amphibians-abundance-2015-2021_2022-01-05.Rdata")

# Load lookup table to site simplication
site.lookup <- read.csv("data/lookup/naming/landcover-lookup_2021-12-06.csv")

# Standardize to the site level
amphibian.occ$Site <- site.lookup$Site_Simple[match(amphibian.occ$Location, site.lookup$Location)]

# Store number of recordings
recordings <- data.frame(Occurrences = c(sum(ifelse(amphibian.occ$BCFR > 0, 1, 0)),
                                         sum(ifelse(amphibian.occ$CATO > 0, 1, 0)),
                                         sum(ifelse(amphibian.occ$WETO > 0, 1, 0)),
                                         sum(ifelse(amphibian.occ$WOFR > 0, 1, 0))))

data.update <- NULL

for (site in unique(amphibian.occ$Site)) {
            
            site.sample <- amphibian.occ[amphibian.occ$Site == site, ] # Subset sites
            species.max <- apply(site.sample[, c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP")], 2, max) # Identify max
            site.sample <- site.sample[1, ] # Select First row
            site.sample[, c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP")] <- species.max # Overwite maximum abundance
            data.update <- rbind(data.update, site.sample) # Rbind
            rm(site.sample)
}

# Calculate unique sites
# Store number of recordings
recordings$nSites <- c(sum(ifelse(data.update$BCFR > 0, 1, 0)),
                       sum(ifelse(data.update$CATO > 0, 1, 0)),
                       sum(ifelse(data.update$WETO > 0, 1, 0)),
                       sum(ifelse(data.update$WOFR > 0, 1, 0)))

# Create species lookup

species <- read.csv("data/lookup/amphibian_lookup_v2020.csv")
species <- species[species$code %in% c("BCFR", "CATO", "WETO", "WOFR"), ]
species <- data.frame(SpeciesID = species$code,
                      ScientificName = species$scientific_name,
                      TSNID = NA,
                      CommonName = species$common_name,
                      ModelNorth = TRUE,
                      ModelSouth = TRUE,
                      UseavailNorth = FALSE,
                      UseavailSouth = FALSE,
                      Occurrences = recordings$Occurrences,
                      nSites = recordings$nSites,
                      SizeNorth = NA,
                      SizeSouth = NA,
                      Nonnative = FALSE,
                      LinkHabitat = "logit",
                      LinkSpclim = "logit",
                      AUCNorth = north.auc,
                      AUCSouth = south.auc,
                      R2North = NA,
                      R2South = NA,
                      Comments = NA,
                      Group = "amphibians")


amphibians <- list(north, south, species)
names(amphibians) <- c("north", "south", "species")

save(amphibians, file = "results/submissions/biobrowser/amphibian-biobrowser-2022-02-28.Rdata")
