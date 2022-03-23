#
# Title: ABMI Amphibian predictions and visualization
# Created: August 8th, 2018
# Last Updated: March 23rd, 2022
# Author: Brandon Allen
# Objective: Intactness analysis for the North model region
# Keywords: Occurrence, Seasonality, Coefficients, Climate Results, Relative abundance, Sector effects, Overlap Region
# Notes: 
#

##############
# Occurrence # 
##############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

library(abmi.themes)
library(ggplot2)
library(ggnewscale)
library(MetBrewer)
library(sf)
library(tidyverse)

load("data/processed/amphibians-abundance-2015-2021_2022-01-05.Rdata")

# Add fonts
add_abmi_fonts()

# Load lookup table to site simplication
site.lookup <- read.csv("data/lookup/naming/landcover-lookup_2021-12-06.csv")

# Standardize to the site level
amphibian.occ$Site <- site.lookup$Site_Simple[match(amphibian.occ$Location, site.lookup$Location)]

data.update <- NULL

for (site in unique(amphibian.occ$Site)) {
    
    site.sample <- amphibian.occ[amphibian.occ$Site == site, ] # Subset sites
    species.max <- apply(site.sample[, c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP")], 2, max) # Identify max
    site.sample <- site.sample[1, ] # Select First row
    site.sample[, c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP")] <- species.max # Overwite maximum abundance
    data.update <- rbind(data.update, site.sample) # Rbind
    rm(site.sample)
}

# Convert abundance to presence/absence
for (spp in c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP")) {
    
    data.update[, spp ] <- ifelse(data.update[, spp ] > 0, "Detected", "Undetected")
    
}

# Convert into shapefile of vorrect projection
data.update <- st_as_sf(x = data.update, 
                 coords = c("Long", "Lat"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Split into ABMI and BU groups
abmi.sites <- data.update[data.update$Organization == "ABMI", ]
bu.sites <- data.update[data.update$Organization == "BU", ]

# Load the shapefile for the natural subregions
boundary.in <- read_sf("data/base/gis/NRNAMEdissolve.shp")

# Combine all natural regions
boundary.in <- boundary.in %>% 
    summarise(Area = sum(Area))

# Make a plain background
boundary.in$Boundary <- "Boundary"

######################
# ABMI Sampling Grid #
######################

png(paste("results/figures/occurrence/ABMI-sampling-grid_", Sys.Date(), ".png", sep = ""),
    height = 1600,
    width = 1600, 
    res = 300)

ggplot() + 
    geom_sf(data = boundary.in, show.legend = FALSE) +
    geom_sf(data = abmi.sites, aes(colour = Project, fill = Project, shape = "21"), show.legend = TRUE) +
    scale_shape_manual(values = c(21), guide = "none") +
    scale_color_manual(values = met.brewer(name = "Renoir", n = 13, type = "continuous")) +
    scale_fill_manual(name = paste0("Project"), values = met.brewer(name = "Renoir", n = 13, type = "continuous"), guide = "none") +
    ggtitle(paste0("ABMI Sampling")) + 
    theme_light() +
    theme_abmi(font = "Montserrat") +
    theme(axis.title = element_text(size=36),
          axis.text.x = element_text(size=36),
          axis.text.y = element_text(size=36),
          title = element_text(size=36), 
          legend.text = element_text(size=20),
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(vjust = -6),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))

dev.off()

####################
# BU Sampling Grid #
####################

png(paste("results/figures/occurrence/BU-sampling-grid_", Sys.Date(), ".png", sep = ""),
    height = 1600,
    width = 1600, 
    res = 300)

ggplot() + 
    geom_sf(data = boundary.in, show.legend = FALSE) +
    geom_sf(data = bu.sites, aes(colour = Project, fill = Project, shape = "21"), show.legend = TRUE) +
    scale_shape_manual(values = c(21), guide = "none") +
    scale_color_manual(values = met.brewer(name = "Renoir", n = 7, type = "continuous")) +
    scale_fill_manual(name = paste0("Project"), values = met.brewer(name = "Renoir", n = 7, type = "continuous"), guide = "none") +
    ggtitle(paste0("BU Sampling")) + 
    theme_light() +
    theme_abmi(font = "Montserrat") +
    theme(axis.title = element_text(size=36),
          axis.text.x = element_text(size=36),
          axis.text.y = element_text(size=36),
          title = element_text(size=36), 
          legend.text = element_text(size=36),
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(vjust = -6),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))

dev.off()

##########################
# Combined Sampling Grid #
##########################

png(paste("results/figures/occurrence/combined-sampling-grid_", Sys.Date(), ".png", sep = ""),
    height = 2400,
    width = 1600, 
    res = 300)

ggplot() + 
    geom_sf(data = boundary.in, mapping = aes(colour = "#000000", fill = "#000000"), show.legend = FALSE) +
    geom_sf(data = data.update, aes(colour = "#E8A396", fill = "#E8A396", shape = "21"), show.legend = FALSE) +
    scale_shape_manual(values = c(21)) +
    scale_color_manual(values = c("#000000", "#000000")) +
    scale_fill_manual(values = alpha(c("#95A09A", "#e89682"), c(0.2, 1))) +
    ggtitle(paste0("All Amphibian Sampling")) + 
    theme_light() +
    theme_abmi(font = "Montserrat") +
    theme(axis.title = element_text(size=36),
          axis.text.x = element_text(size=36),
          axis.text.y = element_text(size=36),
          title = element_text(size=36),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))

dev.off()

###########################
# Single Species sampling #
###########################

# Load species naming information, subset to the four species of interest
species.lookup <- read.csv("data/lookup/amphibian_lookup_v2020.csv")
species.lookup <- species.lookup[species.lookup$code %in% c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP"), ]

for (spp in species.lookup$code) {
    
    png(paste("results/figures/occurrence/", spp, "-sampling_", Sys.Date(), ".png", sep = ""),
        height = 2400,
        width = 1800, 
        res = 300)
    
    print(ggplot() + 
        geom_sf(data = boundary.in, mapping = aes(colour = "#000000", fill = "#000000"), show.legend = FALSE) +
            scale_color_manual(values = c("#000000")) +
            scale_fill_manual(values = alpha(c("#95A09A"), c(0.2))) +
        new_scale_color() +
        new_scale_fill() +
        geom_sf(data = data.update, aes_string(colour = spp, fill = spp, shape = spp)) +
            scale_shape_manual(values = c(21,21)) +
            scale_color_manual(values = c("#000000", "#000000")) +
            scale_fill_manual(values = alpha(c("#e89682", "#2D415B"), c(1, 0.2))) +
            ggtitle(paste0(species.lookup[species.lookup$code == spp, "common_name"])) + 
            theme_light() +
            theme_abmi(font = "Montserrat") +
            theme(axis.title = element_text(size=36),
                  axis.text.x = element_text(size=36),
                  axis.text.y = element_text(size=36),
                  title = element_text(size=36),
                  legend.text = element_text(size=36),
                  legend.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  legend.position = c(0.25, 0.15)))
    
    dev.off()
    
}

###############
# Seasonality # 
###############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear memory
rm(list=ls())
gc()

# Load libraries and data
library(abmi.themes)
library(ggplot2)
library(sf)
library(tidyverse)

load("data/processed/amphibians-abundance-2015-2021_2022-01-05.Rdata")

# Add fonts
add_abmi_fonts()

# Load species naming information, subset to the four species of interest
species.lookup <- read.csv("data/lookup/amphibian_lookup_v2020.csv")
species.lookup <- species.lookup[species.lookup$code %in% c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP"), ]

for (spp in species.lookup$code) {
    
    amphibian.occ[, spp] <- ifelse(amphibian.occ[, spp] > 0, 1, 0)
    
}

amphibian.occ$Stations <- 1

# Aggregate
agg.data <- aggregate(x = amphibian.occ[, c("BCFR", "CATO", "WETO", "WOFR", "CSFR", "GPTO", "NLFR", "PLSP", "Stations")],  by = list(ToY = amphibian.occ$ToY), FUN = sum)

for (spp in species.lookup$code) {
    
    agg.data[, spp] <- agg.data[, spp] / agg.data$Stations
    
}

agg.data$ToY2 <- agg.data$ToY^2

# Visualize
for (spp in species.lookup$code) {
    
    # Subset species
    agg.data$temp <- agg.data[, spp]
    
    # Define maximum value
    max.value <- ifelse(max(agg.data$temp) + 0.1 > 1, 1, max(agg.data$temp) + 0.1)
    
    # Create figure
    png(paste("results/figures/seasonality/all/", spp, "-seasonality_", Sys.Date(), ".png", sep = ""),
        height = 1400,
        width = 2100, 
        res = 300)
    
    print(ggplot(data = agg.data, aes(x = ToY, y = temp)) + 
        geom_point(mapping = aes(color = "1"), show.legend = FALSE) +
        scale_color_manual(values = c("#2D415B")) +
        scale_fill_manual(values = c("#2D415B")) +
        xlab("Julian Calendar Day") +
        ylab("Frequency of Detections") +
        ylim(c(0,max.value)) +
        theme_light() +
            theme_abmi(font = "Montserrat") +
            theme(axis.title = element_text(size=36),
                  axis.text.x = element_text(size=36),
                  axis.text.y = element_text(size=36),
                  title = element_text(size=36),
                  legend.text = element_text(size=36),
                  legend.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1)))

    dev.off()
    
}

################
# Coefficients # 
################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # 21 colors
rm(list=ls())
gc()

library(abmi.themes)
library(ggplot2)
library(ggpubr)
library(raster)

source("src/visualization-functions.R")

# Add fonts
add_abmi_fonts()

##############
# Veg models #
##############

# Load coefficient and identify species list
load("results/coef/veg-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")
species.list <- names(species.boot.results)

# Bootstrap 90% confidence intervals
for (spp in 1:length(species.list)) {
    
    if (spp == 1) {
        
        coef.limits <- apply(species.boot.results[[spp]]$landscape.coef, 2, function(coef.id) quantile(x = coef.id, probs = c(0.10, 0.9)))
        
    } else {
        
        coef.limits <- rbind(coef.limits, apply(species.boot.results[[spp]]$landscape.coef, 2, function(coef.id) quantile(x = coef.id, probs = c(0.10, 0.9))))
        
    }
    
}

coef.limits <- as.data.frame(coef.limits)
rownames(coef.limits) <- c("BCFR Min", "BCFR Max",
                           "CATO Min", "CATO Max",
                           "WETO Min", "WETO Max",
                           "WOFR Min", "WOFR Max")
coef.limits <- t(coef.limits)

# Create the figure
for(spp in species.list) {
    
    coefficient_plots(species = spp, spp.coef = species.boot.results, spp.boot = coef.limits, region = "North")
    
}

rm(coef.limits, species.boot.results, species.list, spp)

###############
# Soil models #
###############

# Load coefficient and identify species list
load("results/coef/soil-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")
species.list <- names(species.boot.results)

# Bootstrap 90% confidence intervals
for (spp in 1:length(species.list)) {
    
    if (spp == 1) {
        
        coef.limits <- apply(species.boot.results[[spp]]$landscape.coef, 2, function(coef.id) quantile(x = coef.id, probs = c(0.10, 0.9)))
        
    } else {
        
        coef.limits <- rbind(coef.limits, apply(species.boot.results[[spp]]$landscape.coef, 2, function(coef.id) quantile(x = coef.id, probs = c(0.10, 0.9))))
        
    }
    
}

coef.limits <- as.data.frame(coef.limits)
rownames(coef.limits) <- c("BCFR Min", "BCFR Max",
                           "CATO Min", "CATO Max",
                           "WETO Min", "WETO Max",
                           "WOFR Min", "WOFR Max")
coef.limits <- t(coef.limits)

# Update the bootstrap limits to the appropriate names

coef.limits <- coef.limits[c("Loamy", "SandyLoam", "RapidDrain", "ClaySub", "ThinBreak", "Blowout", "Other",
                             "Crop", "TameP", "RoughP", "Wellsites", "Rural", "Urban", "Industrial"), ]

rownames(coef.limits) <- factor(c("Loamy", "Sandy/loamy", "Rapid Drain", "Clay", "Thin Break", "Blowout", "Other soil types",
         "Cropland", "Tame pasture", "Rough pasture", "Well sites", "Rural residential", "Urban/Industrial", "Industrial (rural)"), levels = c("Loamy", "Sandy/loamy", "Rapid Drain", "Clay", "Thin Break", "Blowout", "Other soil types",
                                                                                                                                               "Cropland", "Tame pasture", "Rough pasture", "Well sites", "Rural residential", "Urban/Industrial", "Industrial (rural)")) # Rename coefficients and update as factors

# Create model for all species, but don't publish results that are unreasonable.
# Create the figures

for(spp in species.list) {
    
    coefficient_plots(species = spp, spp.coef = species.boot.results, spp.boot = coef.limits, region = "South")
    
}

rm(coef.limits, species.boot.results, species.list, spp)

###################
# Climate Results # 
###################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

library(abmi.themes)
library(ggplot2)

source("src/visualization-functions_2020-10-28.R")

##############
# Veg models #
##############

# Load coefficient and identify species list
load("results/coef/veg-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")
species.list <- names(species.boot.results)

for(spp in species.list) {
    
    # Select climate results
    climate.data <- species.boot.results[[spp]]$climate.coef
    
    coef.limits <- apply(climate.data, 2, function(coef.id) quantile(x = coef.id, probs = c(0.10, 0.9), na.rm = TRUE))
    
    coef.limits <- as.data.frame(t(coef.limits))
    coef.limits$Variable <- rownames(coef.limits)
    
    # Pull out first value
    climate.data <- climate.data[1, ]
    climate.data <- climate.data[!is.na(climate.data)]
    climate.data <- data.frame(Variable = names(climate.data),
                               Value = as.numeric(climate.data))
    climate.data <- merge.data.frame(climate.data, coef.limits, by = "Variable")
    colnames(climate.data)[3:4] <- c("LowerCI", "UpperCI")
    
    # Remove intercept
    climate.data <- climate.data[climate.data$Variable != "Intercept", ]
    
    png(file = paste("results/figures/coefficients/north/", spp, "_climate_", Sys.Date(), ".png", sep = ""),
        width = 1200,
        height = 1200, 
        res = 300)
    
    print(ggplot(data = climate.data, aes(x = Variable, y = Value, fill = Variable, color = Variable), show.legend = FALSE) +
              geom_point(show.legend = FALSE) +
              geom_linerange(aes(ymin = LowerCI, ymax = UpperCI),
                             position = position_dodge(.9), show.legend = FALSE) + 
              scale_fill_abmi(palette = "main") +
              scale_color_abmi(palette = "main") + 
              theme_light() +
              theme_abmi(font = "Montserrat") +
              theme(axis.title = element_text(size=36),
                    axis.text.x = element_text(size=36, angle = 90, hjust = 1),
                    axis.text.y = element_text(size=36),
                    title = element_text(size=36),
                    legend.text = element_text(size=36),
                    legend.title = element_blank(),
                    axis.line = element_line(colour = "black"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1)))
    
    dev.off()
    
}

###############
# Soil models #
###############

# Load coefficient and identify species list
load("results/coef/soil-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")
species.list <- names(species.boot.results)

for(spp in species.list) {
    
    # Select climate results
    climate.data <- species.boot.results[[spp]]$climate.coef
    
    coef.limits <- apply(climate.data, 2, function(coef.id) quantile(x = coef.id, probs = c(0.10, 0.9), na.rm = TRUE))
    
    coef.limits <- as.data.frame(t(coef.limits))
    coef.limits$Variable <- rownames(coef.limits)
    
    # Pull out first value
    climate.data <- climate.data[1, ]
    climate.data <- climate.data[!is.na(climate.data)]
    climate.data <- data.frame(Variable = names(climate.data),
                               Value = as.numeric(climate.data))
    climate.data <- merge.data.frame(climate.data, coef.limits, by = "Variable")
    colnames(climate.data)[3:4] <- c("LowerCI", "UpperCI")
    
    # Remove intercept
    climate.data <- climate.data[climate.data$Variable != "Intercept", ]
    
    png(file = paste("results/figures/coefficients/south/", spp, "_climate_", Sys.Date(), ".png", sep = ""),
        width = 1200,
        height = 1200, 
        res = 300)
    
    print(ggplot(data = climate.data, aes(x = Variable, y = Value, fill = Variable, color = Variable), show.legend = FALSE) +
              geom_point(show.legend = FALSE) +
              geom_linerange(aes(ymin = LowerCI, ymax = UpperCI),
                             position = position_dodge(.9), show.legend = FALSE) + 
              scale_fill_abmi(palette = "main") +
              scale_color_abmi(palette = "main") + 
              theme_light() +
              theme_abmi(font = "Montserrat") +
              theme(axis.title = element_text(size=36),
                    axis.text.x = element_text(size=36, angle = 90, hjust = 1),
                    axis.text.y = element_text(size=36),
                    title = element_text(size=36),
                    legend.text = element_text(size=36),
                    legend.title = element_blank(),
                    axis.line = element_line(colour = "black"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1)))
    
    dev.off()
    
}

######################
# Relative Abundance # 
######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

library(abmi.themes)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(MetBrewer)
library(raster)
library(sf)
library(tidyverse)

# Add fonts
add_abmi_fonts()

# Load the shapefile for the provincial boundary
boundary.in <- read_sf("data/base/gis/NRNAMEdissolve.shp")

# Combine all natural regions
boundary.in <- boundary.in %>% 
    summarise(Area = sum(Area))

# Make a plain background
boundary.in$Boundary <- "Boundary"

# Load kgrid and overlap region
kgrid <- read_sf(dsn = "data/base/gis/ABMI.gdb", layer = "Grid_1KM")
load("data/base/kgrid/overlap-region.Rdata")
overlap.region <- overlap.region[kgrid$GRID_LABEL,] # Align grids

# For each species, load the appropriate predictions file and create figure
species.pred <- c("data/processed/predictions/sectoreffects/BCFR.RData", "data/processed/predictions/sectoreffects/CATO.RData",
                  "data/processed/predictions/sectoreffects/WETO.RData", "data/processed/predictions/sectoreffects/WOFR.RData")
species.names <- c("BCFR", "CATO", "WETO", "WOFR")

for (spp in species.names) {
    
    # Load data
    load(species.pred[grep(spp, species.pred)])
    
    # If north and south regions are available, merge. Otherwise, use what is available
    if(!is.null(north.sector)) {
        
        sector.effect <- north.sector
        
    }
    
    if(!is.null(south.sector)) {
        
        sector.effect <- south.sector
        
    }
    
    if(!(is.null(south.sector)) & !(is.null(south.sector))) {
        
        # averaging comes here
        north.sector <- north.sector[match(rownames(overlap.region), rownames(north.sector)), ]
        north.sector[is.na(north.sector)] <- 0
        rownames(north.sector) <- rownames(overlap.region)
        
        south.sector <- south.sector[match(rownames(overlap.region), rownames(south.sector)), ]
        south.sector[is.na(south.sector)] <- 0
        rownames(south.sector) <- rownames(overlap.region)
        
        sector.effect <- north.sector
        
        for (i in colnames(sector.effect)) {
            
            sector.effect[, i] <- overlap.region$wN * north.sector[, i] + (1-overlap.region$wN) * south.sector[, i]
            
        }
        
    }
    
    # Save the merged file for Joan
    sector.effect <- as.data.frame(sector.effect)
    sector.effect$LinkID <- rownames(sector.effect)
    sector.effect <- sector.effect[, c(ncol(sector.effect), 1:(ncol(sector.effect) - 1))]
    write.csv(sector.effect, file = paste0("results/submissions/2022-02-08/", spp, ".csv"),
              row.names = FALSE)
    
    # Calculate current and reference abundances
    sector.effect <- data.frame(LinkID = rownames(sector.effect),
                                Cur = rowSums(sector.effect[, grep("Cur_", colnames(sector.effect))]),
                                Ref = rowSums(sector.effect[, grep("Ref_", colnames(sector.effect))]))
    rownames(sector.effect) <- sector.effect$LinkID
    
    # Truncate based on 99th percentile
    q <- min(quantile(sector.effect$Ref, 0.999), quantile(sector.effect$Cur, 0.999))
    sector.effect$Ref <- ifelse(sector.effect$Ref > q, q, sector.effect$Ref)
    sector.effect$Cur <- ifelse(sector.effect$Cur > q, q, sector.effect$Cur)
    
    # Define the maximum abundance for the species so the scaling is appropriate
    max.abund <- round(max(sector.effect$Ref, sector.effect$Cur), 2)
    
    # Merge with LinkIDs
    kgrid$Cur <- sector.effect$Cur[match(kgrid$GRID_LABEL, rownames(sector.effect))]
    kgrid$Ref <- sector.effect$Ref[match(kgrid$GRID_LABEL, rownames(sector.effect))]
    
    kgrid$Difference <- kgrid$Cur - kgrid$Ref
    min.diff <- min(kgrid$Difference)
    max.diff <- max(kgrid$Difference)

    # Current Abundance
    png(file = paste0("results/figures/maps/", spp, "_current-abundance_", Sys.Date(), ".png"),
        width = 1800,
        height = 2400, 
        res = 300)
    
    cur.map <- ggplot() +
              geom_sf(data = boundary.in, show.legend = FALSE) +
              geom_sf(data = kgrid, aes(colour = Cur, fill = Cur, shape = "22"), show.legend = TRUE) +
              scale_shape_manual(values = c(22), guide = "none") +
              scale_fill_gradientn(name = paste0("Relative\nAbundance"), colors = rev(met.brewer(name = "Hiroshige", n = 100, type = "continuous")), limits = c(0,max.abund), guide = "colourbar") +
              scale_color_gradientn(colors = rev(met.brewer(name = "Hiroshige", n = 100, type = "continuous")), limits = c(0,max.abund), guide = "none") +
              theme_light() +
              theme_abmi(font = "Montserrat") +
              theme(axis.title = element_text(size=36),
                    axis.text.x = element_text(size=36),
                    axis.text.y = element_text(size=36),
                    title = element_text(size=36),
                    legend.text = element_text(size=36),
                    legend.title = element_text(size=36, lineheight = 0.25, vjust = -6),
                    axis.line = element_line(colour = "black"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    legend.position = c(0.20, 0.17))
    
    print(cur.map)
    
    dev.off()
    
    # Reference Abundance
    png(file = paste0("results/figures/maps/", spp, "_reference-abundance_", Sys.Date(), ".png"),
        width = 1800,
        height = 2400, 
        res = 300)
    
    ref.map <- ggplot() +
              geom_sf(data = boundary.in, show.legend = FALSE) +
              geom_sf(data = kgrid, aes(colour = Ref, fill = Ref, shape = "22"), show.legend = TRUE) +
              scale_shape_manual(values = c(22), guide = "none") +
              scale_fill_gradientn(name = paste0("Relative\nAbundance"), colors = rev(met.brewer(name = "Hiroshige", n = 100, type = "continuous")), limits = c(0,max.abund), guide = "colourbar") +
              scale_color_gradientn(colors = rev(met.brewer(name = "Hiroshige", n = 100, type = "continuous")), limits = c(0,max.abund), guide = "none") +
              theme_light() +
              theme_abmi(font = "Montserrat") +
              theme(axis.title = element_text(size=36),
                    axis.text.x = element_text(size=36),
                    axis.text.y = element_text(size=36),
                    title = element_text(size=36),
                    legend.text = element_text(size=36),
                    legend.title = element_text(size=36, lineheight = 0.25, vjust = -6),
                    axis.line = element_line(colour = "black"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    legend.position = c(0.20, 0.17))
    
    print(ref.map)
    
    dev.off()
    
    # Difference
    png(file = paste0("results/figures/maps/", spp, "_difference_", Sys.Date(), ".png"),
        width = 1800,
        height = 2400, 
        res = 300)
    
    diff.map <- ggplot() +
              geom_sf(data = boundary.in, show.legend = FALSE) +
              geom_sf(data = kgrid, aes(colour = Difference, fill = Difference, shape = "22"), show.legend = TRUE) +
              scale_shape_manual(values = c(22), guide = "none") +
              scale_fill_gradientn(name = paste0("Difference"), colors = (met.brewer(name = "Cassatt2", n = 100, type = "continuous")), limits = c(min.diff,max.diff), guide = "colourbar") +
              scale_color_gradientn(colors = (met.brewer(name = "Cassatt2", n = 100, type = "continuous")), limits = c(min.diff,max.diff), guide = "none") +
              theme_light() +
              theme_abmi(font = "Montserrat") +
              theme(axis.title = element_text(size=36),
                    axis.text.x = element_text(size=36),
                    axis.text.y = element_text(size=36),
                    title = element_text(size=36),
                    legend.text = element_text(size=36),
                    legend.title = element_text(size=36, lineheight = 0.25, vjust = -6),
                    axis.line = element_line(colour = "black"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    legend.position = c(0.20, 0.17))
    
    print(diff.map)
    
    dev.off()

}

##################
# Sector Effects #
##################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#
# Visualize boundaries
#

rm(list=ls())
gc()

library(abmi.themes)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(raster)
library(sf)
library(tidyverse)

# Add fonts
add_abmi_fonts()

# Load the shapefile for the provincial boundary
boundary.in <- read_sf("data/base/gis/NRNAMEdissolve.shp")

# Combine all natural regions
boundary.in <- boundary.in %>% 
    summarise(Area = sum(Area))

# Make a plain background
boundary.in$Boundary <- "Boundary"

# Load the kgrid for Reporting Regions
load("data/base/kgrid/kgrid_table_km.Rdata")

# Define three reporting regions
# North = "#A8AF8C"
# South = "#E8B3A6"
# Rocky = "#2D415B"

kgrid$Boundaries <- ifelse(kgrid$NRNAME %in% c("Grassland", "Parkland") |
                            (kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y < 57), "Prairie", "Forested")
kgrid$Boundaries[kgrid$NRNAME == "Rocky Mountain"] <- "Forested"
kgrid$Boundaries[kgrid$NSRNAME == c("Montane") & kgrid$POINT_X > -114] <- "Prairie"

kgrid <- kgrid[, c("Row_Col", "Boundaries")]

# Load kgrid for mapping
kgrid.map <- read_sf(dsn = "data/base/gis/ABMI.gdb", layer = "Grid_1KM")
kgrid <- kgrid[kgrid.map$GRID_LABEL, ] # Sort grid labels
kgrid.map$Boundaries <- kgrid$Boundaries

png(file = paste0("results/figures/sectoreffects/sector-effect-regions_", Sys.Date(), ".png"),
    width = 2400,
    height = 3600, 
    res = 300)

sector.regions <- ggplot() +
    geom_sf(data = boundary.in, show.legend = FALSE) +
    geom_sf(data = kgrid.map, aes(colour = Boundaries, fill = Boundaries, shape = "22"), show.legend = TRUE) +
    scale_color_manual(values = c("#A8AF8C", "#E8B3A6"), labels = c("Forested", "Prairie")) +
    scale_fill_manual(values = c("#A8AF8C", "#E8B3A6"), labels = c("Forested", "Prairie")) + 
    scale_shape_manual(values = c(22), guide = "none") +
    theme_light() +
    theme_abmi(font = "Montserrat") +
    theme(axis.title = element_text(size=36),
          axis.text.x = element_text(size=36),
          axis.text.y = element_text(size=36),
          title = element_text(size=36),
          legend.text = element_text(size=36),
          legend.title = element_text(size=36, lineheight = 0.25, vjust = -6),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = c(0.25, 0.15))

print(sector.regions)

dev.off()

#
# Visualize sector effects
#

load("results/tables/SectorEffects/regional-sector-effects_HFI2018_2022-02-01.Rdata")

for (name.id in names(regional.results)) {
    
    #####################
    # Species Formating #
    #####################
    
    species.subset <- regional.results[[name.id]]
    
    # Creation of the data frame required for sector effect plotting
    region.se <- data.frame(SpeciesID = rep(species.subset$SpeciesID, 5), 
                            Sector = c(rep("Agriculture", nrow(species.subset)), rep("Energy", nrow(species.subset)), rep("Forestry", nrow(species.subset)), rep("Urban/Industrial", nrow(species.subset)), rep("Transportation", nrow(species.subset))),
                            Regional = c(species.subset$Total_Agriculture, species.subset$Total_Energy, species.subset$Total_Forestry, species.subset$Total_RuralUrban, species.subset$Total_Transportation),
                            UnderHF = c(species.subset$UnderHF_Agriculture, species.subset$UnderHF_Energy, species.subset$UnderHF_Forestry, species.subset$UnderHF_RuralUrban, species.subset$UnderHF_Transportation))
    
    # Need to do a check for the NA values or sector effects over 100%
    region.se[is.na(region.se)] <- 0
    region.se$Sector  <- factor(region.se$Sector, levels = c("Agriculture", "Energy", "Forestry", "Urban/Industrial", "Transportation"))
    region.se$Regional[region.se$Regional > 100] <- 100
    region.se$UnderHF[region.se$UnderHF > 100] <- 100
    
    # Sector effect colors
    colour.pal <- c("#FFD25D", "#836F9A", "#77AA67", "#EB7170", "#4E4E4E")
    
    for (species in unique(region.se$SpeciesID)) {
        
        # Subset species
        species.subset <- region.se[region.se$SpeciesID == species, ]
        
        under.hf <- ggplot(data = species.subset, aes(x = Sector, y = UnderHF, fill = Sector)) +
            geom_bar(stat = "identity", fill = colour.pal, color = "#000000") +
            guides(scale = "none") + 
            labs(x = "Sectors", y = "Under Footprint Sector Effect (%)") +
            scale_y_continuous(limits = c(-100, 100),
                               labels = c(-100, -50, 0, 50, 100)) +
            theme_light() +
            theme_abmi(font = "Montserrat") +
            theme(axis.title = element_text(size=36),
                  axis.text.x = element_text(size=36, angle = 90, hjust = 1),
                  axis.text.y = element_text(size=36),
                  title = element_text(size=36),
                  legend.text = element_text(size=36),
                  legend.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
        
        region.hf <- ggplot(data = species.subset, aes(x = Sector, y = Regional, fill = Sector)) +
            geom_bar(stat = "identity", fill = colour.pal, color = "#000000") +
            guides(scale = "none") + 
            labs(x = "Sectors", y = "Regional Footprint Sector Effect (%)") +
            scale_y_continuous(limits = c(-100, 100),
                               labels = c(-100, -50, 0, 50, 100)) +
            theme_light() +
            theme_abmi(font = "Montserrat") +
            theme(axis.title = element_text(size=36),
                  axis.text.x = element_text(size=36, angle = 90, hjust = 1),
                  axis.text.y = element_text(size=36),
                  title = element_text(size=36),
                  legend.text = element_text(size=36),
                  legend.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
        
        png(paste0("results/figures/sectoreffects/", name.id, "/", species, ".png"),
            height = 1800,
            width = 3600,
            res = 300)
        
        print(ggarrange(region.hf, under.hf, ncol = 2))
        
        dev.off()
        
        
    }
    
    
}

###################
# Linear Features #
###################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#
# Visualize boundaries
#

rm(list=ls())
gc()

library(abmi.themes)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(raster)
library(sf)
library(tidyverse)

# Add fonts
add_abmi_fonts()

#
# Visualize sector effects
#

load("results/tables/sectoreffects/regional-linear-features_HFI2018_2022-02-01.Rdata")

############
# Forested #
############

# Filter for species predicted as TRUE
species.subset <- regional.results$North

# Creation of the data frame required for sector effect plotting
region.se <- data.frame(SpeciesID = rep(species.subset$SpeciesID, 4), 
                        Sector = c(rep("Energy Seismic", nrow(species.subset)), rep("Energy Soft", nrow(species.subset)), rep("Hard", nrow(species.subset)), rep("Transportation Soft", nrow(species.subset))),
                        UnderHF = c(species.subset$UnderHF_EnSeismic, species.subset$UnderHF_EnSoftLin, species.subset$UnderHF_HardLin, species.subset$UnderHF_TrSoftLin))

# Need to do a check for the NA values or sector effects over 100%
region.se[is.na(region.se)] <- 0
region.se$Sector  <- factor(region.se$Sector, levels = c("Energy Seismic", "Energy Soft", "Hard", "Transportation Soft"))
region.se$UnderHF[region.se$UnderHF > 100] <- 100

# Sector effect colors
colour.pal <- c("#802417", "#e8b960", "#646e3b", "#17486f")

for (species in unique(region.se$SpeciesID)) {
    
    # Subset species
    species.subset <- region.se[region.se$SpeciesID == species, ]
    
    # Only report on the under footprint
    under.hf <- ggplot(data = species.subset, aes(x = Sector, y = UnderHF, fill = Sector)) +
        geom_bar(stat = "identity", fill = colour.pal, color = "#000000") +
        guides(scale = "none") + 
        labs(x = "Linear Features", y = "Under Footprint Sector Effect (%)") +
        scale_y_continuous(limits = c(-100, 100),
                           labels = c(-100, -50, 0, 50, 100)) +
        theme_light() +
        theme_abmi(font = "Montserrat") +
        theme(axis.title = element_text(size=36),
              axis.text.x = element_text(size=36, angle = 90, hjust = 1),
              axis.text.y = element_text(size=36),
              title = element_text(size=36),
              legend.text = element_text(size=36),
              legend.title = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    png(paste0("results/figures/sectoreffects/north/", species, "-linear-features.png"),
        height = 1800,
        width = 1800,
        res = 300)
    
    print(under.hf)
    
    dev.off()
    
}

###########
# Prairie #
###########

# Filter for species predicted as TRUE
species.subset <- regional.results$South

# Creation of the data frame required for sector effect plotting
region.se <- data.frame(SpeciesID = rep(species.subset$SpeciesID, 4), 
                        Sector = c(rep("Energy Seismic", nrow(species.subset)), rep("Energy Soft", nrow(species.subset)), rep("Hard", nrow(species.subset)), rep("Transportation Soft", nrow(species.subset))),
                        UnderHF = c(species.subset$UnderHF_EnSeismic, species.subset$UnderHF_EnSoftLin, species.subset$UnderHF_HardLin, species.subset$UnderHF_TrSoftLin))

# Need to do a check for the NA values or sector effects over 100%
region.se[is.na(region.se)] <- 0
region.se$Sector  <- factor(region.se$Sector, levels = c("Energy Seismic", "Energy Soft", "Hard", "Transportation Soft"))
region.se$UnderHF[region.se$UnderHF > 100] <- 100

# Sector effect colors
colour.pal <- c("#802417", "#e8b960", "#646e3b", "#17486f")

for (species in unique(region.se$SpeciesID)) {
    
    # Subset species
    species.subset <- region.se[region.se$SpeciesID == species, ]
    
    # Only report on the under footprint
    under.hf <- ggplot(data = species.subset, aes(x = Sector, y = UnderHF, fill = Sector)) +
        geom_bar(stat = "identity", fill = colour.pal, color = "#000000") +
        guides(scale = "none") + 
        labs(x = "Linear Features", y = "Under Footprint Sector Effect (%)") +
        scale_y_continuous(limits = c(-100, 100),
                           labels = c(-100, -50, 0, 50, 100)) +
        theme_light() +
        theme_abmi(font = "Montserrat") +
        theme(axis.title = element_text(size=36),
              axis.text.x = element_text(size=36, angle = 90, hjust = 1),
              axis.text.y = element_text(size=36),
              title = element_text(size=36),
              legend.text = element_text(size=36),
              legend.title = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    png(paste0("results/figures/sectoreffects/south/", species, "-linear-features.png"),
        height = 1800,
        width = 1800,
        res = 300)
    
    print(under.hf)
    
    dev.off()
    
}

##################
# Overlap Region #
##################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#
# Visualize boundaries
#

rm(list=ls())
gc()

library(abmi.themes)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(raster)
library(sf)
library(tidyverse)

# Add fonts
add_abmi_fonts()

# Load the shapefile for the provincial boundary
boundary.in <- read_sf("data/base/gis/NRNAMEdissolve.shp")

# Combine all natural regions
boundary.in <- boundary.in %>% 
    summarise(Area = sum(Area))

# Make a plain background
boundary.in$Boundary <- "Boundary"

# Load the kgrid with the overlap region
load("data/base/kgrid/overlap-region.Rdata")
overlap.region <- overlap.region[, c("Row_Col", "wN")]

# Load kgrid for mapping
kgrid.map <- read_sf(dsn = "data/base/gis/ABMI.gdb", layer = "Grid_1KM")
overlap.region <- overlap.region[kgrid.map$GRID_LABEL, ] # Sort grid labels
kgrid.map$wN <- overlap.region$wN

png(file = paste0("results/figures/overlap/prediction-overlap_", Sys.Date(), ".png"),
    width = 2400,
    height = 3600, 
    res = 300)

overlap.regions <- ggplot() +
    geom_sf(data = boundary.in, show.legend = FALSE) +
    geom_sf(data = kgrid.map, aes(colour = wN, fill = wN, shape = "22"), show.legend = TRUE) +
    scale_shape_manual(values = c(22), guide = "none") +
    scale_fill_gradientn(name = paste0("North Weight"), colors = met.brewer(name = "Homer2", n = 200, type = "continuous")[101:200], limits = c(0,1), guide = "colourbar") +
    scale_color_gradientn(colors = met.brewer(name = "Homer2", n = 200, type = "continuous")[101:200], limits = c(0,1), guide = "none") +
    theme_light() +
    theme_abmi(font = "Montserrat") +
    theme(axis.title = element_text(size=36),
          axis.text.x = element_text(size=36),
          axis.text.y = element_text(size=36),
          title = element_text(size=36),
          legend.text = element_text(size=36),
          legend.title = element_text(size=36, lineheight = 0.25, vjust = -4),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = c(0.25, 0.15))

print(overlap.regions)

dev.off()
