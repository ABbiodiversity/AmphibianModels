#
# Title: Standardize coefficients
# Created: January 29th, 2021
# Last Updated: December 15th, 2021
# Author: Brandon Allen
# Objectives: Standardize the coefficients so they are the same format as all other species
# Keywords: Northern Models
# Notes: 
#

###################
# Northern Models #
###################

# Clear memory
rm(list=ls())
gc()

load("results/coef/veg-coefficients-amphibians-detection-IVW_2022-01-20.Rdata")
species.list <- names(species.boot.results)

for (spp in species.list) {
        
        # Grass needs to be standardized to GrassHerb
        colnames(species.boot.results[[spp]]$landscape.coef) <- gsub("Grass", "GrassHerb", colnames(species.boot.results[[spp]]$landscape.coef))
        colnames(species.boot.results[[spp]]$landscape.se) <- gsub("Grass", "GrassHerb", colnames(species.boot.results[[spp]]$landscape.se))
        
        # UrbInd needs to be standardized as two separate coefficients (Urban; Industrial)
        colnames(species.boot.results[[spp]]$landscape.coef) <- gsub("UrbInd", "Industrial", colnames(species.boot.results[[spp]]$landscape.coef))
        colnames(species.boot.results[[spp]]$landscape.se) <- gsub("UrbInd", "Industrial", colnames(species.boot.results[[spp]]$landscape.se))
        
        species.boot.results[[spp]]$landscape.coef <- cbind(species.boot.results[[spp]]$landscape.coef, species.boot.results[[spp]]$landscape.coef[, "Industrial"])
        colnames(species.boot.results[[spp]]$landscape.coef)[ncol(species.boot.results[[spp]]$landscape.coef)] <- "Urban"
        
        species.boot.results[[spp]]$landscape.se <- cbind(species.boot.results[[spp]]$landscape.se, species.boot.results[[spp]]$landscape.se[, "Industrial"])
        colnames(species.boot.results[[spp]]$landscape.se)[ncol(species.boot.results[[spp]]$landscape.se)] <- "Urban"
        
        # Need to add Mine coefficient of 0
        species.boot.results[[spp]]$landscape.coef <- cbind(species.boot.results[[spp]]$landscape.coef, 0)
        colnames(species.boot.results[[spp]]$landscape.coef)[ncol(species.boot.results[[spp]]$landscape.coef)] <- "Mine"
        
        species.boot.results[[spp]]$landscape.se <- cbind(species.boot.results[[spp]]$landscape.se, 0)
        colnames(species.boot.results[[spp]]$landscape.se)[ncol(species.boot.results[[spp]]$landscape.se)] <- "Mine"
        
        # Need to add MineV coefficient of 0
        species.boot.results[[spp]]$landscape.coef <- cbind(species.boot.results[[spp]]$landscape.coef, 0)
        colnames(species.boot.results[[spp]]$landscape.coef)[ncol(species.boot.results[[spp]]$landscape.coef)] <- "MineV"
        
        species.boot.results[[spp]]$landscape.se <- cbind(species.boot.results[[spp]]$landscape.se, 0)
        colnames(species.boot.results[[spp]]$landscape.se)[ncol(species.boot.results[[spp]]$landscape.se)] <- "MineV"
        
        # Rename Blackspruce to TreedBog
        colnames(species.boot.results[[spp]]$landscape.coef) <- gsub("BlackSpruce", "TreedBog", colnames(species.boot.results[[spp]]$landscape.coef))
        colnames(species.boot.results[[spp]]$landscape.se) <- gsub("BlackSpruce", "TreedBog", colnames(species.boot.results[[spp]]$landscape.se))
        
        # Expand TreedFen for age classes
        colnames(species.boot.results[[spp]]$landscape.coef) <- gsub("TreedFen", "TreedFenR", colnames(species.boot.results[[spp]]$landscape.coef))
        colnames(species.boot.results[[spp]]$landscape.se) <- gsub("TreedFen", "TreedFenR", colnames(species.boot.results[[spp]]$landscape.se))
        
        for (x in c(1:8)) {
                
                species.boot.results[[spp]]$landscape.coef <- cbind(species.boot.results[[spp]]$landscape.coef, species.boot.results[[spp]]$landscape.coef[, "TreedFenR"])
                colnames(species.boot.results[[spp]]$landscape.coef)[ncol(species.boot.results[[spp]]$landscape.coef)] <- paste0("TreedFen", x)
                
                species.boot.results[[spp]]$landscape.se <- cbind(species.boot.results[[spp]]$landscape.se, species.boot.results[[spp]]$landscape.se[, "TreedFenR"])
                colnames(species.boot.results[[spp]]$landscape.se)[ncol(species.boot.results[[spp]]$landscape.se)] <- paste0("TreedFen", x)
                
        }

}

save(species.boot.results, file = "results/coef/veg-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")

###################
# Southern Models #
###################

# Clear memory
rm(list=ls())
gc()

load("results/coef/soil-coefficients-amphibians-detection-IVW_2022-01-20.Rdata")

species.list <- names(species.boot.results)

for (spp in species.list) {
        
        # UrbInd needs to be standardized as two separate coefficients (Urban; Industrial)
        colnames(species.boot.results[[spp]]$landscape.coef) <- gsub("UrbInd", "Industrial", colnames(species.boot.results[[spp]]$landscape.coef))
        colnames(species.boot.results[[spp]]$landscape.se) <- gsub("UrbInd", "Industrial", colnames(species.boot.results[[spp]]$landscape.se))
        
        species.boot.results[[spp]]$landscape.coef <- cbind(species.boot.results[[spp]]$landscape.coef, species.boot.results[[spp]]$landscape.coef[, "Industrial"])
        colnames(species.boot.results[[spp]]$landscape.coef)[ncol(species.boot.results[[spp]]$landscape.coef)] <- "Urban"
        
        species.boot.results[[spp]]$landscape.se <- cbind(species.boot.results[[spp]]$landscape.se, species.boot.results[[spp]]$landscape.se[, "Industrial"])
        colnames(species.boot.results[[spp]]$landscape.se)[ncol(species.boot.results[[spp]]$landscape.se)] <- "Urban"
        
        # Need to add Mine coefficient of 0
        species.boot.results[[spp]]$landscape.coef <- cbind(species.boot.results[[spp]]$landscape.coef, 0)
        colnames(species.boot.results[[spp]]$landscape.coef)[ncol(species.boot.results[[spp]]$landscape.coef)] <- "Mine"
        
        species.boot.results[[spp]]$landscape.se <- cbind(species.boot.results[[spp]]$landscape.se, 0)
        colnames(species.boot.results[[spp]]$landscape.se)[ncol(species.boot.results[[spp]]$landscape.se)] <- "Mine"
        
        # Need to add MineV coefficient of 0
        species.boot.results[[spp]]$landscape.coef <- cbind(species.boot.results[[spp]]$landscape.coef, 0)
        colnames(species.boot.results[[spp]]$landscape.coef)[ncol(species.boot.results[[spp]]$landscape.coef)] <- "MineV"
        
        species.boot.results[[spp]]$landscape.se <- cbind(species.boot.results[[spp]]$landscape.se, 0)
        colnames(species.boot.results[[spp]]$landscape.se)[ncol(species.boot.results[[spp]]$landscape.se)] <- "MineV"
        
}

save(species.boot.results, file = "results/coef/soil-coefficients-amphibians-detection-IVW-standardized_2022-01-20.Rdata")

rm(list=ls())
gc()




