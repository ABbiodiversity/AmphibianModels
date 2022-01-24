#
# Title: Functions used calculate intactness using kgrid predictions
# Created: January 19th, 2022
# Last Updated: January 24th, 2022
# Author: Brandon Allen
# Objectives: Load functions required for calculating intactness
# Keywords: Intactness + SE table
#

#########################
# Intactness + SE table # 
#########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

si_se_summary_table <- function(species.list.2018, se.areas.2018, region.id, species.lookup, status.lookup, threshold, results.out) {
        
        for( species in species.lookup ) {
                
                # Define species results table
                results.temp <- as.data.frame(matrix(nrow = 1, ncol = 3, dimnames = list(c(species), c("SpeciesID", "2018_SI",
                                                                                                       "2018_SI2"))))
                results.temp[species, "SpeciesID"] <- species
                
                ###################
                # Intactness 2018 #
                ###################
                
                # Load the provincial results and identify maximum abundance
                load.name <- species.list.2018[grep(as.character(species), species.list.2018)]
                load.name <- load.name[nchar(load.name) == min(nchar(load.name))]
                load(load.name)
                rm(load.name)
                
                # Merge predictions with sector areas
                sector.effect <- as.data.frame(sector.effect)
                sector.effect$LinkID <- rownames(sector.effect)
                sector.effect <- merge.data.frame(sector.effect, se.areas.2018, by = "LinkID")
                rownames(sector.effect) <- sector.effect$LinkID
                
                # Identify column names for the current, reference, and area groups
                cur.names <- colnames(sector.effect)[grep("Cur_", colnames(sector.effect))]
                ref.names <- colnames(sector.effect)[grep("Ref_", colnames(sector.effect))]
                area.names <- colnames(sector.effect)[grep("Area_", colnames(sector.effect))]
                area.names <- area.names[-grep("Native", area.names)]
                
                # Calculate Ref and Current sums
                sector.effect["Ref"] <- as.numeric(rowSums(sector.effect[, ref.names]))
                sector.effect["Cur_2018"] <- as.numeric(rowSums(sector.effect[, cur.names]))
                
                # Truncate to 99%
                q <- min(quantile(sector.effect$Ref, 0.999), quantile(sector.effect$Cur_2018, 0.999))
                total.ref <- sector.effect$Ref
                total.cur <- sector.effect$Cur_2018
                
                for (se.name in ref.names) {
                        
                        sector.effect[, se.name] <- ifelse(total.ref >= q , (sector.effect[, se.name] * (q/total.ref)), sector.effect[, se.name])
                        
                }
                
                for (se.name in cur.names) {
                        
                        sector.effect[, se.name] <- ifelse(total.cur >= q, (sector.effect[, se.name] * (q/total.cur)), sector.effect[, se.name])
                        
                }
                
                sector.effect$Cur_Native <- 0.5 * (sector.effect$Cur_Native + sector.effect$Ref_Native)
                sector.effect$Ref_Native <- sector.effect$Cur_Native 
                
                sector.effect["Ref"] <- as.numeric(rowSums(sector.effect[, ref.names]))
                sector.effect["Cur_2018"] <- as.numeric(rowSums(sector.effect[, cur.names]))
                
                # Calculate species max
                species.max <- as.numeric(quantile(sector.effect$Ref, 0.99)) * threshold
                
                # Check which LinkID in the region have predictions
                region.subset <- rownames(sector.effect)[rownames(sector.effect) %in% region.id]
                
                # Subset to the region of interest and define temporary data frame
                sector.effect <- sector.effect[region.subset, ]
                
                # Calculate intactness based on 2018 HFI
                results.temp[species, "2018_SI"] <- ifelse(sum(sector.effect[, "Cur_2018"]) >= sum(sector.effect[, "Ref"]), 
                                                           100 * (sum(sector.effect[, "Ref"]) / sum(sector.effect[, "Cur_2018"])), 
                                                           100 * (sum(sector.effect[, "Cur_2018"]) / sum(sector.effect[, "Ref"])))
                
                results.temp[species, "2018_SI2"] <- ifelse(100 * (sum(sector.effect[, "Cur_2018"]) / sum(sector.effect[, "Ref"])) >= 100, 
                                                            200 - results.temp[species, "2018_SI"], 
                                                            100 * (sum(sector.effect[, "Cur_2018"]) / sum(sector.effect[, "Ref"])))
                
                # Identify if species is an increaser or decreaser and if prediction should be kept
                
                results.temp[species, "Increaser"] <- ifelse(results.temp[species, "2018_SI2"] >= 100, "Increaser", "Decreaser")
                results.temp[species, "Pred_Keep"] <- ifelse(mean(sector.effect[, "Ref"]) >= species.max, TRUE, FALSE)
                
                gc()
                
                #################
                # Sector Effect #
                #################
                
                # Sector effect information is already loaded
                # Regional
                regional.names <- NULL # Store the regional names for organizing
                for (se.name in cur.names[-grep("Native", cur.names)]) {
                        
                        update.name <- gsub("Cur", "Total", se.name)
                        results.temp[update.name] <- ((sum(sector.effect[, se.name]) - sum(sector.effect[, gsub("Cur", "Ref", se.name)])) / sum(sector.effect[, ref.names])) * 100
                        regional.names <- c(regional.names, update.name)
                        rm(update.name)
                        
                }
                
                # Under HF
                under.hf.names <- NULL # Store the under HF names for organizing
                for (se.name in cur.names[-grep("Native", cur.names)]) {
                        
                        update.name <- gsub("Cur", "UnderHF", se.name)
                        results.temp[update.name] <- ((sum(sector.effect[, se.name]) - sum(sector.effect[, gsub("Cur", "Ref", se.name)])) / sum(sector.effect[, gsub("Cur", "Ref", se.name)])) * 100
                        under.hf.names <- c(under.hf.names, update.name)
                        rm(update.name)
                        
                }
                
                # Area
                # Calculate area based on the m2 results
                results.temp[area.names] <- colSums(sector.effect[, area.names]) / sum(sector.effect[, c(area.names, "Area_Native")]) * 100
                
                # Bind results
                results.out <- rbind(results.out, results.temp)
                rm(sector.effect, results.temp)
                
        }
        
        # Merge the results with the lookup table where info is available. 
        results.out <- merge.data.frame(status.lookup, results.out, by = "SpeciesID", all.y = TRUE)
        results.out <- results.out[, c("SpeciesID", "ScientificName", "CommonName", "Taxon",
                                       "Model_north", "Model_south", "Pred_Keep", "Increaser", "2018_SI", "2018_SI2",
                                       area.names, regional.names, under.hf.names)]
        
        return(results.out)
        
}
