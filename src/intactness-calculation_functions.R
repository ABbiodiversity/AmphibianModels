#
# Title: Functions used calculate intactness using kgrid predictions
# Created: January 19th, 2022
# Last Updated: February 1st, 2022
# Author: Brandon Allen
# Objectives: Load functions required for calculating intactness
# Keywords: Custom reporting (single), Regional reporting
#

####################
# Custom reporting # Single region
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bii_se_summary <- function(species.list, se.areas, region.id, species.lookup, status.lookup, overlap, threshold, results.out) {
        
        for( species in species.lookup ) {
                
                # Define species results table
                results.temp <- as.data.frame(matrix(nrow = 1, ncol = 5, dimnames = list(c(species), c("SpeciesID", "2010_SI",
                                                                                                       "2010_SI2", "2018_SI",
                                                                                                       "2018_SI2"))))
                
                ###################
                # Intactness 2010 #
                ###################
                
                # Load the provincial results
                load.name <- species.list.2010[grep(as.character(species), species.list.2010)]
                load.name <- load.name[nchar(load.name) == min(nchar(load.name))]
                load(load.name)
                rm(load.name)
                
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
                
                rm(north.sector, south.sector)
                gc()
                
                # Identify column names for the current, reference, and area groups
                cur.names <- colnames(sector.effect)[grep("Cur_", colnames(sector.effect))]
                ref.names <- colnames(sector.effect)[grep("Ref_", colnames(sector.effect))]
                
                # Calculate Ref and Current sums
                sector.effect <- as.data.frame(sector.effect)
                sector.effect["Ref"] <- as.numeric(rowSums(sector.effect[, ref.names]))
                sector.effect["Cur_2010"] <- as.numeric(rowSums(sector.effect[, cur.names]))
                
                # Truncate to 99%
                q <- min(quantile(sector.effect$Ref, 0.999), quantile(sector.effect$Cur_2010, 0.999))
                sector.effect$Ref <- ifelse(sector.effect$Ref > q, q, sector.effect$Ref)
                sector.effect$Cur_2010 <- ifelse(sector.effect$Cur_2010 > q, q, sector.effect$Cur_2010)
                
                # Check which LinkID in the region have predictions
                region.subset <- rownames(sector.effect)[rownames(sector.effect) %in% region.id]
                
                # If there are no predictions made for a species within a boundary of interest, default values to NA
                if(length(region.subset) == 0) {
                        
                        results.temp <- as.data.frame(matrix(ncol = 25, nrow = 1))
                        colnames(results.temp) <- c("SpeciesID", "2010_SI", "2010_SI2", "2018_SI",
                                                    "2018_SI2", "Increaser", "Pred_Keep", "Total_Agriculture",
                                                    "Total_Energy", "Total_Forestry", "Total_Misc", "Total_RuralUrban",
                                                    "Total_Transportation", "UnderHF_Agriculture",
                                                    "UnderHF_Energy", "UnderHF_Forestry", "UnderHF_Misc", "UnderHF_RuralUrban",
                                                    "UnderHF_Transportation", "Area_Agriculture",
                                                    "Area_Energy", "Area_Forestry", "Area_Misc", "Area_RuralUrban",
                                                    "Area_Transportation")
                        results.temp[, "SpeciesID"] <- species
                        results.out <- rbind(results.out, results.temp)
                        next()
                        
                }
                
                # Subset to the region of interest and define temporary data frame
                sector.effect <- sector.effect[region.subset, ]
                results.temp[species, "SpeciesID"] <- species
                
                # Calculate intactness based on 2010 HFI
                results.temp[species, "2010_SI"] <- ifelse(sum(sector.effect[, "Cur_2010"]) >= sum(sector.effect[, "Ref"]), 
                                                           100 * (sum(sector.effect[, "Ref"]) / sum(sector.effect[, "Cur_2010"])), 
                                                           100 * (sum(sector.effect[, "Cur_2010"]) / sum(sector.effect[, "Ref"])))
                
                results.temp[species, "2010_SI2"] <- ifelse(100 * (sum(sector.effect[, "Cur_2010"]) / sum(sector.effect[, "Ref"])) >= 100, 
                                                            200 - results.temp[species, "2010_SI"], 
                                                            100 * (sum(sector.effect[, "Cur_2010"]) / sum(sector.effect[, "Ref"])))
                
                rm(sector.effect, cur.names, ref.names)
                
                ###################
                # Intactness 2018 #
                ###################
                
                # Load the provincial results and identify maximum abundance
                load.name <- species.list.2018[grep(as.character(species), species.list.2018)]
                load.name <- load.name[nchar(load.name) == min(nchar(load.name))]
                load(load.name)
                rm(load.name)
                
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
                
                rm(north.sector, south.sector)
                gc()
                
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
        results.out <- results.out[, c("SpeciesID", "ScientificName", "CommonNameIC", "TSNID", "Taxon",
                                       "Model_north", "Model_south", "Native", "OF_Birds",
                                       "OF_Species_ALPAC", "COSEWIC_Status", "SARA_Schedule",
                                       "AB_ESCC_2010", "AEP_GS_2015", "Pred_Keep", "Species_Inclusion",
                                       "Habitat_Assoc", "Increaser", "2010_SI", "2010_SI2", "2018_SI", "2018_SI2",
                                       area.names, regional.names, under.hf.names)]
        
        return(results.out)
        
}

######################
# Regional reporting # Creates sector effect information for the two modeling regions (forested vs prairie)
######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

regional_reports <- function(species.list, species.lookup, summary.region, overlap.region) {
            
            # Define blank results
            north.results <- NULL
            south.results <- NULL
            results.out <- list()
            
            for( species in species.lookup$SpeciesID ) {
                        
                        # Load the provincial results and identify maximum abundance
                        load.name <- species.list[grep(as.character(species), species.list)]
                        load.name <- load.name[nchar(load.name) == min(nchar(load.name))]
                        load(load.name)
                        rm(load.name)
                        
                        # Define north and south 
                        north.temp <- NULL
                        south.temp <- NULL
                        model.region <- NULL
                        
                        # If north and south regions are available, merge. Otherwise, use what is available
                        if(!is.null(north.sector)) {
                                    
                                    sector.effect <- north.sector
                                    model.region <- "N"
                                    
                        }
                        
                        if(!is.null(south.sector)) {
                                    
                                    sector.effect <- south.sector
                                    model.region <- "S"
                                    
                        }
                        
                        if(!(is.null(south.sector)) & !(is.null(south.sector))) {
                                    
                                    model.region <- "B"
                                    
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
                        
                        rm(north.sector, south.sector)
                        gc()
                        
                        # Merge predictions with sector areas
                        sector.effect <- as.data.frame(sector.effect)
                        sector.effect$LinkID <- rownames(sector.effect)
                        rownames(sector.effect) <- sector.effect$LinkID
                        
                        # Identify column names for the current, reference, and area groups
                        cur.names <- colnames(sector.effect)[grep("Cur_", colnames(sector.effect))]
                        ref.names <- colnames(sector.effect)[grep("Ref_", colnames(sector.effect))]
                        
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
                        
                        #################
                        # Sector Effect #
                        #################
                        
                        # North
                        if(model.region %in% c("B", "N")) {
                                    
                                    # Create vector for storing results
                                    north.temp <- data.frame(SpeciesID = species)
                                    
                                    # Filter to summary region
                                    region.id <- summary.region[summary.region$Regions == "Forested", "LinkID"]
                                    region.id <- region.id[region.id %in% rownames(sector.effect)] # match IDs 
                                    forested.sector <- sector.effect[region.id, ]
                                    
                                    # Sector effect information is already loaded
                                    # Regional
                                    regional.names <- NULL # Store the regional names for organizing
                                    for (se.name in cur.names[-grep("Native", cur.names)]) {
                                                
                                                update.name <- gsub("Cur", "Total", se.name)
                                                north.temp[update.name] <- ((sum(forested.sector[, se.name]) - sum(forested.sector[, gsub("Cur", "Ref", se.name)])) / sum(forested.sector[, ref.names])) * 100
                                                regional.names <- c(regional.names, update.name)
                                                rm(update.name)
                                                
                                    }
                                    
                                    # Under HF
                                    under.hf.names <- NULL # Store the under HF names for organizing
                                    for (se.name in cur.names[-grep("Native", cur.names)]) {
                                                
                                                update.name <- gsub("Cur", "UnderHF", se.name)
                                                north.temp[update.name] <- ((sum(forested.sector[, se.name]) - sum(forested.sector[, gsub("Cur", "Ref", se.name)])) / sum(forested.sector[, gsub("Cur", "Ref", se.name)])) * 100
                                                under.hf.names <- c(under.hf.names, update.name)
                                                rm(update.name)
                                                
                                    }
                                    
                        }
                        
                        # South
                        if(model.region %in% c("B", "S")) {
                                    
                                    # Create vector for storing results
                                    south.temp <- data.frame(SpeciesID = species)
                                    
                                    # Filter to summary region
                                    region.id <- summary.region[summary.region$Regions == "Prairie", "LinkID"]
                                    region.id <- region.id[region.id %in% rownames(sector.effect)] # match IDs 
                                    prairie.sector <- sector.effect[region.id, ]
                                    
                                    # Sector effect information is already loaded
                                    # Regional
                                    regional.names <- NULL # Store the regional names for organizing
                                    for (se.name in cur.names[-grep("Native", cur.names)]) {
                                                
                                                update.name <- gsub("Cur", "Total", se.name)
                                                south.temp[update.name] <- ((sum(prairie.sector[, se.name]) - sum(prairie.sector[, gsub("Cur", "Ref", se.name)])) / sum(prairie.sector[, ref.names])) * 100
                                                regional.names <- c(regional.names, update.name)
                                                rm(update.name)
                                                
                                    }
                                    
                                    # Under HF
                                    under.hf.names <- NULL # Store the under HF names for organizing
                                    for (se.name in cur.names[-grep("Native", cur.names)]) {
                                                
                                                update.name <- gsub("Cur", "UnderHF", se.name)
                                                south.temp[update.name] <- ((sum(prairie.sector[, se.name]) - sum(prairie.sector[, gsub("Cur", "Ref", se.name)])) / sum(prairie.sector[, gsub("Cur", "Ref", se.name)])) * 100
                                                under.hf.names <- c(under.hf.names, update.name)
                                                rm(update.name)
                                                
                                    }
                                    
                        }
                        
                        
                        # Combine species results
                        north.results <- rbind(north.results, north.temp)
                        south.results <- rbind(south.results, south.temp)
                        
                        # Store the results
                        results.out[["North"]] <- north.results
                        results.out[["South"]] <- south.results
                        
            }
            
            # Combine into dataframe
            results.out[["North"]] <- merge.data.frame(x = species.lookup, y = results.out[["North"]], by = "SpeciesID")
            results.out[["South"]] <- merge.data.frame(x = species.lookup, y = results.out[["South"]], by = "SpeciesID")
            
            return(results.out)
            
}

######################
# Regional reporting # Creates sector effect information for the two modeling regions (forested vs prairie)
######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# regional_reports_old <- function(species.list, species.lookup, summary.region) {
#         
#         # Define blank results
#         north.results <- NULL
#         south.results <- NULL
#         results.out <- list()
#         
#         for( species in species.lookup$SpeciesID ) {
#                 
#                 # Load the provincial results and identify maximum abundance
#                 load.name <- species.list[grep(as.character(species), species.list)]
#                 load.name <- load.name[nchar(load.name) == min(nchar(load.name))]
#                 load(load.name)
#                 rm(load.name)
#                 
#                 # Define north and south 
#                 north.temp <- NULL
#                 south.temp <- NULL
#                 
#                 # Northern predictions
#                 if(!is.null(north.sector)) {
#                         
#                         # Identify column names for the current, reference, and area groups
#                         cur.names <- colnames(north.sector)[grep("Cur_", colnames(north.sector))]
#                         ref.names <- colnames(north.sector)[grep("Ref_", colnames(north.sector))]
#                         
#                         # Calculate Ref and Current sums
#                         north.sector <- cbind(north.sector, Ref = as.numeric(rowSums(north.sector[, ref.names])))
#                         north.sector <- cbind(north.sector, Cur = as.numeric(rowSums(north.sector[, cur.names])))
#                         
#                         # Truncate to 99%
#                         q <- min(quantile(north.sector[, "Ref"], 0.999), quantile(north.sector[, "Cur"], 0.999))
#                         total.ref <- north.sector[, "Ref"]
#                         total.cur <- north.sector[, "Cur"]
#                         
#                         for (se.name in ref.names) {
#                                 
#                                 north.sector[, se.name] <- ifelse(total.ref >= q , (north.sector[, se.name] * (q/total.ref)), north.sector[, se.name])
#                                 
#                         }
#                         
#                         for (se.name in cur.names) {
#                                 
#                                 north.sector[, se.name] <- ifelse(total.cur >= q, (north.sector[, se.name] * (q/total.cur)), north.sector[, se.name])
#                                 
#                         }
#                         
#                         north.sector[, "Cur_Native"] <- 0.5 * (north.sector[, "Cur_Native"] + north.sector[, "Ref_Native"])
#                         north.sector[, "Ref_Native"] <- north.sector[, "Cur_Native"] 
#                         
#                         # Filter to summary region
#                         region.id <- summary.region[summary.region$Regions == "Forested", "LinkID"]
#                         region.id <- region.id[region.id %in% rownames(north.sector)] # match IDs 
#                         north.sector <- north.sector[region.id, ]
#                         
#                         ##################
#                         # Sector Effects #
#                         ##################
#                         
#                         # Create vector for storing results
#                         north.temp <- data.frame(SpeciesID = species)
#                         
#                         # Regional
#                         regional.names <- NULL # Store the regional names for organizing
#                         for (se.name in cur.names[-grep("Native", cur.names)]) {
#                                 
#                                 update.name <- gsub("Cur", "Total", se.name)
#                                 north.temp[update.name] <- ((sum(north.sector[, se.name]) - sum(north.sector[, gsub("Cur", "Ref", se.name)])) / sum(north.sector[, ref.names])) * 100
#                                 regional.names <- c(regional.names, update.name)
#                                 rm(update.name)
#                                 
#                         }
#                         
#                         # Under HF
#                         under.hf.names <- NULL # Store the under HF names for organizing
#                         for (se.name in cur.names[-grep("Native", cur.names)]) {
#                                 
#                                 update.name <- gsub("Cur", "UnderHF", se.name)
#                                 north.temp[update.name] <- ((sum(north.sector[, se.name]) - sum(north.sector[, gsub("Cur", "Ref", se.name)])) / sum(north.sector[, gsub("Cur", "Ref", se.name)])) * 100
#                                 under.hf.names <- c(under.hf.names, update.name)
#                                 rm(update.name)
#                                 
#                         }
#                         
#                         
#                 }
#                 
#                 # Southern predictions
#                 if(!is.null(south.sector)) {
#                         
#                         # Identify column names for the current, reference, and area groups
#                         cur.names <- colnames(south.sector)[grep("Cur_", colnames(south.sector))]
#                         ref.names <- colnames(south.sector)[grep("Ref_", colnames(south.sector))]
#                         
#                         # Calculate Ref and Current sums
#                         south.sector <- cbind(south.sector, Ref = as.numeric(rowSums(south.sector[, ref.names])))
#                         south.sector <- cbind(south.sector, Cur = as.numeric(rowSums(south.sector[, cur.names])))
#                         
#                         # Truncate to 99%
#                         q <- min(quantile(south.sector[, "Ref"], 0.999), quantile(south.sector[, "Cur"], 0.999))
#                         total.ref <- south.sector[, "Ref"]
#                         total.cur <- south.sector[, "Cur"]
#                         
#                         for (se.name in ref.names) {
#                                 
#                                 south.sector[, se.name] <- ifelse(total.ref >= q , (south.sector[, se.name] * (q/total.ref)), south.sector[, se.name])
#                                 
#                         }
#                         
#                         for (se.name in cur.names) {
#                                 
#                                 south.sector[, se.name] <- ifelse(total.cur >= q, (south.sector[, se.name] * (q/total.cur)), south.sector[, se.name])
#                                 
#                         }
#                         
#                         south.sector[, "Cur_Native"] <- 0.5 * (south.sector[, "Cur_Native"] + south.sector[, "Ref_Native"])
#                         south.sector[, "Ref_Native"] <- south.sector[, "Cur_Native"] 
#                         
#                         # Filter to summary region
#                         region.id <- summary.region[summary.region$Regions == "Prairie", "LinkID"]
#                         region.id <- region.id[region.id %in% rownames(south.sector)] # match IDs 
#                         south.sector <- south.sector[region.id, ]
#                         
#                         ##################
#                         # Sector Effects #
#                         ##################
#                         
#                         # Create vector for storing results
#                         south.temp <- data.frame(SpeciesID = species)
#                         
#                         # Regional
#                         regional.names <- NULL # Store the regional names for organizing
#                         for (se.name in cur.names[-grep("Native", cur.names)]) {
#                                 
#                                 update.name <- gsub("Cur", "Total", se.name)
#                                 south.temp[update.name] <- ((sum(south.sector[, se.name]) - sum(south.sector[, gsub("Cur", "Ref", se.name)])) / sum(south.sector[, ref.names])) * 100
#                                 regional.names <- c(regional.names, update.name)
#                                 rm(update.name)
#                                 
#                         }
#                         
#                         # Under HF
#                         under.hf.names <- NULL # Store the under HF names for organizing
#                         for (se.name in cur.names[-grep("Native", cur.names)]) {
#                                 
#                                 update.name <- gsub("Cur", "UnderHF", se.name)
#                                 south.temp[update.name] <- ((sum(south.sector[, se.name]) - sum(south.sector[, gsub("Cur", "Ref", se.name)])) / sum(south.sector[, gsub("Cur", "Ref", se.name)])) * 100
#                                 under.hf.names <- c(under.hf.names, update.name)
#                                 rm(update.name)
#                                 
#                         }
#                         
#                 }
#                 
#                 # Combine species results
#                 north.results <- rbind(north.results, north.temp)
#                 south.results <- rbind(south.results, south.temp)
#                 
#                 # Store the results
#                 results.out[["North"]] <- north.results
#                 results.out[["South"]] <- south.results
#                 
#         }
#         
#         # Combine into dataframe
#         results.out[["North"]] <- merge.data.frame(x = species.lookup, y = results.out[["North"]], by = "SpeciesID")
#         results.out[["South"]] <- merge.data.frame(x = species.lookup, y = results.out[["South"]], by = "SpeciesID")
#         
#         return(results.out)
#         
# }

