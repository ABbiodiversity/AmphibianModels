#
# Title: Functions to execute species habitat models
# Created: June 13th, 2019
# Last Updated: January 24th, 2022
# Author: Brandon Allen
# Objectives: Soil and land facet habitat model functions
# Keywords: Presence/Absence models 3.0
# Notes: 

###############################
# Presence/Absence models 3.0 # This version required time of year and temperature information for estimating the effects of detection.
###############################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pa_models_3.0 <- function (data.analysis, results.store, landscape.models, prediction.matrix, region, age.class, species.ID, weight.method, detection, coef.adjust, boot) {
        
        #
        # 0.0 Cleaning
        #
        
        # Remove stations that were not sampled within the seasonality window 
        data.analysis <- data.analysis[!is.na(data.analysis[, species.ID]), ]
        data.analysis["nQuadrant"] <- 1
        
        # Create new column which is the observed count data (used for weighting in the age relationships), one for the converted probabilities
        data.analysis["Count"] <- as.integer(ifelse(data.analysis[, species.ID] > 0, 1, 0))
        data.analysis["pcount"] <- as.integer(ifelse(data.analysis[, species.ID] > 0, 1, 0) / data.analysis$nQuadrant)
        
        # Since we are treating each recording individually, we fix visit to 1 as we no longer need to apply this weight.
        data.analysis$visit <- 1 
        
        #
        # 0.1 Bootstrap
        # If the boot index is 1, use the complete data set, otherwise perform subsampling with replacement based on block ID
        # Added catch to maintain sample size greater than 20
        #
        
        if (boot != 1) {
                
                unique.blocks <- unique(data.analysis$Block)
                n.site <- 0
                
                while (n.site < 20) {
                        
                        data.in <- NULL
                        for (sample_block in 1:length(unique.blocks)) {
                                
                                data.block <- data.analysis[grep(unique.blocks[sample_block], data.analysis$Block), ]
                                data.block <- data.block[sample(1:nrow(data.block), nrow(data.block), replace=TRUE), ]
                                data.in <- rbind(data.in, data.block)
                                
                        }
                        n.site <- sum(ifelse(data.in[, "Count"] > 0, 1, 0))
                        
                }
                
        } else {
                
                data.in <- data.analysis
                
        }
        
        #
        # 1. Landscape models
        # Loop through each model in the model list, calculate AIC scores, and weight
        #
        
        space.climate.store <- landscape.store <- list(NULL) 
        
        for (model in 1:length(landscape.models)) {
                
                landscape.store[[model]] <- try(glm(landscape.models[[model]], family = "binomial", 
                                                    data = data.in, 
                                                    maxit = 250))
                
        }
        
        # 1.01 Time of year adjustment
        
        if (detection == TRUE) {
                
                for (model in 1:length(landscape.models)) {
                        
                        landscape.store[[model]] <- try(update(landscape.store[[model]], .~.+ ToY + ToY2 + Hourly_Temp))
                        
                }
                
        } 
        
        nModels <- length(landscape.store)
        
        # AIC calculation  (I'm using AIC here, because this is primarily for prediction, rather than finding a minimial best model)
        aic.ta <- rep(999999999, (nModels))
        
        for (i in 1:(nModels)) {
                
                if (!is.null(landscape.store[[i]]) & class(landscape.store[[i]])[1] != "try-error") {  # last part is to not used non-converged models, unless none converged
                        aic.ta[i] <- AICc(landscape.store[[i]])
                }
                
        }
        
        aic.delta <- aic.ta - min(aic.ta)
        aic.exp <- exp(-1 / 2 * aic.delta)
        aic.wt.ta <- aic.exp / sum(aic.exp)
        
        # # MODEL CHECKING
        # # Base model check (AUC)
        # for (model in 1:length(landscape.models)) {
        # 
        #         print(auc(ifelse(data.in[, "Count"] > 0, 1, 0), predict(landscape.store[[model]])))
        # 
        # }
        
        # 1.02 Add ToY, ToY2, and Hourly_Temp columns/rows to the prediction matrix. 
        # Default the value to 0, then update by pulling the coefficient/se from the model output directly
        
        if (detection == TRUE) {
                
                prediction.matrix[, "ToY"] <- 0
                prediction.matrix[, "ToY2"] <- 0
                prediction.matrix[, "Hourly_Temp"] <- 0
                
                prediction.matrix["ToY", ] <- 0
                prediction.matrix["ToY2", ] <- 0
                prediction.matrix["Hourly_Temp", ] <- 0
                
                prediction.matrix["ToY", "ToY"] <- 1
                prediction.matrix["ToY2", "ToY2"] <- 1
                prediction.matrix["Hourly_Temp", "Hourly_Temp"] <- 1
                

        } 
        
        # 1.1 Store coefficients
        
        if (region == "North") {
                
                p1 <- p1.se <- array(0, c(nModels, nrow(prediction.matrix))) # Predictions for each model for each soil and HF type type.  These are the main coefficients
                p.site1 <- p.site1.se <- array(0, c(nModels, nrow(data.in)))  # Predictions for each model and each site.  These are used as the offsets in the next stage
                
                for (i in 1:nModels) {
                        
                        if (class(landscape.store[[i]])[1] != "try-error") {   # Prediction is 0 if model failed, but this is not used because AIC wt would equal 0
                                
                                p <- predict(landscape.store[[i]], newdata = data.frame(prediction.matrix), se.fit = TRUE)  # For each type.  Predictions made at 0% Aspen.  All predictions made with new protocol.  Aspen effect added later, and plotted as separate points
                                p1[i,] <- p$fit
                                p1.se[i,] <- p$se.fit
                                
                                # Update non landcover coefficients
                                if (detection == TRUE) {

                                        p1[i, c((ncol(p1) - 2):ncol(p1))] <- summary(landscape.store[[i]])$coefficients[c("ToY", "ToY2", "Hourly_Temp"), 1]
                                        p1.se[i, c((ncol(p1) - 2):ncol(p1))] <- summary(landscape.store[[i]])$coefficients[c("ToY", "ToY2", "Hourly_Temp"), 2]


                                }

                                model.fit <- predict(landscape.store[[i]], se.fit = TRUE)
                                p.site1[i,] <- model.fit$fit  # For each site (using the original d.sp data frame)
                                p.site1.se[i,] <- model.fit$se.fit  # For each site SE
                        }
                        
                        rm(p, model.fit)
                        
                }
                
                # Choose if the AIC or Inverse variance method is used for weighting coefficients
                if (weight.method == "AIC") {
                        
                        tTypeMean <- tTypeVar <- NULL  # Logit-scaled model averages for each veg type
                        
                        aic.wt.ta.adj <- ifelse(aic.wt.ta < 0.01, 0, aic.wt.ta)  # Use adjusted weight to avoid poor fitting models with extreme values for certain types (usually numerically equivalent to 0's or 1's)
                        aic.wt.ta.adj <- aic.wt.ta.adj / sum(aic.wt.ta.adj)
                        
                        for (i in 1:nrow(prediction.matrix)) {
                                
                                tTypeMean[i] <- sum(aic.wt.ta.adj * p1[, i])  # Mean of nModels for each veg type, on transformed scale
                                tTypeVar[i] <- sum(aic.wt.ta.adj * sqrt(p1.se[, i]^2 + (rep(tTypeMean[i], nModels) - p1[, i])^2))^2  # AIC-weighted variance of mean...
                                
                        }
                        
                        names(tTypeMean) <- names(tTypeVar) <- rownames(prediction.matrix) # Rename the final intercept coefficient to productive
                        data.in$prediction <- colSums(p.site1 * aic.wt.ta.adj)  # Add logit-scaled model average prediction to data frame
                        
                        # Put coefficients in Coef matrix (and SE's)
                        results.store$landscape.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))] <- plogis(tTypeMean[colnames(results.store$landscape.coef)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))]])  # On ordinal scale
                        results.store$landscape.se[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))] <- sqrt(tTypeVar[colnames(results.store$landscape.se)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))]])  # On logit scale
                        
                        # If detection is present, store results
                        if (detection == TRUE) {
                                
                                results.store$detection.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))] <- tTypeMean[colnames(results.store$detection.coef)[na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))]]  # On ordinal scale
                                
                        }

                } 
                
                if (weight.method == "IVW") {
                        
                        tTypeMean <- tTypeVar <- NULL   # Logit-scaled IVW for each veg type
                        
                        mod.converged <- NULL
                        for (i in 1:nModels) {mod.converged[i] <- landscape.store[[i]]$converged}
                        
                        for (i in 1:nrow(prediction.matrix)) {
                                
                                tTypeMean[i] <- sum(p1[mod.converged ,i] / p1.se[mod.converged ,i]^2 )  / sum(1/p1.se[mod.converged ,i]^2 )# IVW mean of  nModels for each veg type, on transformed scale
                                tTypeVar[i]<-  1/sum(1/p1.se[mod.converged ,i]^2 ) # IVWvariance of mean...
                                
                        }
                        
                        names(tTypeMean) <- names(tTypeVar) <- rownames(prediction.matrix)
                        data.in$prediction <- colSums(p.site1/ p.site1.se^2) / colSums(1/p.site1.se^2)  # Add logit-scaled model average prediction to data frame
                        
                        # Put coefficients in Coef matrix (and SE's)
                        results.store$landscape.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))] <- plogis(tTypeMean[colnames(results.store$landscape.coef)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))]])  # On ordinal scale
                        results.store$landscape.se[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))] <- sqrt(tTypeVar[colnames(results.store$landscape.se)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))]])  # On logit scale
                        
                        # If detection is present, store results
                        if (detection == TRUE) {
                                
                                results.store$detection.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))] <-tTypeMean[colnames(results.store$detection.coef)[na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))]]  # On ordinal scale
                                
                        }
                        
                } 
                
                
        }
        
        if (region == "South") {
                
                # South model needs to account for the addition of the probability of aspen coefficient, otherwise the processes is similar. 
                p1 <- p1.se <- array(0, c(nModels, nrow(prediction.matrix))) # Predictions for each model for each soil and HF type type.  These are the main coefficients
                p.site1 <- p.site1.se <- array(0, c(nModels, nrow(data.in)))  # Predictions for each model and each site.  These are used as the offsets in the next stage
                p1.aspen <- p1.se.aspen <- rep(0, nModels)  # Predictions for each model for aspen coefficient
                
                for (i in 1:nModels) {
                        
                        if (class(landscape.store[[i]])[1] != "try-error") {   # Prediction is 0 if model failed, but this is not used because AIC wt would equal 0
                                
                                p <- predict(landscape.store[[i]], newdata = data.frame(prediction.matrix, paspen = 0), se.fit = TRUE)  # For each type.  Predictions made at 0% Aspen.  All predictions made with new protocol.  Aspen effect added later, and plotted as separate points
                                p1[i,] <- p$fit
                                p1.se[i,] <- p$se.fit
                                
                                # Update non landcover coefficients
                                if (detection == TRUE) {
                                        
                                        p1[i, c((ncol(p1) - 2):ncol(p1))] <- summary(landscape.store[[i]])$coefficients[c("ToY", "ToY2", "Hourly_Temp"), 1]
                                        p1.se[i, c((ncol(p1) - 2):ncol(p1))] <- summary(landscape.store[[i]])$coefficients[c("ToY", "ToY2", "Hourly_Temp"), 2]
                                        
                                        
                                } 
                                
                                model.fit <- predict(landscape.store[[i]], se.fit = TRUE)
                                p.site1[i,] <- model.fit$fit  # For each site (using the original d.sp data frame)
                                p.site1.se[i,] <- model.fit$se.fit  # For each site SE
                                
                                if ("paspen" %in% names(coef(landscape.store[[i]]))) {
                                        p1.aspen[i] <- coefficients(summary(landscape.store[[i]]))["paspen","Estimate"]
                                        p1.se.aspen[i] <- coefficients(summary(landscape.store[[i]]))["paspen","Std. Error"]
                                }
                                
                                
                        }
                        
                        rm(p, model.fit)
                        
                }
                
                # Choose if the AIC or Inverse variance method is used for weighting coefficients
                if (weight.method == "AIC") { 
                        
                        tTypeMean <- tTypeVar <- NULL  # Logit-scaled model averages for each veg type
                        
                        aic.wt.ta.adj <- ifelse(aic.wt.ta < 0.01, 0, aic.wt.ta)  # Use adjusted weight to avoid poor fitting models with extreme values for certain types (usually numerically equivalent to 0's or 1's)
                        aic.wt.ta.adj <- aic.wt.ta.adj / sum(aic.wt.ta.adj)
                        
                        tpAspenMean <- sum(aic.wt.ta * p1.aspen)
                        tpAspenVar <- sum(aic.wt.ta * sqrt(p1.se.aspen^2 + (rep(tpAspenMean, nModels) - p1.aspen)^2))^2  # AIC-weighted variance of mean...
                        
                        for (i in 1:nrow(prediction.matrix)) {
                                
                                tTypeMean[i] <- sum(aic.wt.ta.adj * p1[, i])  # Mean of nModels for each veg type, on transformed scale
                                tTypeVar[i] <- sum(aic.wt.ta.adj * sqrt(p1.se[, i]^2 + (rep(tTypeMean[i], nModels) - p1[, i])^2))^2  # AIC-weighted variance of mean...
                                
                        }
                        
                        names(tTypeMean) <- names(tTypeVar)<- rownames(prediction.matrix) # Rename the final intercept coefficient to productive
                        
                        data.in$prediction <- colSums(p.site1 * aic.wt.ta.adj)  # Add logit-scaled model average prediction to data frame
                        
                        # Put coefficients in Coef matrix (and SE's)
                        results.store$landscape.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))] <- plogis(tTypeMean[colnames(results.store$landscape.coef)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))]])  # On ordinal scale
                        results.store$landscape.se[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))] <- sqrt(tTypeVar[colnames(results.store$landscape.se)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))]])  # On logit scale
                        results.store$paspen.coef[boot, "paspen"] <- tpAspenMean
                        results.store$paspen.se[boot, "paspen"] <- sqrt(tpAspenVar)
                        
                        # If detection is present, store results
                        if (detection == TRUE) {
                                
                                results.store$detection.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))] <- tTypeMean[colnames(results.store$detection.coef)[na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))]]  # On ordinal scale
                                
                        }
                        
                }
                
                if (weight.method == "IVW") { 
                        
                        tTypeMean <- tTypeVar <- NULL   # Logit-scaled IVW for each veg type
                        
                        mod.converged <- NULL
                        
                        # If paspen coefficients / SE is 0, remove
                        p1.aspen <- p1.aspen[p1.aspen > 0]
                        p1.se.aspen <- p1.se.aspen[p1.se.aspen > 0]
                        
                        tpAspenMean <- sum(p1.aspen / p1.se.aspen^2) / sum(1 / p1.se.aspen^2)
                        tpAspenVar <-  1/sum(1/p1.se.aspen^2 ) # IVWvariance of mean...
                        
                        for (i in 1:nModels) {mod.converged[i] <- landscape.store[[i]]$converged} # Only use models that have converged
                        
                        for (i in 1:nrow(prediction.matrix)) {
                                
                                tTypeMean[i]<- sum(p1[mod.converged ,i] / p1.se[mod.converged ,i]^2 )  / sum(1/p1.se[mod.converged ,i]^2 )# IVW mean of  nModels for each veg type, on transformed scale
                                tTypeVar[i]<-  1/sum(1/p1.se[mod.converged ,i]^2 ) # IVWvariance of mean...
                                
                        }
                        
                        names(tTypeMean) <- names(tTypeVar) <- rownames(prediction.matrix)
                        data.in$prediction <-  data.in$prediction <- colSums(p.site1/ p.site1.se^2) / colSums(1/p.site1.se^2) # Add logit-scaled model average prediction to data frame
                        
                        # Put coefficients in Coef matrix (and SE's)
                        results.store$landscape.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))] <- plogis(tTypeMean[colnames(results.store$landscape.coef)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))]])  # On ordinal scale
                        results.store$landscape.se[boot, na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))] <- sqrt(tTypeVar[colnames(results.store$landscape.se)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))]])  # On logit scale
                        results.store$paspen.coef[boot, "paspen"] <- tpAspenMean
                        results.store$paspen.se[boot, "paspen"] <- sqrt(tpAspenVar)
                        
                        # If detection is present, store results
                        if (detection == TRUE) {
                                
                                results.store$detection.coef[boot, na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))] <- tTypeMean[colnames(results.store$detection.coef)[na.omit(match(names(tTypeMean), colnames(results.store$detection.coef)))]]  
                                
                        }
                        
                }
                
        }
        
        # # TESTING CREATION OF PREDICTION THROUGH PREDICT VERSUS MANUAL
        # # Predict function
        # auc(ifelse(data.in[, "Count"] > 0, 1, 0), data.in$prediction)
        # 
        # # Manual function
        # # Landcover
        # p1 <- colSums(plogis(tTypeMean[-c(45,46,47)]) * t(data.in[ ,names(tTypeMean)[-c(45,46,47)]]))  # Multiply coefficients by landscape type proportion for each site to get offset
        # 
        # # Add the detection information
        # p1 <- plogis(qlogis(p1) + colSums(tTypeMean[c(45,46,47)] * t(data.in[ ,names(tTypeMean)[c(45,46,47)]])))
        # 
        # p1 <- colSums(tTypeMean * t(data.in[ ,names(tTypeMean)]))  # Multiply coefficients by landscape type proportion for each site to get offset
        # 
        # 
        # auc(ifelse(data.in[, "Count"] > 0, 1, 0), p1)
        
        # 
        # 2.0 Run models for age within each stand type
        #
        
        if (region == "North" & age.class == TRUE) {
                
                # 2.1. Set up separate dataframes for sites containing a minimum amount of each broad stand type
                stand.data <- list(WhiteSpruce = NULL, Pine = NULL, 
                                   Deciduous = NULL, Mixedwood = NULL, 
                                   BlackSpruce = NULL)
                
                cutoff.pc.for.age <- 0.1  # Set the cut-off for the minimum proportion of a stand type for a site to be included in the age analysis for a stand type
                pSoftLin <- c(0.049, 0.0893, 0.434, 0.396)  # Values based on overlap of softlin and the stand types in the 1km summary file.  Lowland spruce includes other wetlands. Updated proportions 2020-11-17
                pSoftLin <- pSoftLin / sum(pSoftLin)  # Make sure they sum to zero
                aic.age <- array(NA, c(1, 3)) 
                stand.name <- c("WhiteSpruce", "Pine", "Deciduous", "Mixedwood", "BlackSpruce")
                
                for (i in 0:8) {  # Add sites with each age class of the stand type to a separate data frames for each stand type
                        
                        i1 <- ifelse(i == 0, "R", i)  # For variable name
                        i2 <- ifelse(i == 0, 0.5, i)  # For twenty-year age
                        
                        for (stand.id in 1:length(stand.name)) {
                                
                                temp.name <- paste(stand.name[stand.id], i1, sep = "")  # Col name
                                if (sum(data.in[ ,temp.name] > cutoff.pc.for.age) > 0) stand.data[[stand.name[stand.id]]] <- rbind(stand.data[[stand.name[stand.id]]], data.frame(pCount = data.in[data.in[, temp.name] > cutoff.pc.for.age, "Count"] / data.in$nQuadrant[data.in[, temp.name] > cutoff.pc.for.age], age = i2, wt1 = data.in[data.in[, temp.name] > cutoff.pc.for.age, temp.name]*data.in$visit[data.in[, temp.name] > cutoff.pc.for.age]*data.in$nQuadrant[data.in[, temp.name] > cutoff.pc.for.age],  p = data.in$prediction[data.in[, temp.name] > cutoff.pc.for.age]))  # Weight is the proportion of the site of that age class and stand type, multiplied by the original weight (which accounts for revisited sites) and the number of binomial trials
                                
                        }
                        
                }
                
                rm(stand.name)
                
                # For upland forests, add grass sites as age 1, but with half the normal weight
                d.grass <- data.in[data.in[, "Grass"] > cutoff.pc.for.age, ]  # Select sites with grass (to make lines below less unruly)
                pWhiteSpruce <- d.grass$WhiteSpruce / (d.grass$WhiteSpruce + d.grass$Pine + d.grass$Deciduous + d.grass$Mixedwood + d.grass$BlackSpruce + 0.01)
                pPine <- d.grass$Pine / (d.grass$WhiteSpruce + d.grass$Pine + d.grass$Deciduous + d.grass$Mixedwood + d.grass$BlackSpruce + 0.01)
                pDeciduous <- d.grass$Deciduous / (d.grass$WhiteSpruce + d.grass$Pine + d.grass$Deciduous + d.grass$Mixedwood + d.grass$BlackSpruce + 0.01)
                pMixedwood <- d.grass$WhiteSpruce / (d.grass$WhiteSpruce + d.grass$Pine + d.grass$Deciduous + d.grass$Mixedwood + d.grass$BlackSpruce + 0.01)
                
                # If there are no sites to add, don't bind
                if(sum(pWhiteSpruce) != 0) {stand.data$WhiteSpruce <- rbind(stand.data$WhiteSpruce, data.frame(pCount = d.grass[pWhiteSpruce > 0, "Count"] / d.grass$nQuadrant[pWhiteSpruce > 0], age = 0.5, wt1 = d.grass[pWhiteSpruce > 0, "Grass"]/2*d.grass$visit[pWhiteSpruce > 0]*d.grass$nQuadrant[pWhiteSpruce > 0]*pWhiteSpruce[pWhiteSpruce > 0], p = d.grass$prediction[pWhiteSpruce > 0]))}
                if(sum(pPine) != 0) {stand.data$Pine <- rbind(stand.data$Pine, data.frame(pCount = d.grass[pPine > 0, "Count"] / d.grass$nQuadrant[pPine>0], age = 0.5, wt1 = d.grass[pPine>0,"Grass"]/2*d.grass$visit[pPine>0]*d.grass$nQuadrant[pPine>0]*pPine[pPine>0], p = d.grass$prediction[pPine>0]))}
                if(sum(pDeciduous) != 0) {stand.data$Deciduous <- rbind(stand.data$Deciduous, data.frame(pCount = d.grass[pDeciduous > 0, "Count"] / d.grass$nQuadrant[pDeciduous > 0], age = 0.5, wt1 = d.grass[pDeciduous > 0, "Grass"] / 2*d.grass$visit[pDeciduous > 0]*d.grass$nQuadrant[pDeciduous > 0]*pDeciduous[pDeciduous > 0], p = d.grass$prediction[pDeciduous > 0]))}
                if(sum(pMixedwood) != 0) {stand.data$Mixedwood <- rbind(stand.data$Mixedwood, data.frame(pCount = d.grass[pMixedwood > 0, "Count"] / d.grass$nQuadrant[pMixedwood > 0], age = 0.5, wt1 = d.grass[pMixedwood > 0, "Grass"] / 2*d.grass$visit[pMixedwood > 0]*d.grass$nQuadrant[pMixedwood > 0]*pMixedwood[pMixedwood > 0], p = d.grass$prediction[pMixedwood > 0]))}
                
                # And same for sites with shrubs
                d.shrub <- data.in[data.in[, "Shrub"] > cutoff.pc.for.age,]  # Select sites with shrubs (to make lines below less unruly)
                pWhiteSpruce <- d.shrub$WhiteSpruce / (d.shrub$WhiteSpruce + d.shrub$Pine + d.shrub$Deciduous + d.shrub$Mixedwood + d.shrub$BlackSpruce + 0.01)
                pPine <- d.shrub$Pine / (d.shrub$WhiteSpruce + d.shrub$Pine + d.shrub$Deciduous + d.shrub$Mixedwood + d.shrub$BlackSpruce + 0.01)
                pDeciduous <- d.shrub$Deciduous / (d.shrub$WhiteSpruce + d.shrub$Pine + d.shrub$Deciduous + d.shrub$Mixedwood + d.shrub$BlackSpruce + 0.01)
                pMixedwood <- d.shrub$WhiteSpruce / (d.shrub$WhiteSpruce + d.shrub$Pine + d.shrub$Deciduous + d.shrub$Mixedwood + d.shrub$BlackSpruce + 0.01)
                
                if(sum(pWhiteSpruce) != 0) {stand.data$WhiteSpruce <- rbind(stand.data$WhiteSpruce, data.frame(pCount = d.shrub[pWhiteSpruce > 0 , "Count"] / d.shrub$nQuadrant[pWhiteSpruce > 0], age = 1, wt1 = d.shrub[pWhiteSpruce > 0, "Shrub"] / 2*d.shrub$visit[pWhiteSpruce > 0]*d.shrub$nQuadrant[pWhiteSpruce > 0]*pWhiteSpruce[pWhiteSpruce > 0], p = d.shrub$prediction[pWhiteSpruce > 0]))}
                if(sum(pPine) != 0) {stand.data$Pine <- rbind(stand.data$Pine, data.frame(pCount = d.shrub[pPine > 0, "Count"] / d.shrub$nQuadrant[pPine > 0], age = 1, wt1 = d.shrub[pPine > 0, "Shrub"] / 2*d.shrub$visit[pPine > 0]*d.shrub$nQuadrant[pPine > 0]*pPine[pPine > 0], p = d.shrub$prediction[pPine > 0]))}
                if(sum(pDeciduous) != 0) {stand.data$Deciduous <- rbind(stand.data$Deciduous, data.frame(pCount = d.shrub[pDeciduous > 0, "Count"] / d.shrub$nQuadrant[pDeciduous > 0], age = 1, wt1 = d.shrub[pDeciduous > 0, "Shrub"] / 2*d.shrub$visit[pDeciduous > 0]*d.shrub$nQuadrant[pDeciduous > 0]*pDeciduous[pDeciduous > 0], p = d.shrub$prediction[pDeciduous > 0]))}
                if(sum(pMixedwood) != 0) {stand.data$Mixedwood <- rbind(stand.data$Mixedwood, data.frame(pCount = d.shrub[pMixedwood > 0, "Count"] / d.shrub$nQuadrant[pMixedwood > 0], age = 1, wt1 = d.shrub[pMixedwood > 0, "Shrub"] / 2*d.shrub$visit[pMixedwood > 0]*d.shrub$nQuadrant[pMixedwood > 0]*pMixedwood[pMixedwood > 0], p = d.shrub$prediction[pMixedwood > 0]))}
                
                # Combined data frames for models using more than one stand type
                stand.data$Conifpine <- rbind(cbind(stand.data$WhiteSpruce, Sp = "WhiteSpruce"), cbind(stand.data$Pine, Sp = "Pine"))
                stand.data$Decidmixed <- rbind(cbind(stand.data$Pine, Sp = "Decid"), cbind(stand.data$Pine, Sp = "Mixed"))
                stand.data$All <- rbind(cbind(stand.data$WhiteSpruce, Sp = "WhiteSpruce"), cbind(stand.data$Pine, Sp = "Pine"), cbind(stand.data$Deciduous, Sp = "Decid"), cbind(stand.data$Mixedwood, Sp = "Mixed"), cbind(stand.data$BlackSpruce, Sp = "BlackSpruce"))
                
                # 2.2. Then fit models of age functions
                m.age <- list(NULL)  # Age models for each of the 5 stand types plus 3 combinations
                if (sum(sign(stand.data$WhiteSpruce$pCount)) > 4) {
                        m.age[[1]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$WhiteSpruce, family = "binomial", weights = wt1)  # Fit spline through age data for that stand type if enough records
                } else {
                        m.age[[1]] <- gam(pCount~1 + offset(p), data = stand.data$WhiteSpruce, family = "binomial", weights = wt1)  # Fit constant only if too few records in that stand type
                }
                if (sum(sign(stand.data$Pine$pCount)) > 4) {
                        m.age[[2]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$Pine, family = "binomial", weights = wt1)
                } else {
                        m.age[[2]] <- gam(pCount~1 + offset(p), data = stand.data$Pine, family = "binomial", weights = wt1)
                }
                if (sum(sign(stand.data$Pine$pCount)) > 4) {
                        m.age[[3]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$Pine, family = "binomial", weights = wt1)
                } else {
                        m.age[[3]] <- gam(pCount~1 + offset(p), data = stand.data$Pine, family = "binomial", weights = wt1)
                }
                if (sum(sign(stand.data$Mixedwood$pCount)) > 4) {
                        m.age[[4]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$Mixedwood, family = "binomial", weights = wt1)
                } else {
                        m.age[[4]] <- gam(pCount~1 + offset(p), data = stand.data$Mixedwood, family = "binomial", weights = wt1)
                }
                if (sum(sign(stand.data$BlackSpruce$pCount)) > 4) {
                        m.age[[5]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$BlackSpruce, family = "binomial", weights = wt1)
                } else {
                        m.age[[5]] <- gam(pCount~1 + offset(p), data = stand.data$BlackSpruce, family = "binomial", weights = wt1)
                }
                if (sum(sign(stand.data$Conifpine$pCount)) > 4) {
                        m.age[[6]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$Conifpine, family = "binomial", weights = wt1)
                } else {
                        m.age[[6]] <- gam(pCount~1 + offset(p), data = stand.data$Conifpine, family = "binomial", weights = wt1)
                }
                if (sum(sign(stand.data$Decidmixed$pCount)) > 4) {
                        m.age[[7]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$Decidmixed, family = "binomial", weights = wt1)
                } else {
                        m.age[[7]] <- gam(pCount~1 + offset(p), data = stand.data$Decidmixed, family = "binomial", weights = wt1)
                }
                if (sum(sign(stand.data$All$pCount)) > 4) {
                        m.age[[8]] <- gam(pCount~s(sqrt(age),k = 3,m = 2) + offset(p), data = stand.data$All, family = "binomial", weights = wt1)
                } else {
                        m.age[[8]] <- gam(pCount~1 + offset(p), data = stand.data$All, family = "binomial", weights = wt1)
                }
                # AIC for options - each separate, pairwise combinations, combine all
                aic.age[1,1] <- AIC(m.age[[1]]) + AIC(m.age[[2]]) + AIC(m.age[[3]]) + AIC(m.age[[4]]) + AIC(m.age[[5]])-4  # Minus 4 for the 4 extra variances
                aic.age[1,2] <- AIC(m.age[[5]]) + AIC(m.age[[6]]) + AIC(m.age[[7]])-2
                aic.age[1,3] <- AIC(m.age[[8]])
                
                # 2.3. Predict coefficients (and logit-scale SE) for each age class for each stand type
                if (which.min(aic.age[1, ]) == 1) {
                        
                        p3 <- predict(m.age[[1]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["WhiteSpruce"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "WhiteSpruceR"):which(colnames(results.store$landscape.coef) == "WhiteSpruce8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "WhiteSpruceR"):which(colnames(results.store$landscape.coef) == "WhiteSpruce8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[2]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["Pine"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "PineR"):which(colnames(results.store$landscape.coef) == "Pine8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "PineR"):which(colnames(results.store$landscape.coef) == "Pine8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[3]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["Deciduous"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "DeciduousR"):which(colnames(results.store$landscape.coef) == "Deciduous8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "DeciduousR"):which(colnames(results.store$landscape.coef) == "Deciduous8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[4]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["Mixedwood"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "MixedwoodR"):which(colnames(results.store$landscape.coef) == "Mixedwood8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "MixedwoodR"):which(colnames(results.store$landscape.coef) == "Mixedwood8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[5]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["BlackSpruce"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "BlackSpruceR"):which(colnames(results.store$landscape.coef) == "BlackSpruce8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "BlackSpruceR"):which(colnames(results.store$landscape.coef) == "BlackSpruce8")] <- p3$se.fit  # SE still on the logit scale
                        
                }
                
                if (which.min(aic.age[1, ]) == 2) {
                        
                        p3 <- predict(m.age[[6]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["WhiteSpruce"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "WhiteSpruceR"):which(colnames(results.store$landscape.coef) == "WhiteSpruce8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "WhiteSpruceR"):which(colnames(results.store$landscape.coef) == "WhiteSpruce8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[6]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["Pine"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "PineR"):which(colnames(results.store$landscape.coef) == "Pine8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "PineR"):which(colnames(results.store$landscape.coef) == "Pine8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[7]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["Deciduous"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "DeciduousR"):which(colnames(results.store$landscape.coef) == "Deciduous8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "DeciduousR"):which(colnames(results.store$landscape.coef) == "Deciduous8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[7]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["Mixedwood"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "MixedwoodR"):which(colnames(results.store$landscape.coef) == "Mixedwood8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "MixedwoodR"):which(colnames(results.store$landscape.coef) == "Mixedwood8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[5]], newdata = data.frame(age = c(0.5,1:8), p = tTypeMean["BlackSpruce"]), se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "BlackSpruceR"):which(colnames(results.store$landscape.coef) == "BlackSpruce8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "BlackSpruceR"):which(colnames(results.store$landscape.coef) == "BlackSpruce8")] <- p3$se.fit  # SE still on the logit scale
                        
                }
                
                if (which.min(aic.age[1, ]) == 3) {
                        
                        p3 <- predict(m.age[[8]], newdata = data.frame(age = c(0.5,1:8),  p = tTypeMean["WhiteSpruce"]),  se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "WhiteSpruceR"):which(colnames(results.store$landscape.coef) == "WhiteSpruce8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot,which(colnames(results.store$landscape.coef) == "WhiteSpruceR"):which(colnames(results.store$landscape.coef) == "WhiteSpruce8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[8]], newdata = data.frame(age = c(0.5,1:8),  p = tTypeMean["Pine"]),  se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "PineR"):which(colnames(results.store$landscape.coef) == "Pine8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "PineR"):which(colnames(results.store$landscape.coef) == "Pine8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[8]], newdata = data.frame(age = c(0.5,1:8),  p = tTypeMean["Deciduous"]),  se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "DeciduousR"):which(colnames(results.store$landscape.coef) == "Deciduous8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "DeciduousR"):which(colnames(results.store$landscape.coef) == "Deciduous8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[8]], newdata = data.frame(age = c(0.5,1:8),  p = tTypeMean["Mixedwood"]),  se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "MixedwoodR"):which(colnames(results.store$landscape.coef) == "Mixedwood8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "MixedwoodR"):which(colnames(results.store$landscape.coef) == "Mixedwood8")] <- p3$se.fit  # SE still on the logit scale
                        
                        p3 <- predict(m.age[[8]], newdata = data.frame(age = c(0.5,1:8),  p = tTypeMean["BlackSpruce"]),  se.fit = T)
                        results.store$landscape.coef[boot, which(colnames(results.store$landscape.coef) == "BlackSpruceR"):which(colnames(results.store$landscape.coef) == "BlackSpruce8")] <- plogis(p3$fit)  # Coefficient on ordinal scale
                        results.store$landscape.se[boot, which(colnames(results.store$landscape.coef) == "BlackSpruceR"):which(colnames(results.store$landscape.coef) == "BlackSpruce8")] <- p3$se.fit  # SE still on the logit scale
                        
                }
                
                # 2.4 Adjust SE for estimates that are numerically equivalent to 0 - use exact binomial CI
                i4 <- which(is.na(results.store$landscape.coef[boot, ]) | results.store$landscape.coef[boot, ] < 0.0001)
                results.store$landscape.coef[boot, i4] <- 0.0001
                ntrials <- (colSums(data.in[, match(colnames(results.store$landscape.coef), names(data.in))])*2) + 1  # *2 because many rare types (which have the 0's are estimated from groups of types),  + 1 to avoid failure for old cutblocks, 
                for (i5 in i4) results.store$landscape.se[boot, i5] <-  (qlogis(binom.confint(0, ntrials[i5], conf.level = 0.9)[5,6]) - qlogis(0.0001)) / 1.65 # That conf/level gives 1SE.  Value is on logit scale, assuming estimate of 0.0001.  This is set up so that when 90% CI's are ploted below, they will return the exact 90% binom confints
                
                # 2.5 Converge CC to natural trajectory
                # Smooth CC trajectories into natural ones
                results.store$landscape.coef[boot, "CCWhiteSpruce2"] <- results.store$landscape.coef[boot, "WhiteSpruce2"]*0.50 + results.store$landscape.coef[boot, "CCWhiteSpruce2"]*(1-0.50)  # WhiteSpruce CC 20-40yr is 50% of way to 20-40yr natural forest - see document on convergence rates
                results.store$landscape.coef[boot, "CCWhiteSpruce3"] <- results.store$landscape.coef[boot, "WhiteSpruce3"]*0.849 + results.store$landscape.coef[boot, "CCWhiteSpruce3"]*(1-0.849)  # WhiteSpruce CC 40-60yr is 84.9% of way to 40-60yr natural forest
                results.store$landscape.coef[boot, "CCWhiteSpruce4"] <- results.store$landscape.coef[boot, "WhiteSpruce4"]*0.96 + results.store$landscape.coef[boot, "CCWhiteSpruce4"]*(1-0.96)  # WhiteSpruce CC 60-80yr is 96% of way to 60-80yr natural forest
                results.store$landscape.coef[boot, "CCPine2"] <- results.store$landscape.coef[boot, "Pine2"]*0.50 + results.store$landscape.coef[boot, "CCPine2"]*(1-0.50)  # Pine CC 20-40 is 50% of way to 20-40yr natural forest - see document on convergence rates
                results.store$landscape.coef[boot, "CCPine3"] <- results.store$landscape.coef[boot, "Pine3"]*0.849 + results.store$landscape.coef[boot, "CCPine3"]*(1-0.849)  # Pine CC 40-60 is 84.9% of way to 40-60yr natural forest
                results.store$landscape.coef[boot, "CCPine4"] <- results.store$landscape.coef[boot, "Pine4"]*0.96 + results.store$landscape.coef[boot, "CCPine4"]*(1-0.96)  # Pine CC 60-80 is 96% of way to 60-80yr natural forest
                results.store$landscape.coef[boot, "CCDeciduous2"] <- results.store$landscape.coef[boot, "Deciduous2"]*0.705 + results.store$landscape.coef[boot, "CCDeciduous2"]*(1-0.705)  # Decid CC 20-40yr is 70.5% of way to 20-40yr natural forest - see document on convergence rates
                results.store$landscape.coef[boot, "CCDeciduous3"] <- results.store$landscape.coef[boot, "Deciduous3"]*0.912 + results.store$landscape.coef[boot, "CCDeciduous3"]*(1-0.912)  # Decid CC 40-60yr is 91.2% of way to 40-60yr natural forest
                results.store$landscape.coef[boot, "CCDeciduous4"] <- results.store$landscape.coef[boot, "Deciduous4"]*0.97 + results.store$landscape.coef[boot, "CCDeciduous4"]*(1-0.97)  # Decid CC 60-80yr is 97% of way to 60-80yr natural forest
                results.store$landscape.coef[boot, "CCMixedwood2"] <- results.store$landscape.coef[boot, "Mixedwood2"]*0.705 + results.store$landscape.coef[boot, "CCMixedwood2"]*(1-0.705)  # Mixedwood CC 20-40yr is 70.5% of way to 20-40yr natural forest - see document on convergence rates
                results.store$landscape.coef[boot, "CCMixedwood3"] <- results.store$landscape.coef[boot, "Mixedwood3"]*0.912 + results.store$landscape.coef[boot, "CCMixedwood3"]*(1-0.912)  # Mixedwood CC 40-60yr is 91.2% of way to 40-60yr natural forest
                results.store$landscape.coef[boot, "CCMixedwood4"] <- results.store$landscape.coef[boot, "Mixedwood4"]*0.97 + results.store$landscape.coef[boot, "CCMixedwood4"]*(1-0.97)  # Mixedwood CC 60-80yr is 97% of way to 60-80yr natural forest
                
                # And same for SE's
                results.store$landscape.se[boot, "CCWhiteSpruce2"] <- results.store$landscape.se[boot, "WhiteSpruce2"]*0.50 + results.store$landscape.se[boot, "CCWhiteSpruce2"]*(1-0.50)  
                results.store$landscape.se[boot, "CCWhiteSpruce3"] <- results.store$landscape.se[boot, "WhiteSpruce3"]*0.849 + results.store$landscape.se[boot, "CCWhiteSpruce3"]*(1-0.849)
                results.store$landscape.se[boot, "CCWhiteSpruce4"] <- results.store$landscape.se[boot, "WhiteSpruce4"]*0.96 + results.store$landscape.se[boot, "CCWhiteSpruce4"]*(1-0.96) 
                results.store$landscape.se[boot, "CCPine2"] <- results.store$landscape.se[boot, "Pine2"]*0.50 + results.store$landscape.se[boot, "CCPine2"]*(1-0.50)
                results.store$landscape.se[boot, "CCPine3"] <- results.store$landscape.se[boot, "Pine3"]*0.849 + results.store$landscape.se[boot, "CCPine3"]*(1-0.849)
                results.store$landscape.se[boot, "CCPine4"] <- results.store$landscape.se[boot, "Pine4"]*0.96 + results.store$landscape.se[boot, "CCPine4"]*(1-0.96)
                results.store$landscape.se[boot, "CCDeciduous2"] <- results.store$landscape.se[boot, "Deciduous2"]*0.705 + results.store$landscape.se[boot, "CCDeciduous2"]*(1-0.705)
                results.store$landscape.se[boot, "CCDeciduous3"] <- results.store$landscape.se[boot, "Deciduous3"]*0.912 + results.store$landscape.se[boot, "CCDeciduous3"]*(1-0.912)
                results.store$landscape.se[boot, "CCDeciduous4"] <- results.store$landscape.se[boot, "Deciduous4"]*0.97 + results.store$landscape.se[boot, "CCDeciduous4"]*(1-0.97)
                results.store$landscape.se[boot, "CCMixedwood2"] <- results.store$landscape.se[boot, "Mixedwood2"]*0.705 + results.store$landscape.se[boot, "CCMixedwood2"]*(1-0.705)
                results.store$landscape.se[boot, "CCMixedwood3"] <- results.store$landscape.se[boot, "Mixedwood3"]*0.912 + results.store$landscape.se[boot, "CCMixedwood3"]*(1-0.912)
                results.store$landscape.se[boot, "CCMixedwood4"] <- results.store$landscape.se[boot, "Mixedwood4"]*0.97 + results.store$landscape.se[boot, "CCMixedwood4"]*(1-0.97)                
        }
        
        #
        # 3.0 Adjust coefficients that are often estimated poorly.
        #
        
        if (coef.adjust == TRUE) { # If true, adjust coefficients, otherwise leave as is. Can help to determine how much the adjustment is needed
                
                results.store$landscape.coef[boot, "HardLin"] <- plogis( (qlogis(results.store$landscape.coef[boot, "HardLin"]) / results.store$landscape.se[boot, "HardLin"]^2 + qlogis(results.store$landscape.coef[boot, "UrbInd"])/results.store$landscape.se[boot, "UrbInd"]^2) / (1 / results.store$landscape.se[boot, "HardLin"]^2 + 1 / results.store$landscape.se[boot, "UrbInd"]^2) )  # Inverse-variance weighted, done on logit scale then converted back
                results.store$landscape.se[boot, "HardLin"] <- sqrt(1 / (1 / results.store$landscape.se[boot, "HardLin"]^2 + 1 / results.store$landscape.se[boot, "UrbInd"]^2))
                
                if (region == "North" & age.class == TRUE) { 
                        
                        young <- (pSoftLin[1]*results.store$landscape.coef[boot, "CCWhiteSpruceR"] + pSoftLin[2]*results.store$landscape.coef[boot, "CCPineR"] + pSoftLin[3]*results.store$landscape.coef[boot, "CCDeciduousR"] + pSoftLin[4]*results.store$landscape.coef[boot, "BlackSpruce1"]) / sum(pSoftLin)  
                        young.se <- sqrt( ( pSoftLin[1]*results.store$landscape.se[boot, "CCWhiteSpruceR"]^2 + pSoftLin[2]*results.store$landscape.se[boot,"CCPineR"]^2 + pSoftLin[3]*results.store$landscape.se[boot, "CCDeciduousR"]^2 + pSoftLin[4]*results.store$landscape.se[boot, "BlackSpruce1"]^2 ) / sum(pSoftLin) )  # Weighted se (on logit scale) - does not include component of variance among types, because fixed weighting
                        
                        # Adjust for the three soft linear type coefficients (EnSoftLin, EnSeismic, TrSoftLin)
                        results.store$landscape.coef[boot, "EnSoftLin"] <- plogis( (qlogis(results.store$landscape.coef[boot, "EnSoftLin"])/results.store$landscape.se[boot, "EnSoftLin"]^2 + qlogis(young) / young.se^2) / (1 / results.store$landscape.se[boot, "EnSoftLin"]^2 + 1 / young.se^2) )  # Inverse-variance weighted, done on logit scale then converted back
                        results.store$landscape.se[boot, "EnSoftLin"] <- sqrt(1 / (1 / results.store$landscape.se[boot, "EnSoftLin"]^2 + 1 / young.se^2))
                        
                        results.store$landscape.coef[boot, "EnSeismic"] <- plogis( (qlogis(results.store$landscape.coef[boot, "EnSeismic"])/results.store$landscape.se[boot, "EnSeismic"]^2 + qlogis(young) / young.se^2) / (1 / results.store$landscape.se[boot, "EnSeismic"]^2 + 1 / young.se^2) )  # Inverse-variance weighted, done on logit scale then converted back
                        results.store$landscape.se[boot, "EnSeismic"] <- sqrt(1 / (1 / results.store$landscape.se[boot, "EnSeismic"]^2 + 1 / young.se^2))
                        
                        results.store$landscape.coef[boot, "TrSoftLin"] <- plogis( (qlogis(results.store$landscape.coef[boot, "TrSoftLin"])/results.store$landscape.se[boot, "TrSoftLin"]^2 + qlogis(young) / young.se^2) / (1 / results.store$landscape.se[boot, "TrSoftLin"]^2 + 1 / young.se^2) )  # Inverse-variance weighted, done on logit scale then converted back
                        results.store$landscape.se[boot, "TrSoftLin"] <- sqrt(1 / (1 / results.store$landscape.se[boot, "TrSoftLin"]^2 + 1 / young.se^2))
                        
                        rm(young, young.se)
                        
                }
                
        } 
        
        #
        # 4.0 Residual variation due to space and climate
        #
        
        # Landcover
        p1 <- colSums(results.store$landscape.coef[boot, ] * t(data.in[ ,colnames(results.store$landscape.coef)]))  # Multiply coefficients by landscape type proportion for each site to get offset
        
        # Add the detection information
        p1 <- plogis(qlogis(p1) + colSums(results.store$detection.coef[boot, ] * t(data.in[ ,colnames(results.store$detection.coef)]))) 
        
        # auc(ifelse(data.in[, "Count"] > 0, 1, 0), p1)
        
        if (region == "South") {
                
                p1 <- plogis(qlogis(p1) + results.store$paspen.coef[boot, ] * data.in$paspen)  # And add in effect of pAspen on logit scale
                
        }
        
        # 4.1 Climate
        space.climate.store[[1]] <- try(glm(pcount ~ offset(qlogis(0.998 * p1 + 0.001)), data = data.in, family = "binomial"))  # Prediction can be 0 (water), so need to scale a bit to avoid infinities in logit offset. Protocol not included here, because already accounted for in offset values 
        space.climate.store[[2]] <- try(update(space.climate.store[[1]], .~.+ PET))
        space.climate.store[[3]] <- try(update(space.climate.store[[1]], .~.+ AHM))
        space.climate.store[[4]] <- try(update(space.climate.store[[1]], .~.+ MAT))
        space.climate.store[[5]] <- try(update(space.climate.store[[1]], .~.+ FFP))
        space.climate.store[[6]] <- try(update(space.climate.store[[1]], .~.+ MAP + FFP))
        space.climate.store[[7]] <- try(update(space.climate.store[[1]], .~.+ MAP +FFP + MAPFFP))
        space.climate.store[[8]] <- try(update(space.climate.store[[1]], .~.+ MAT + MAP + PET + AHM))
        space.climate.store[[9]] <- try(update(space.climate.store[[1]], .~.+ MAT + MAP + PET + AHM + MAPPET + MATAHM))
        space.climate.store[[10]] <- try(update(space.climate.store[[1]] ,.~.+ MAT + MAP))
        space.climate.store[[11]] <- try(update(space.climate.store[[1]] ,.~.+ MWMT + MCMT))
        space.climate.store[[12]] <- try(update(space.climate.store[[1]] ,.~.+ AHM + PET))
        space.climate.store[[13]] <- try(update(space.climate.store[[1]] ,.~.+ MAT + MAT2 + MWMT + MWMT2))
        space.climate.store[[14]] <- try(update(space.climate.store[[1]] ,.~.+ MWMT + MCMT + FFP + MAT))
        
        # 4.2 Spatial + Water
        for (i in 1:14) space.climate.store[[i + 14]] <- try(update(space.climate.store[[i]], .~.+ Lat + Long + LatLong)) # Adjusted Lat in the northern model so southern sites are treated as further north.
        for (i in 1:14) space.climate.store[[i + 28]] <- try(update(space.climate.store[[i]], .~.+ Lat + Long + LatLong + Lat2 + Long2))
        for (i in 1:42) space.climate.store[[i + 42]] <- try(update(space.climate.store[[i]], .~.+ Waterkm + Waterkm2)) # Including the amount of surrounding water within the buffer

        nModels.sc<-length(space.climate.store)
        # BIC calculation to select best covariate set Uses BIC for more conservative variable set
        bic.sc <- rep(999999999, nModels.sc)
        for (i in 1:(nModels.sc)) {
                
                if (!is.null(space.climate.store[[i]]) & class(space.climate.store[[i]])[1] != "try-error") {
                        bic.sc[i] <- BIC(space.climate.store[[i]])
                }
                
        }
        
        bic.delta.sc <- bic.sc - min(bic.sc)
        bic.exp.sc <- exp(-1 / 2 * bic.delta.sc)
        bic.wt.sc <- bic.exp.sc / sum(bic.exp.sc)
        best.model.sc <- space.climate.store[[which.max(bic.wt.sc)]]
        
        # And the subset of climate and/or spatial variables
        vnames <- names(coef(space.climate.store[[which.max(bic.wt.sc)]]))
        vnames1 <- ifelse(vnames == "(Intercept)", "Intercept", vnames)  # Correct the intercept name
        
        # Storing of climate coefficients
        coef.match <- match(colnames(results.store$climate.coef), vnames1, nomatch = 0)
        results.store$climate.coef[boot, coef.match != 0] <- coef(best.model.sc)[coef.match]
        
        # 4.2 Store detections, surveys, and AIC information
        results.store$fit[boot, "dect"] <- sum(data.in$Count) 
        results.store$fit[boot, "survey"] <- sum(data.in$nQuadrant)
        
        # Adjusting the hard linear coefficient (default to 0 as amphibians are not breeding on hard linear features)
        # Soft Linear coefficients were divided into three categories, but are no longer adjusted.
        
        temp.coef <- results.store$landscape.coef[boot, ]
        temp.coef["HardLin"] <- 0
        
        ## No longer perform this correction to soft linear as of 2020-11-17
        # if (region == "North") {
        # 
        #         temp.coef["SoftLin"] <- mean(temp.coef[c("Grass", "Shrub")])
        # 
        # } else {
        # 
        #         temp.coef["SoftLin"] <- mean(temp.coef[c("Productive", "Clay", "Saline", "RapidDrain")])
        # 
        # }
        
        # Vegetation prediction
        
        km2.pveg.curr <- colSums(temp.coef*t(data.in[ ,names(temp.coef)])) # Prediction based on veg types only.
        
        # Calculate cimate and spatial prediction for each species
        
        km2.pres.curr <- colSums(results.store$climate.coef[boot, c(FALSE, as.logical(!is.na(results.store$climate.coef[boot, -1])))] * t(data.in[ ,colnames(results.store$climate.coef)[c(FALSE, as.logical(!is.na(results.store$climate.coef[boot, -1])))]])) # Prediction of residual (climate and spatial) effect
        
        # Take the difference between unqiue species intercept and the predicted climate residuals
        
        km2.pres.curr <- km2.pres.curr + results.store$climate.coef[boot, 1]
        
        # If the region is the South, add the pAspen prediction
        if (region == "South") {
                
                km2.pAspen <- colSums(results.store$paspen.coef[boot, "paspen"] * t(data.in[, "paspen"])) 
                km2.p.curr <- plogis(qlogis(0.998*km2.pveg.curr+0.001) + km2.pres.curr + km2.pAspen) # Apply same transformation as used in fitting residual model
                
        } else {
                
                km2.p.curr <- plogis(qlogis(0.998*km2.pveg.curr+0.001) + km2.pres.curr) # Apply same transformation as used in fitting residual model
        }
        
        # Add detection prediction
        km2.pveg.curr2 <- plogis(qlogis(km2.pveg.curr) + colSums(results.store$detection.coef[boot, ] * t(data.in[ ,colnames(results.store$detection.coef)])))
        km2.p.curr2 <- plogis(qlogis(km2.p.curr) + colSums(results.store$detection.coef[boot, ] * t(data.in[ ,colnames(results.store$detection.coef)])))
        
        # AIC calcualtion
        # results.store$fit[boot, "auc_LC"] <- auc(ifelse(data.in[, "Count"] > 0, 1, 0), km2.pveg.curr) # AUC Calculation with landcover only
        # results.store$fit[boot, "auc_both"] <- auc(ifelse(data.in[, "Count"] > 0, 1, 0), km2.p.curr) # AUC Calculation with both
        # 
        results.store$fit[boot, "auc_LC"] <- auc(ifelse(data.in[, "Count"] > 0, 1, 0), km2.pveg.curr2) # AUC Calculation with landcover only
        results.store$fit[boot, "auc_both"] <- auc(ifelse(data.in[, "Count"] > 0, 1, 0), km2.p.curr2) # AUC Calculation with both
        
        # Return results
        return(results.store)
        
}

