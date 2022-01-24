#
# Title: Functions to visualize coefficients and predictions
# Created: August 8th, 2018
# Last Updated: January 24th, 2022
# Author: Brandon Allen
# Objective: Function for recreating coefficient plots
# Keywords: Coefficients
# 

################
# Coefficients #
################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coefficient_plots <- function (species, spp.coef, spp.boot, estimate.id, region) {
    
    require(ggplot2)
    require(abmi.themes)
    
    #
    # 0.1 Subset species coefficients and identify the model with the complete data set (row 1)
    #
    
    spp.coef <- spp.coef[[species]]
    
    if (region == "South") {
        
        #
        # 1. Soil models
        #
        
        # 1.1 Organize coefficients
        main.coef <- as.data.frame(as.matrix(spp.coef$landscape.coef[1, c("Loamy", "SandyLoam", "RapidDrain", "ClaySub", "ThinBreak", "Blowout", "Other",
                                                                          "Crop", "TameP", "RoughP", "Wellsites", "Rural", "Urban", "Industrial")]))
        main.coef <- cbind.data.frame(rownames(main.coef), main.coef[,1])
        colnames(main.coef) <- c("Name", "Coef")
        main.coef$Name <- factor(c("Loamy", "Sandy/loamy", "Rapid Drain", "Clay", "Thin Break", "Blowout", "Other soil types",
                                   "Cropland", "Tame pasture", "Rough pasture", "Well sites", "Rural residential", "Urban/Industrial", "Industrial (rural)"), levels = c("Loamy", "Sandy/loamy", "Rapid Drain", "Clay", "Thin Break", "Blowout", "Other soil types",
                                                                                                                                                                         "Cropland", "Tame pasture", "Rough pasture", "Well sites", "Rural residential", "Urban/Industrial", "Industrial (rural)")) # Rename coefficients and update as factors
        
        # 1.2 Estimating 90% confidence intervals
        main.boot <- spp.boot[as.character(main.coef$Name), ]
        main.coef["LowerCI"] <- main.boot[, paste(species, "Min")]
        main.coef["UpperCI"] <- main.boot[, paste(species, "Max")]
        
        # 1.3 Estimating coefficients modified by paspen
        main.coef["Mod_Coef"] <- plogis(qlogis(main.coef$Coef) + spp.coef$paspen.coef[1, ])
        main.coef["Mod_LowerCI"] <- main.coef$Mod_Coef*main.coef$LowerCI/main.coef$Coef
        main.coef["Mod_UpperCI"] <- main.coef$Mod_Coef*main.coef$UpperCI/main.coef$Coef
        main.coef[is.na(main.coef)] <- 0
        main.coef[, 2:ncol(main.coef)][main.coef[, 2:ncol(main.coef)] > 1] <- 1
        
        colour.pal <- c("#DAD157", "#663301", "#60AE9F", "#663301",
                        "#63A70C", "#748838", "#FE9929",
                        "#C89222", "#60A52A", "#DCCF63", "#9D350B", 
                        "#532E8C", "#448278", "#1D3991")
        
        # 1.4 Non-treed
        png(file = paste("results/figures/coefficients/south/", species, "-non-treed_", Sys.Date(), ".png", sep = ""),
            width = 1500,
            height = 1500, 
            res = 300)
        
        non.treed <- ggplot(data = main.coef, aes(x = Name, y = Coef, fill = Name)) +
            geom_bar(stat="identity", fill = colour.pal) +
            geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width=.2,
                          position=position_dodge(.9)) + 
            guides(scale = "none") + 
            ylim(c(0,ifelse(max(1.1*main.coef$Mod_UpperCI) >= 1, 1, max(1.1*main.coef$Mod_UpperCI)))) +
            labs(x = "Non-Treed Coefficient", y = "Relative Abundance") +
            geom_text(aes(x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
                          label = paste("Detections:", spp.coef$fit[1, "dect"], "; Surveys:", spp.coef$fit[1, "survey"], "; AUC:", round(spp.coef$fit[1, "auc_both"], 3))), 
                      size = 10) +
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
       
        print(non.treed)
        
        dev.off()
        
        # 1.5 Treed
        png(file = paste("results/figures/coefficients/south/", species, "-treed_", Sys.Date(), ".png", sep = ""),
            width = 1500,
            height = 1500, 
            res = 300)
        
        treed <- ggplot(data = main.coef, aes(x = Name, y = Mod_Coef, fill = Name)) +
            geom_bar(stat="identity", fill = colour.pal) +
            geom_errorbar(aes(ymin = Mod_LowerCI, ymax = Mod_UpperCI), width=.2,
                          position=position_dodge(.9)) + 
            guides(scale = "none") + 
                ylim(c(0,ifelse(max(1.1*main.coef$Mod_UpperCI) >= 1, 1, max(1.1*main.coef$Mod_UpperCI)))) +
            labs(x = "Treed Coefficient", y = "Relative Abundance") +
            geom_text(aes(x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
                          label = paste("Detections:", spp.coef$fit[1, "dect"], "; Surveys:", spp.coef$fit[1, "survey"], "; AUC:", round(spp.coef$fit[1, "auc_both"], 3))), 
                      size = 10) +
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
        
        print(treed)
        dev.off()
        
        # 1.6 Combined
        
        png(file = paste("results/figures/coefficients/south/", species, "-combined_", Sys.Date(), ".png", sep = ""),
            width = 3000,
            height = 1500, 
            res = 300)
        
        combined <- ggarrange(non.treed, treed,
                            ncol = 2, nrow = 1)
        
        print(combined)
        
        dev.off()
        
        # 1.7 Linear features
        
        linear.data <- data.frame(Name = c("Soft Linear", "Hard Linear", "Soft Linear", "Hard Linear"), 
                                  Coef = NA, Type = factor(c("None", "None", "10% Linear", "10% Linear"), levels = c("None", "10% Linear")))
        
        linear.data$Coef[linear.data$Type == "None"] <- mean(spp.coef$landscape.coef[1, c("Loamy", "SandyLoam", "RapidDrain", "ClaySub", "ThinBreak", "Blowout", "Other")])
        linear.data[3, "Coef"] <- linear.data[1, "Coef"] * 0.9 + mean(spp.coef$landscape.coef[1, c("EnSoftLin", "TrSoftLin")]) * 0.1
        linear.data[4, "Coef"] <- linear.data[1, "Coef"] * 0.9 + spp.coef$landscape.coef[1, "HardLin"] * 0.1
        
        y.limit <- max(linear.data$Coef,2*linear.data[1, "Coef"])*1.03
        
        png(file = paste("results/figures/coefficients/south/", species, "_linear-features_", Sys.Date(), ".png", sep = ""),
            width = 1800,
            height = 1800, 
            res = 300)
        
        print(ggplot(data = linear.data, aes(x = Type, y = Coef)) +
                  geom_line(aes(group = Name)) +
                  geom_point(aes(col = Name, shape = Name), size = 3) +
                  labs(x = "Human Footprint", y = "Relative Abundance") + 
                  scale_shape_manual(values=c(15, 19), name = "Coefficients") +
                  scale_color_manual(values=c(abmi_pal("main")(2)), name = "Coefficients") +
                  ylim(c(0, y.limit)) +
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
    
    if (region == "North") {
        
        # 2.1 Organize coefficients
        main.coef <- as.data.frame(as.matrix(spp.coef$landscape.coef[1, c("WhiteSpruceR", "WhiteSpruce1", "WhiteSpruce2", "WhiteSpruce3", "WhiteSpruce4", "WhiteSpruce5", "WhiteSpruce6", "WhiteSpruce7", "WhiteSpruce8", 
                                                                          "PineR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6", "Pine7", "Pine8",
                                                                          "DeciduousR", "Deciduous1", "Deciduous2", "Deciduous3", "Deciduous4", "Deciduous5", "Deciduous6", "Deciduous7", "Deciduous8",
                                                                          "MixedwoodR", "Mixedwood1", "Mixedwood2", "Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7", "Mixedwood8",
                                                                          "TreedBogR", "TreedBog1", "TreedBog2", "TreedBog3", "TreedBog4", "TreedBog5", "TreedBog6", "TreedBog7", "TreedBog8",
                                                                          "TreedFenR", "TreedSwamp", "ShrubbySwamp", "ShrubbyBog", "ShrubbyFen", "GraminoidFen", "Marsh", "Shrub", "GrassHerb", 
                                                                          "Crop", "TameP", "RoughP", "Wellsites", "Rural", "Urban", "Industrial")]))
        main.coef <- cbind.data.frame(rownames(main.coef), main.coef[,1])
        colnames(main.coef) <- c("Name", "Coef")
        main.coef$Name <- factor(main.coef$Name, levels = main.coef$Name) # Define factor order for proper plotting
        
        cc.coef <- as.data.frame(as.matrix(spp.coef$landscape.coef[1, c("CCWhiteSpruceR", "CCWhiteSpruce1", "CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4",
                                                                        "CCPineR", "CCPine1", "CCPine2", "CCPine3", "CCPine4",
                                                                        "CCDeciduousR", "CCDeciduous1", "CCDeciduous2", "CCDeciduous3", "CCDeciduous4",
                                                                        "CCMixedwoodR", "CCMixedwood1", "CCMixedwood2", "CCMixedwood3", "CCMixedwood4")])) # Clear cuts
        cc.coef <- cbind.data.frame(rownames(cc.coef), cc.coef[,1])
        colnames(cc.coef) <- c("Name", "Coef")
        cc.coef$Name <- factor(cc.coef$Name, levels = cc.coef$Name) # Define factor order for proper plotting
        
        # 2.2 Estimating 90% confidence intervals
        main.boot <- spp.boot[as.character(main.coef$Name), ]
        cc.boot <- spp.boot[c(as.character(cc.coef$Name)), ]
        
        main.coef["LowerCI"] <- main.boot[, paste(species, "Min")]
        main.coef["UpperCI"] <- main.boot[, paste(species, "Max")]
        
        cc.coef["LowerCI"] <- cc.boot[, paste(species, "Min")] # Clear cuts
        cc.coef["UpperCI"] <- cc.boot[, paste(species, "Max")]
        cc.coef.low <- c(cc.coef$LowerCI[1:5], rep(NA, 4), cc.coef$LowerCI[6:10], rep(NA, 4), cc.coef$LowerCI[11:15], rep(NA, 4), cc.coef$LowerCI[16:20], rep(NA, 29))
        cc.coef.up <- c(cc.coef$UpperCI[1:5], rep(NA, 4), cc.coef$UpperCI[6:10], rep(NA, 4), cc.coef$UpperCI[11:15], rep(NA, 4), cc.coef$UpperCI[16:20], rep(NA, 29))
        
        main.coef[is.na(main.coef)] <- 0
        main.coef[, 2:4][main.coef[, 2:4] > 1] <- 1
        cc.coef[is.na(cc.coef)] <- 0
        cc.coef[, 2:4][cc.coef[, 2:4] > 1] <- 1
        
        # Rename the old vegetation coefficients
        main.coef$Name <- c("White Spruce 0-9", "White Spruce 10-19", "White Spruce 20-39", "White Spruce 40-59", "White Spruce 60-79", "White Spruce 80-99", "White Spruce 100-119", "White Spruce 120-139", "White Spruce 140+", 
                              "Pine 0-9", "Pine 10-19", "Pine 20-39", "Pine 40-59", "Pine 60-79", "Pine 80-99", "Pine 100-119", "Pine 120-139", "Pine 140+",
                              "Deciduous 0-9", "Deciduous 10-19", "Deciduous 20-39", "Deciduous 40-59", "Deciduous 60-79", "Deciduous 80-99", "Deciduous 100-119", "Deciduous 120-139", "Deciduous 140+",
                              "Mixedwood 0-9", "Mixedwood 10-19", "Mixedwood 20-39", "Mixedwood 40-59", "Mixedwood 60-79", "Mixedwood 80-99", "Mixedwood 100-119", "Mixedwood 120-139", "Mixedwood 140+",
                              "Treed Bog 0-9", "Treed Bog 10-19", "Treed Bog 20-39", "Treed Bog 40-59", "Treed Bog 60-79", "Treed Bog 80-99", "Treed Bog 100-119", "Treed Bog 120-139", "Treed Bog 140+",
                              "Treed Fen", "Treed Swamp", "Shrubby Swamp", "Shrubby Bog", "Shrubby Fen", "Graminoid Fen", "Marsh", "Shrub", "Grass", 
                              "Crop", "Tame Pasture", "Rough Pasture", "Wellsites", "Rural", "Urban/Industrial", "Industrial (Rural)")
        main.coef$Name <- factor(main.coef$Name, levels = main.coef$Name) # Define factor order for proper plotting
        
        # 2.3 Vegetation 
        
        # ABMI Biobrowser colour palette
        colour.pal <- c(rep("#9A9723", 9), rep("#74863F", 9), rep("#316413", 9), rep("#649869", 9), rep("#676514", 9), 
                        "#609ACA", "#7DBA4E", "#CF9823", "#4A902F", "#2B6797",
                        "#F17122", "#2E3996", "#38A963", "#306537",
                        "#C89222", "#60A52A", "#DCCF63", "#9D350B", 
                        "#532E8C", "#448278", "#1D3991")
        
        png(file = paste("results/figures/coefficients/north/", species, "_", Sys.Date(), ".png", sep = ""),
            width = 5400,
            height = 1800, 
            res = 300)
        
        print(ggplot(data = main.coef, aes(x = Name, y = Coef, fill = Name)) +
                  geom_bar(stat = "identity", fill = colour.pal) +
                  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width=.2,
                                position = position_dodge(.9)) + 
                  geom_errorbar(aes(ymin = cc.coef.low, ymax = cc.coef.up), width=.2,
                                position = position_dodge(.8), color = "#A8AF8C")  +
                  guides(scale = "none") + 
                  labs(x = "Coefficient", y = "Relative Abundance") +
                  geom_text(aes(x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
                                label = paste("Detections:", spp.coef$fit[1, "dect"], "; Surveys:", spp.coef$fit[1, "survey"], "; AUC:", round(spp.coef$fit[1, "auc_both"], 3))), 
                            size = 14) +
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
        
        # 2.4 Linear features
        
        linear.data <- data.frame(Name = c("Soft Linear", "Hard Linear", "Soft Linear", "Hard Linear"), 
                                  Coef = NA, Type = factor(c("None", "None", "10% Linear", "10% Linear"), levels = c("None", "10% Linear")))

        linear.data$Coef[linear.data$Type == "None"] <- mean(spp.coef$landscape.coef[1, -c(51:72)])
        linear.data[3, "Coef"] <- linear.data[1, "Coef"] * 0.9 + mean(spp.coef$landscape.coef[1, c("TrSoftLin", "EnSoftLin", "EnSeismic")]) * 0.1
        linear.data[4, "Coef"] <- linear.data[1, "Coef"] * 0.9 + spp.coef$landscape.coef[1, c("HardLin")] * 0.1
        
        y.limit <- max(linear.data$Coef,2*linear.data[1, "Coef"])*1.03
        
        png(file = paste("results/figures/coefficients/north/", species, "_linear-features_", Sys.Date(), ".png", sep = ""),
            width = 1800,
            height = 1800, 
            res = 300)
        
        print(ggplot(data = linear.data, aes(x = Type, y = Coef)) +
                  geom_line(aes(group = Name)) +
                  geom_point(aes(col = Name, shape = Name), size = 3) +
                  labs(x = "Human Footprint", y = "Relative Abundance") + 
                  scale_shape_manual(values=c(15, 19), name = "Coefficients") +
                  scale_color_manual(values=c(abmi_pal("main")(2)), name = "Coefficients") +
                  ylim(c(0, y.limit)) +
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
    
}
