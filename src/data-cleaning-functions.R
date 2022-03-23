#
# Title: Data cleaning functions
# Created: June 13th, 2019
# Last Updated: January 24th, 2022
# Author: Brandon Allen
# Objectives: Functions required to summarize the landscape long-form data
# Keywords: Landscape summaries, Kgrid Initialization
#

#######################
# Landscape summaries # 
#######################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

landscape_hf_summary <- function(data.in, landscape.lookup, class.in, class.out) {
        
        # Matching of lookup tables and merging native features
        landscape.lookup <- landscape.lookup[landscape.lookup[, class.in] %in% colnames(data.in), ]
        landscape.clean <- matrix(nrow = nrow(data.in), ncol = length(unique(landscape.lookup[, class.out])))
        
        for (abmi.coef in 1:length(unique(landscape.lookup[, class.out]))) {
                
                coef.temp <- as.character(landscape.lookup[landscape.lookup[, class.out] %in% as.character(unique(landscape.lookup[, class.out]))[abmi.coef], class.in])
                
                if(length(coef.temp) == 1) {
                        
                        landscape.clean[, abmi.coef] <- data.in[, coef.temp]
                        
                } else {
                        
                        landscape.clean[, abmi.coef] <- rowSums(data.in[, coef.temp])
                        
                }
                
        }
        
        colnames(landscape.clean) <- as.character(unique(landscape.lookup[, class.out]))
        rownames(landscape.clean) <- rownames(data.in)
        
        return(landscape.clean)
        
}

########################
# Kgrid Initialization #
########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create_kgrid <- function(kgrid.data) {
        
        load(kgrid.data)
        names(kgrid)[which(names(kgrid)=="POINT_X")]<-"Long"
        names(kgrid)[which(names(kgrid)=="POINT_Y")]<-"Lat"
        return(data.frame(LinkID=rownames(kgrid),
                          Lat=kgrid$Lat,
                          Long=kgrid$Long,
                          NR=kgrid$NRNAME,
                          NSR=kgrid$NSRNAME,
                          LUF=kgrid$LUF_NAME,
                          Col=kgrid$Col,
                          Row=kgrid$Row))
        gc()
        
}