######################################################################################
## Build 3 .rdata with
## 1. Ignitions
## 2. Coordinates of the study area
## 3. Orographic variables
######################################################################################

read.static.vars <- function(){
  
  library(raster)
  library(tidyverse)
  
  cat("Read coordinates, orography and ignitions", "\n")
  
  ## Build a coordinates data frame from MASK raster
  load("inputlyrs/rdata/mask.rdata")
  coord <- data.frame(cell.id=1:ncell(MASK), coordinates(MASK), mask=MASK[])
  coord <- filter(coord, !is.na(mask)) %>% select(-mask)
  save(coord, file="inputlyrs/rdata/coordinates.rdata") 
  
  ## Read orography variables, build and save the data frame
  ELEVATION <- raster("inputlyrs/asc/DEM_100m_31N-ETRS89.asc")
  ASPECT <- raster("inputlyrs/asc/Aspect_100m_31N-ETRS89.asc")
  orography <- data.frame(cell.id=1:ncell(MASK), elev=ELEVATION[], aspect=ASPECT[])
  orography <- orography[!is.na(MASK[]),]
  orography$elev[is.na(orography$elev)] <- rnorm(sum(is.na(orography$elev)),
                                                 mean(orography$elev, na.rm=T), sd(orography$elev, na.rm=T))
  orography$aspect[is.na(orography$aspect)] <- sample(1:4, sum(is.na(orography$aspect)), replace=T)
  save(orography, file="inputlyrs/rdata/orography.rdata")
  
  ## Read ignitions, build and save the data frame
  IGNIS <- raster("inputlyrs/asc/FireIgnis8912_31N-ETRS89.asc")
  ignis <- data.frame(cell.id=1:ncell(MASK), fire.id=IGNIS[])
  ignis <- ignis[!is.na(MASK[]),]
  save(ignis, file="inputlyrs/rdata/ignitions.rdata")
  
}
