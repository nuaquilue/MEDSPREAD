rm(list = ls())
library(raster)
library(tidyverse)

## Load ignitions
load("inputlyrs/rdata/ignitions.rdata")

## Fire perimetres 1989 - 1999
load("inputlyrs/rdata/mask.rdata")
perim.id <- data.frame(cell.id=1:ncell(MASK[]), mask=MASK[])
for(y in 1989:1999){
  cat(paste("Fire perimeters year", y), "\n")
  ## Read fire perimeter, project to 31N-ETRS89 and crop to the default extent
  FIREIDed50 <- raster(paste0("c:/work/MEDMOD/InputLayers_MEDSPREAD/FirePerimetersAscii/firesID", y))
  crs(FIREIDed50) <- CRS("+init=epsg:23031")
  FIREIDetrs <- projectRaster(FIREIDed50, res=100, crs=CRS("+init=epsg:25831"), method="ngb")
  FIREID <- crop(FIREIDetrs, FIREIDed50, snap="near")
  ## Filter ignis of the current year
  ids <- unique(FIREID[])[-(1:2)]
  ignis.y <- filter(ignis, id %in% ids)
  ## Assign the perim.id of perimeters to a new variable X, and make sure the ignition is 
  ## in the fire perimeter
  perim.id$x <- FIREID[]
  perim.id <- left_join(perim.id, ignis.y, by="cell.id")
  perim.id$x[!is.na(perim.id$id)] <- perim.id$id[!is.na(perim.id$id)]
  perim.id <- select(perim.id, -id)
  names(perim.id)[ncol(perim.id)] <- paste0("y", y)
}
perim.id <- filter(perim.id, !is.na(mask)) %>% select(-mask)
save(perim.id, file="inputlyrs/rdata/fireperim.89-99.rdata")

# New MASK
LCF <- raster("inputlyrs/asc/ForestMapSpp00_31N-ETRS89.asc")
MASK <- LCF
MASK[!is.na(MASK[])] <- 1
crs(MASK) <- CRS("+init=epsg:25831")
perim.id <- data.frame(cell.id=1:ncell(MASK[]), mask=MASK[])
for(y in 2000:2012){
  cat(paste("Fire perimeters year", y), "\n")
  ## Read fire perimeter, project to 31N-ETRS89 and crop to the default extent
  FIREIDed50 <- raster(paste0("c:/work/MEDMOD/InputLayers_MEDSPREAD/FirePerimetersAscii/firesID", y))
  crs(FIREIDed50) <- CRS("+init=epsg:23031")
  FIREIDetrs <- projectRaster(FIREIDed50, res=100, crs=CRS("+init=epsg:25831"), method="ngb")
  FIREID <- crop(FIREIDetrs, FIREIDed50, snap="near")
  ## Filter ignis of the current year
  ids <- unique(FIREID[])[-(1:2)]
  ignis.y <- filter(ignis, id %in% ids)
  ## Assign the perim.id of perimeters to a new variable X, and make sure the ignition is 
  ## in the fire perimeter
  perim.id$x <- FIREID[]
  perim.id <- left_join(perim.id, ignis.y, by="cell.id")
  perim.id$x[!is.na(perim.id$id)] <- perim.id$id[!is.na(perim.id$id)]
  perim.id <- select(perim.id, -id)
  names(perim.id)[ncol(perim.id)] <- paste0("y", y)
}
perim.id <- filter(perim.id, !is.na(mask)) %>% select(-mask)
save(perim.id, file="inputlyrs/rdata/fireperim.00-12.rdata")
