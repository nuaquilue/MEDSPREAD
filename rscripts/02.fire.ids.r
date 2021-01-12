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
  writeRaster(FIREID, paste0("c:/work/MEDMOD/InputLayers_MEDSPREAD/Fires8912_31N-ETRS89/FiresID", y),
              format="GTiff", overwrite=T)
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
  writeRaster(FIREID, paste0("c:/work/MEDMOD/InputLayers_MEDSPREAD/Fires8912_31N-ETRS89/FiresID", y),
              format="GTiff", overwrite=T)
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


###### Build 3 layers, one per each fire spread type
## Function to select items not in a vector
`%notin%` <- Negate(`%in%`)
##
for(fst in c("Wind", "Conv", "Topo")){
  fire.ignis <- read.table(paste0("inputfiles/FireIgnitions", fst, ".txt"), header=T)
  load("inputlyrs/rdata/fireperim.89-99.rdata")
  for(y in 1:11){
    ids <- unlist(filter(fire.ignis, year==y) %>% select(fire.id))
    aux <- perim.id[,y+1] 
    aux[aux %notin% ids] <- 0
    # aux[aux %in% ids] <- paste0(y+1988, "_", aux[aux %in% ids])
    perim.id[,y+1] <- aux
  }
  perim.id_8999 <- perim.id
  load("inputlyrs/rdata/fireperim.00-12.rdata")
  for(y in 1:13){
    ids <- unlist(filter(fire.ignis, year==y+11) %>% select(fire.id))
    aux <- perim.id[,y+1] 
    aux[aux %notin% ids] <- 0
    # aux[aux %in% ids] <- paste0(y+1999, "_", aux[aux %in% ids])
    perim.id[,y+1] <- aux
  }
  ## All perims (from 1989 to 2012)
  perim.id.all <- left_join(perim.id, perim.id_8999, by="cell.id")
  rm(perim.id_8999); rm(perim.id); rm(aux); gc()
  fire.id <- data.frame(cell.id=perim.id.all$cell.id, id=apply(perim.id.all[,-1], 1, max))
  ## Translate the info into a Raster
  load("inputlyrs/rdata/mask.00-12.rdata")
  MASK[!is.na(MASK[])] <- fire.id$id
  writeRaster(MASK, paste0("c:/work/MEDMOD/InputLayers_MEDSPREAD/Fires8912_31N-ETRS89/Fires", fst),
              format="GTiff", overwrite=T)
}

