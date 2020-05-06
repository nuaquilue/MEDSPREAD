######################################################################################
## Build 5 .rdata with
## 1. Raster MASK of the study area 
## 2. Coordinates of the study area
## 3. Orographic variables
## 4. UTM grid as a data frame
## 5. Fire regime related variables
######################################################################################

read.static.vars <- function(work.path){
  
  library(raster)
  library(RANN)
  library(tidyverse)
  
  cat("Reading orographyic, utm, fire regime variables, harvesting-restrictions", "\n")
  
  ## MASK of the study area
  MASK <- raster(paste0(work.path, "/inputlyrs/asc/LCFspp_100m_31N-ETRS89.asc"))
  MASK[!is.na(MASK[])] <- 1
  crs(MASK) <- CRS("+init=epsg:25831")
  save(MASK, file="inputlyrs/rdata/mask.rdata") 
  
  ## Build a coordinates data frame and then load model static variables in a data.frame format
  coord <- data.frame(cell.id=1:ncell(MASK), coordinates(MASK), mask=MASK[])
  coord <- filter(coord, !is.na(mask)) %>% select(-mask)
  save(coord, file="inputlyrs/rdata/coordinates.rdata") 
  
  ## Read initial state vars,  build and save the data frame
  ELEVATION <- raster(paste0(work.path, "/inputlyrs/asc/DEM_100m_31N-ETRS89.asc"))
  ASPECT <- raster(paste0(work.path, "/inputlyrs/asc/Aspect_100m_31N-ETRS89.asc"))
  SLOPE <- raster(paste0(work.path, "/inputlyrs/asc/SlopeDegree_100m_31N-ETRS89.asc"))
  ROAD <- raster(paste0(work.path, "/inputlyrs/asc/DensRoad_100m_31N-ETRS89.asc"))
  ROAD[is.na(ROAD[])] <- 0  ## why there are so many NA is ROAD layer?
  orography <- data.frame(cell.id=1:ncell(MASK), elev=ELEVATION[], aspect=ASPECT[], slope=SLOPE[], road=ROAD[])
  orography <- orography[!is.na(MASK[]),]
  save(orography, file="inputlyrs/rdata/orography.rdata")
  
  ## UTM layer
  UTM <- raster(paste0(work.path, "/inputlyrs/asc/UTM1k_100m_31N-ETRS89.asc"))
  ## Are there NAs in the UTM layer within CAT?
  dta <- data.frame(cell.id=1:ncell(UTM), m=MASK[], coordinates(UTM), z=UTM[]) %>% filter(!is.na(m))
  na.var <- filter(dta, is.na(z))
  for(id in na.var$cell.id){
    neighs <- nn2(select(dta, x, y), filter(na.var, cell.id==id)%>% select(x,y), 
                  searchtype="priority", k=9)
    phago <- mean(dta$z[neighs$nn.idx], na.rm=T)
    if(!is.na(phago))
      na.var$z[na.var$cell.id==id] <- phago
  }
  na.var2 <- filter(na.var, is.na(z))
  dta$z[is.na(dta$z)] <- na.var$z
  UTM[!is.na(MASK[])] <- dta$z
  utm <- data.frame(cell.id=1:ncell(UTM),  utm=UTM[])
  save(utm, file="inputlyrs/rdata/utm.rdata")
}


