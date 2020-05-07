######################################################################################
## Build land.rdata with the initialization of the 6 state variables:
## Land-cover / Forest species
## Biomass (m2/ha x 10 for forest, kg/ha * 10 for shrub)
## Forest species age 
## Time since last disturbance
## Type of last disturbance
######################################################################################

read.state.vars <- function(year){
  
  library(raster)
  
  cat("Read state variables", "\n")
  
  ## Read initial state varsas
  LCF <- raster(paste0("inputlyrs/asc/ForestMapSpp", year, "_31N-ETRS89.asc"))
  TSDIST <- raster(paste0("inputlyrs/asc/TimeSinceFire", year, "_31N-ETRS89.asc"))
  
  ## Build data frame
  land <- data.frame(cell.id=1:ncell(LCF), spp=LCF[], tsdist=TSDIST[])
  land <- land[!is.na(land$spp),]
  land$tsdist[is.na(land$tsdist)] <- 200
  
  ## Save it
  save(land, file="inputlyrs/rdata/land.rdata")
  
  ## MASK of the study area
  MASK <- LCF
  MASK[!is.na(MASK[])] <- 1
  crs(MASK) <- CRS("+init=epsg:25831")
  save(MASK, file="inputlyrs/rdata/mask.rdata") 
  
}


