######################################################################################
## Build land.rdata with the initialization of the 6 state variables:
## Land-cover / Forest species
## Biomass (m2/ha x 10 for forest, kg/ha * 10 for shrub)
## Forest species age 
## Time since last disturbance
## Type of last disturbance
######################################################################################

read.state.vars <- function(work.path){
  
  library(raster)
  
  cat("Reading initial state variables", "\n")
  
  ## Read initial state varsas
  LCF <- raster(paste0(work.path, "/inputlyrs/asc/LCFspp_100m_31N-ETRS89.asc"))
  TSDIST <- raster(paste0(work.path, "/inputlyrs/asc/TSDisturb_100m_31N-ETRS89.asc"))
  
  ## Build data frame
  land <- data.frame(cell.id=1:ncell(LCF), spp=LCF[], tsdist=TSDIST[])
  land <- land[!is.na(land$spp),]
  
  ## Save it
  save(land, file="inputlyrs/rdata/land.rdata")
   
}