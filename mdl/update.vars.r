update.vars <- function(year, MASK){
  
  library(raster)
  
  cat("Update state variables", "\n")
  
  ## Read initial state varsas
  LCF <- raster(paste0("inputlyrs/asc/ForestMapSpp", year, "_31N-ETRS89.asc"))
  TSDIST <- raster(paste0("inputlyrs/asc/TimeSinceFire", year, "_31N-ETRS89.asc"))
  
  ## Build data frame
  land <- data.frame(cell.id=1:ncell(LCF), spp=LCF[], tsdist=TSDIST[])
  land <- land[!is.na(land$spp),]
  land$tsdist[is.na(land$tsdist)] <- 200
  
  ## Build a coordinates data frame from MASK raster
  coord <- data.frame(cell.id=1:ncell(MASK), coordinates(MASK), mask=MASK[])
  coord <- filter(coord, !is.na(mask)) %>% select(-mask)
  
  ## Read orography variables, build and save the data frame
  ELEVATION <- raster("inputlyrs/asc/DEM_100m_31N-ETRS89.asc")
  ASPECT <- raster("inputlyrs/asc/Aspect_100m_31N-ETRS89.asc")
  orography <- data.frame(cell.id=1:ncell(MASK), elev=ELEVATION[], aspect=ASPECT[])
  orography <- orography[!is.na(MASK[]),]
  orography$elev[is.na(orography$elev)] <- rnorm(sum(is.na(orography$elev)),
                                                 mean(orography$elev, na.rm=T), sd(orography$elev, na.rm=T))
  orography$aspect[is.na(orography$aspect)] <- sample(1:4, sum(is.na(orography$aspect)), replace=T)

  ## Read ignitions, build and save the data frame
  IGNIS <- raster("inputlyrs/asc/FireIgnis8912_31N-ETRS89.asc")
  ignis <- data.frame(cell.id=1:ncell(MASK), id=IGNIS[])
  ignis <- ignis[!is.na(MASK[]),]
  
    
  return(list(land=land, coord=coord, orography=orography, ignis=ignis))
  
}
