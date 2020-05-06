######################################################################################
##
######################################################################################

land.dyn.mdl <- function(scn.name){
  
  ## Load required packages and functions 
  suppressPackageStartupMessages({
    library(tictoc)
    library(sp)
    library(raster)  
    library(RANN)  # for nn2()
    library(Rcpp)
    library(tidyverse)
  })
  source("mdl/growth.r")
  source("mdl/fire.regime.r")
  source("mdl/post.fire.r")
  sourceCpp("mdl/is.in.cpp")
  
  
  ## Load scenario definition (global variables and scenario parameters)
  ## and customized scenario parameters
  source(paste0("outputs/", scn.name, "/scn.def.r"))
  if(file.exists(paste0("outputs/", scn.name, "/scn.custom.def.r")))
    source(paste0("outputs/", scn.name, "/scn.custom.def.r"))
  
  
  ## Load:
  ## 1. Mask of the study area (raster)
  ## 2. Data frame with cell.id and coordinates x, y
  ## 3. Data frame of the model static variables 
  ## 4. Data frame with interface value
  load("inputlyrs/rdata/mask.rdata")
  load("inputlyrs/rdata/coordinates.rdata")
  load("inputlyrs/rdata/orography.rdata")
  
  
  ## Set the directory for writing spatial outputs (create it, if it does not exist yet) 
  if(write.sp.outputs){      
    if(!file.exists(paste0(out.path, "/lyr")))
      dir.create(file.path(getwd(), out.path, "/lyr"), showWarnings = F) 
  }

  
  ## List the name of the forest species
  species <- c("phalepensis", "pnigra", "ppinea", "psylvestris", "ppinaster", "puncinata",
               "aalba", "qilex", "qsuber", "qfaginea", "qhumilis", "fsylvatica", "other")
                  

  ## Tracking data.frames
  track.fire <-  data.frame(run=NA, year=NA, swc=NA, clim.sever=NA, fire.id=NA, fst=NA, wind=NA, atarget=NA, 
                            aburnt.highintens=NA, aburnt.lowintens=NA, asupp.fuel=NA, asupp.sprd=NA)
  track.fire.spp <-  data.frame(run=NA, year=NA, fire.id=NA, spp=NA, aburnt=NA, bburnt=NA)
  track.post.fire <- data.frame(run=NA, year=NA, spp.out=NA, Var2=NA, Freq=NA)
  track.afforest <- data.frame(run=NA, year=NA, Var1=NA, Freq=NA)
  track.land <- data.frame(run=NA, year=NA, spp=NA, area=NA, vol=NA, volbark=NA, carbon=NA)
  
  
  ## Start the simulations   
  for(irun in 1:nrun){
    
    ## Copy the schedulings in auxiliar vectors (only for those processes included in the current version)
    temp.fire.schedule <- seq(1, time.horizon, fire.step)
    temp.post.fire.schedule <- seq(1, time.horizon, post.fire.step)
    
    ## Load initial spatial dynamic state variables in a data.frame format
    load("inputlyrs/rdata/land.rdata")
    
    ## Start the discrete time sequence 
    for(t in time.seq){
      
      ## Track scenario, replicate and time step
      cat(paste0("scn: ", scn.name," - run: ", irun, "/", nrun, " - time: ", t, "/", time.horizon), "\n")
      
      ## FIRES
      ## Tracking variables to be re-initialized each time step
      ## Out of the "if(fires)" in case only prescribed burns are applied
      burnt.cells <- integer()
      fire.ids <- integer()
      id.fire <- annual.burnt <- 0
      fire.out <- fire.regime(land, coord, orography, t)
      burnt.cells <- fire.out[[1]]; fintensity <- fire.out[[2]]; 
      fire.ids <- fire.out[[3]]; id.fire <- id.fire+nrow(fire.out[[4]])
      # track fire events and total annual burnt area
      if(nrow(fire.out[[4]])>0)
        track.fire <- rbind(track.fire, data.frame(run=irun, fire.out[[4]]))
      # track spp and biomass burnt
      aux <- data.frame(cell.id=burnt.cells, fire.id=fire.ids, fintensity) %>% 
             left_join(select(land, cell.id, spp, biom), by="cell.id") %>%
             mutate(bburnt=ifelse(fintensity>fire.intens.th, biom, biom*(1-fintensity))) %>%
             group_by(fire.id, spp) %>% summarize(aburnt=length(spp), bburnt=round(sum(bburnt, na.rm=T),1))
      if(nrow(aux)>0)
        track.fire.spp <-  rbind(track.fire.spp, data.frame(run=irun, year=t, aux)) 
      # Done with fires! 
      temp.fire.schedule <- temp.fire.schedule[-1] 
      rm(fire.out); rm(aux)
      
      
      ## POST-FIRE REGENERATION
      ## Forest transition of tree species burnt in high intensity
      aux  <- post.fire(land, coord, orography)
      spp.out <- land$spp[land$cell.id %in% aux$cell.id]
      land$spp[land$cell.id %in% aux$cell.id] <- aux$spp
      track.post.fire <- rbind(track.post.fire, data.frame(run=irun, year=t, table(spp.out, aux$spp)))  
      land$tsdist[land$cell.id %in% burnt.cells] <- 0
      temp.post.fire.schedule <- temp.post.fire.schedule[-1] 
      
      
      ## AGING
      land$tsdist <- pmin(land$tsdist+1,600)
      
      
      ## Print maps every time step with ignition and low/high intenstiy burnt
      if(write.sp.outputs){
        MAP <- MASK
        cat("... writing output layers", "\n")
        nfire <- sum(track.fire$year==t, na.rm=T)
        sizes <- filter(track.fire, year==t) %>% group_by(swc, fire.id) %>% summarise(ab=aburnt.highintens+aburnt.lowintens)
        # Ignitions' cell.id 
        igni.id <- burnt.cells[c(1,cumsum(sizes$ab)[1:(nfire-1)]+1)] 
        MAP[!is.na(MASK[])] <- land$distype*(land$tsdist==1)
        MAP[igni.id] <- 9
        writeRaster(MAP, paste0(out.path, "/lyr/DistType_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
      }
      
      ## Deallocate memory
      gc(verbose=F)  
      cat("\n")
      
    } # time
  
  } # run
  
  cat("... writing outputs", "\n")
  track.fire$rem <- pmax(0,track.fire$atarget-track.fire$aburnt.highintens-track.fire$aburnt.lowintens)
  write.table(track.fire[-1,], paste0(out.path, "/Fires.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.fire.spp[-1,], paste0(out.path, "/FiresSpp.txt"), quote=F, row.names=F, sep="\t")
  names(track.post.fire)[4:5] <- c("spp.in", "ha")
  write.table(track.post.fire[-1,], paste0(out.path, "/PostFire.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.land[-1,], paste0(out.path, "/Land.txt"), quote=F, row.names=F, sep="\t")

}
