######################################################################################
##
######################################################################################

land.dyn.mdl <- function(scn.name){
  
  ## Load required packages and functions 
  suppressPackageStartupMessages({
    library(assertthat)
    library(raster)  
    library(RANN)
    library(Rcpp)
    library(tidyverse)
  })
  source("mdl/read.state.vars.r")
  source("mdl/fire.regime.r")
  source("mdl/post.fire.r")
  source("mdl/auxiliars.r")
  source("mdl/update.vars.r")
  sourceCpp("mdl/is.in.cpp")
  
  ## Load global variables and customized scenario parameters
  source(paste0("outputs/", scn.name, "/scn.def.r"))
  if(file.exists(paste0("outputs/", scn.name, "/scn.custom.def.r")))
    source(paste0("outputs/", scn.name, "/scn.custom.def.r"))
  
  ## Set the directory for writing spatial outputs (create it, if it does not exist yet) 
  if(write.sp.outputs){      
    if(!file.exists(paste0(out.path, "/lyr")))
      dir.create(file.path(getwd(), out.path, "/lyr"), showWarnings = F) 
  }
  
  ## Load ignitions
  load("inputlyrs/rdata/ignitions.rdata")
  
  ## Tracking data.frames
  track.fire <-  data.frame(run=NA, year=NA, fire.id=NA, fst=NA, wind=NA, atarget=NA, aburnt=NA)
  track.post.fire <- data.frame(run=NA, year=NA, spp.out=NA, Var2=NA, Freq=NA)
  track.sprd <- data.frame(run=NA, year=NA, fire.id=NA, step=NA, cell.id=NA, slope=NA, wind=NA,
                           flam=NA, aspc=NA, fuel=NA, sr=NA, fi=NA, pb=NA, burn=NA)
  
  ## Dataframe to record validation
  if(validation)
    result <- data.frame(scn=NA, run=NA, year=NA, fire.id=NA, am=NA)

  ## Set up time sequence
  time.seq <- seq(1, time.horizon, 1)
  
  ## Start the simulations   
  for(irun in 1:nrun){
  
    ## Load:
    ## 1. Mask of the study area (raster)
    ## 2. Data frame with cell.id and coordinates x, y
    ## 3. Data frame of the model static variables 
    ## 4. Data frame with dynamic state variables
    load("inputlyrs/rdata/mask.89-99.rdata")
    load("inputlyrs/rdata/coordinates.rdata")
    load("inputlyrs/rdata/orography.rdata")
    load("inputlyrs/rdata/land.rdata")
    if(validation)
      load("inputlyrs/rdata/fireperim.89-99.rdata")
    
    ## Schedule of fires and post-fire
    fire.schedule <- seq(1, time.horizon, 1)
    post.fire.schedule <- seq(1, time.horizon, 1)
    
    ## Start the discrete time sequence 
    for(t in time.seq){  
      
      ## Track scenario, replicate and time step
      cat(paste0("scn: ", scn.name," - run: ", irun, "/", nrun, " - time: ", t, "/", time.horizon), "\n")
      
      ## UPDATE LAND
      if(t==12){
        load("inputlyrs/rdata/mask.00-12.rdata")
        out <- update.vars(year="00", MASK)
        land <- out[[1]]
        coord <- out[[2]]
        orography <- out[[3]]
        ignis <- out[[4]]
        rm(out); gc(verbose=F)
        if(validation)
          load("inputlyrs/rdata/fireperim.00-12.rdata")
      }
      
      ## FIRES
      fire.out <- fire.regime(land, ignis, coord, orography, t)
      burnt.cells <- fire.out[[1]]; fire.ids <- fire.out[[2]] 
      # track fire events 
      if(nrow(fire.out[[3]])>0)
        track.fire <- rbind(track.fire, data.frame(run=irun, fire.out[[3]]))
      # track fire.spread
      if(nrow(fire.out[[4]])>0)
        track.sprd<- rbind(track.sprd, data.frame(run=irun, fire.out[[4]]))
              # # track spp burnt
              # aux <- data.frame(cell.id=burnt.cells, fire.id=fire.ids) %>%
              #        left_join(select(land, cell.id, spp), by="cell.id") %>%
              #        group_by(fire.id, spp) %>% summarize(aburnt=length(spp))
              # if(nrow(aux)>0)
              #   track.fire.spp <-  rbind(track.fire.spp, data.frame(run=irun, year=t, aux))
      # Done with fires!
      fire.schedule <- fire.schedule[-1] 
      rm(fire.out)
      
      ## POST-FIRE REGENERATION
      if(length(burnt.cells)>0){
        aux  <- post.fire(land, coord, orography, burnt.cells)
        if(!is.atomic(aux)){
          spp.out <- land$spp[land$cell.id %in% aux$cell.id]
          land$spp[land$cell.id %in% aux$cell.id] <- aux$spp
          track.post.fire <- rbind(track.post.fire, data.frame(run=irun, year=t, table(spp.out, aux$spp)))  
        }
      }
      post.fire.schedule <- post.fire.schedule[-1] 
      
      ## AGING
      cat("Aging", "\n")
      land$tsdist[land$cell.id %in% burnt.cells] <- 0
      land$tsdist <- pmin(land$tsdist+1,600)
      
      ## Print maps of fire.id in burnt locations
      if(write.sp.outputs & length(burnt.cells)>0){
        cat("... writing output layers", "\n")
        MAP <- MASK
        aux <- data.frame(cell.id=burnt.cells, fire.id=fire.ids)
        aux <- data.frame(cell.id=land$cell.id) %>% left_join(aux, by="cell.id")
        MAP[!is.na(MASK[])] <- aux$fire.id
        writeRaster(MAP, paste0("outputs/", scn.name, "/lyr/FireID_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
      }
      
      ## VALIDATION
      if(validation & length(burnt.cells)>0){
        perim.y <- perim.id[,c(1,t+1)]; names(perim.y)[2] <- "perim"
        ## Build fire layer and transform into a dataframe
        MAP <- MASK
        aux <- data.frame(cell.id=burnt.cells, fire.id=fire.ids)
        aux <- data.frame(cell.id=land$cell.id) %>% left_join(aux, by="cell.id")
        MAP[!is.na(MASK[])] <- aux$fire.id
        dta <- data.frame(cell.id=1:ncell(MASK[]), mask=MASK[], fire.id=MAP[])
        ## Match simulated and observed fires
        aux <- filter(dta, !is.na(mask)) %>% select(-mask) %>% left_join(perim.y, by="cell.id") %>%
          filter(!is.na(fire.id)) %>% mutate(dif=fire.id-perim) %>% filter(dif==0) %>% 
          group_by(fire.id) %>% summarize(am=length(fire.id))
        if(nrow(aux)>0)
          result <- rbind(result, data.frame(scn=scn.name, run=irun, year=t, aux))    
      }
      
      
      ## Deallocate memory
      gc(verbose=F)  
      cat("\n")
      
    } # time
  } # run
  
  cat("... writing outputs", "\n")
  track.fire$extra <- track.fire$atarget-track.fire$aburnt
  track.fire$pextra <- round(track.fire$extra/track.fire$atarget*100,1)
  write.table(track.fire[-1,], paste0(out.path, "/Fires.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.sprd[-1,], paste0(out.path, "/FiresSprd.txt"), quote=F, row.names=F, sep="\t")
  names(track.post.fire)[4:5] <- c("spp.in", "ha")
  write.table(track.post.fire[-1,], paste0(out.path, "/PostFire.txt"), quote=F, row.names=F, sep="\t")
  if(validation)
    write.table(result[-1,], paste0(out.path, "/Validation.txt"), quote=F, row.names=F, sep="\t")
}
