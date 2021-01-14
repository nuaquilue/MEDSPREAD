######################################################################################
###  fire.regime()
###
######################################################################################

wildfires <- function(land, ignis, coord, orography, t, MASK, facc, rpb, fire.intens.th, print.maps, 
                      irun, pb.lower.th, pb.upper.th, fuel.opt, wwind, wslope){
                        
  cat("Fires", "\n") 
  
  ## Function to select items not in a vector
  `%notin%` <- Negate(`%in%`)
  
  ## To plot MAP of fire.ids and fire.step
  MAP <- MASK
  map <- data.frame(cell.id=land$cell.id, id=NA, step=NA)
  
  ## Read and load input data (això podria estar a land.dyn.mdl)
  fire.ignis <- read.table(paste0("inputfiles/", file.fire.ignis, ".txt"), header=T)
  fst.sprd.weight <- read.table(paste0("inputfiles/", file.sprd.weight, ".txt"), header=T)
  spp.ages <- read.table("inputfiles/SppAges.txt", header=T)
  shrub.fuel <- read.table("inputfiles/ShrubFuel.txt", header=T)
  
  ## To be sure that non-burnable covers do not burn (water, rock, urban)
  subland <- land[land$spp<=13,]          # 3.043.317
  suborography <- orography[land$spp<=13,]
  
  ## Reset TrackFires data frame each run 
  track.fire <- data.frame(fire.id=NA, fst=NA, wind=NA, atarget=NA, aburnt.highintens=NA, aburnt.lowintens=NA)
  track.burnt.cells <- data.frame(fire.id=NA, cell.id=NA, igni=NA, fintensity=NA)
  track.step <- data.frame(fire.id=NA, step=NA, nneigh=NA, nneigh.in=NA, nburn=NA, ncell.ff=NA)
  track.sprd <- data.frame(fire.id=NA, step=NA, cell.id=NA, dif.elev=NA, slope=NA, dif.wind=NA, wind=NA, fuel=NA, sr=NA, fi=NA, pb=NA, burn=NA)  
  
  ## Fuel; Fuel values of GrassCrop, Shrubs, Saplings, YoungForest, MatureForest
  fuel.val <- rbind(data.frame(opt="A", x=c(0.2,0.9,0.2,0.4,0.95)),
                    data.frame(opt="B", x=c(0.1,0.9,0.2,0.4,0.4,0.95,0.95)),
                    data.frame(opt="C", x=c(0.1,1,0.2,0.4,04,0.95,0.95)),
                    data.frame(opt="D", x=c(0.1,0.9,0.2,0.5,0.3,1,0.8)), # GC, SH, Sapling, ConifYoung, DecidYoung, ConifMature, DecidMature
                    data.frame(opt="E", x=c(0.1,0.9,0.2,0.6,0.3,1,0.7)), # GC, SH, Sapling, ConifYoung, DecidYoung, ConifMature, DecidMature
                    data.frame(opt="F", x=c(0.1,0.9,0.2,0.8,0.3,1,0.5)), # GC, SH, Sapling, ConifYoung, DecidYoung, ConifMature, DecidMature
                    data.frame(opt="G", x=c(0.1,0.9,0.2,0.6,0.3,1.6,0.8)), # GC, SH, Sapling, ConifYoung, DecidYoung, ConifMature, DecidMature
                    data.frame(opt="H", x=c(0.1,0.9,0.2,2,0.3,4,0.6))) # GC, SH, Sapling, ConifYoung, DecidYoung, ConifMature, DecidMature
  fuel.val <- fuel.val$x[fuel.val$opt==fuel.opt]
  aux <- left_join(subland, spp.ages, by="spp") %>% 
         mutate(fuel=ifelse(spp %in% c(11,12,13), fuel.val[1], # grass+crop
                       ifelse(spp==10, fuel.val[2],   ## shrub, max in MEDFIRE is 0.89
                          ifelse(tsdist<=7, fuel.val[3],  # saplings
                              ifelse(tsdist<=young & spp %in% c(1:4,8), fuel.val[4],  # conifer young
                                 ifelse(tsdist<=young & spp %in% c(5:7,9), fuel.val[5], # decid young
                                     ifelse(tsdist<=young & spp %in% c(1:4,8), fuel.val[6], # conifer mature
                                          fuel.val[7])))))))  # decidous  mature
  aux2 <- filter(aux, spp==10, tsdist<=24) %>% select(-fuel) %>%
          left_join(shrub.fuel, by="tsdist")
  aux$fuel[aux$cell.id %in% aux2$cell.id] <- aux2$fuel
  subland$fuel <- aux$fuel
  rm(aux); rm(aux2)

  
  ## Start with the 8 neigbours of the ignition
  ## Wind direction is coded as 0-N, 45-NE, 90-E, 135-SE, 180-S, 225-SW, 270-W, 315-NE
  default.neigh <- data.frame(x=c(-1,1,2900,-2900,2899,-2901,2901,-2899),
                              windir=c(270,90,180,0,225,315,135,45),
                              dist=c(100,100,100,100,141.421,141.421,141.421,141.42))
  default.nneigh <- nrow(default.neigh)

  
  ## Fire ids
  ids <- unlist(filter(fire.ignis, year==t) %>% select(fire.id))
  if(is_empty(ids))
    return(list(track.fire=track.fire[-1,], track.burnt.cells=track.burnt.cells[-1,], track.step=track.step[-1,],
              track.sprd=track.sprd[-1,]))
  
  
  ## Start burning the fires of the current time step
  burnt.cells <- visit.cells <- numeric()
  for(i in ids){
    
    ## cell.id of the ignition point
    igni.id <- unlist(filter(ignis, id==i) %>% select(cell.id))
    step <- 1
    
    ## Get the target area, the fire spread type and the wind direction
    fire.size.target <- unlist(filter(fire.ignis, fire.id==i) %>% select(area))
    fire.spread.type <- unlist(filter(fire.ignis, fire.id==i) %>% select(fst))
    fire.wind <- unlist(filter(fire.ignis, fire.id==i) %>% select(wind))
      
    # ## According to the fire spread type, look at the weights of each factor on spread rate.
    # ## Any other parameter depends on the fire.spread.type ¿¿
    # wwind <- fst.sprd.weight[1,fire.spread.type+1]
    # wslope <- fst.sprd.weight[2,fire.spread.type+1]

    ## Controls of the fire shape
    ## Max number of cells in the fire front
    mx.ncell.ff <- ifelse(fire.size.target<=500, 12, ifelse(fire.size.target<=1500, 20, ifelse(fire.size.target<=5000, 30, 40)))
    ## Min number of cells in the fire front, no sé si per incendis de me´s de 5000 o me´s d 10000
    mn.ncell.ff <- ifelse(fire.size.target<=5000, 8, 16)  # not sure if 16 for wind and 12 for convective
    ## wnsource --> per seguir direccionalitat vent
    wnsource <- ifelse(fire.spread.type==1, 100, 1)
    # threshodl ratio.burnt to be aplied, no sé si per incendis de me´s de 5000 o me´s d 10000
    thruky <- ifelse(fire.size.target<=5000, 0.85, 0.95)
    
    ## Initialize tracking variables
    ## Ignition always burnt, and it does in high intensity when no-PB
    fire.front <- igni.id
    cumul.source <- 1  # 0 si nsource=a+b
    aburnt.highintens <- 1
    aburnt.lowintens <- 0
    visit.cells <- c(visit.cells, igni.id)
    burnt.cells <- c(burnt.cells, igni.id) 
    map$id[map$cell.id==igni.id] <- i
    map$step[map$cell.id==igni.id] <- step
    track.burnt.cells <- rbind(track.burnt.cells, data.frame(fire.id=i, cell.id=igni.id, igni=T, fintensity=1))
    track.sprd <- rbind(track.sprd, data.frame(fire.id=i, step=1, cell.id=igni.id, dif.elev=0, slope=0, dif.wind=0, wind=0, fuel=0, sr=0, fi=0, pb=1, burn=T))
    track.step <- rbind(track.step, data.frame(fire.id=i, step=1, nneigh=1, nneigh.in=1, nburn=1, ncell.ff=1))
      
    ## Start speading from active cells (i.e. the fire front)
    while((aburnt.highintens+aburnt.lowintens)<fire.size.target){
      
      ## Increment step
      step=step+1
      
      ## Build a data frame with the theoretical 12 (=default.nneigh) neighbours of cells in fire.front, 
      ## and add the per definition wind direction and the distance.
      ## Filter cells thathave not been visited yet.
      neigh.id <- data.frame(cell.id=as.integer(rep(fire.front, each=default.nneigh)+
                                                rep(default.neigh$x, length(fire.front))),
                             source.id=rep(fire.front, each=default.nneigh),
                             position=rep(cumul.source, each=default.nneigh),
                             dist=rep(default.neigh$dist, length(fire.front)),
                             windir=rep(default.neigh$windir, length(fire.front))) %>%
                  filter(cell.id %notin% burnt.cells)
      
      ## Now find those neighbours that are currenty in Catalonia
      ## is_inCpp returns the position of neigh.id$cell.id in the 'land' data.frame (not the cell.id)!
      neighs.in.land <- is_inCpp(neigh.id$cell.id, subland$cell.id)
      i.land.in.neigh <- unique(neighs.in.land[which(neighs.in.land!=-1)])
      ## If all the available neighbours are out of Catalonia, stop spreading
      if(length(i.land.in.neigh)==0)
        break
      
      ## Retrive the current neighs
      neigh.land <- subland[i.land.in.neigh,] 
      
      ## Now, add to i.land.in.neigh, the indexes (positions) of fire.front cells.
      ## Further on, we'll need to know the elevation of the fire.front cells.
      i.land.in.neigh <- c(i.land.in.neigh, is_inCpp(fire.front, subland$cell.id)) 
      
      ## Retrieve the orography variables for fire.front and neigbhour cells, and compute aspect factor
      neigh.orography <- suborography[i.land.in.neigh,]
      
      ## Get spread rate by:
      ## Joining to the neig.id data.frame the neigh.land and keep only burnable neighs 
      ## Joining to this df, the neigh.orography to get the elevation of the source cells
      ## Joining to this df, the neigh.orography to get the elevation of the neighbour cells
      ## Computing slope and wind factors
      sprd.rate <- left_join(neigh.land, neigh.id, by="cell.id") %>%
                   left_join(select(neigh.orography, cell.id, elev), by=c("source.id"="cell.id")) %>%
                   left_join(select(neigh.orography, cell.id, elev), by="cell.id") %>% 
                   mutate(dif.elev = elev.y-elev.x, 
                          dif.wind = abs(windir-fire.wind),
                          slope = pmax(pmin(dif.elev,0.5),-0.5)+0.5,  
                          wind = ifelse(dif.wind==0, 0, ifelse(dif.wind %in% c(45,315), 0.25, 
                                   ifelse(dif.wind %in% c(90,270), 0.5, ifelse(dif.wind %in% c(135,225), 0.75, 1)))) ) %>% 
                   mutate(sr=wslope*slope+wwind*wind, fi=sr*fuel)
      sprd.rate$pb <- 1-exp(-facc*sprd.rate$fi) + runif(nrow(sprd.rate), -rpb, rpb)   
      sprd.rate <- group_by(sprd.rate, cell.id) %>% 
                   summarize(fire.id=i, dif.elev=max(dif.elev), slope=max(slope), dif.wind=max(dif.wind), wind=max(wind),  
                      fuel=max(fuel), sr=max(sr), fi=max(fi),  pb=max(pb), a=length(i), b=sum(position), nsource=b) 
      
      
      ## Now compute actual burning state (T or F) according to pb:
      sprd.rate$burn <- sprd.rate$pb >= runif(nrow(sprd.rate), pb.lower.th, pb.upper.th)
      if(nrow(sprd.rate)>0)
        track.sprd <- rbind(track.sprd,
                            data.frame(fire.id=i, step=step, cell.id=sprd.rate$cell.id, dif.elev=sprd.rate$dif.elev,
                                       slope=sprd.rate$slope, dif.wind=sprd.rate$dif.wind, wind=sprd.rate$wind, 
                                       fuel=sprd.rate$fuel, 
                                       sr=sprd.rate$sr, fi=sprd.rate$fi, pb=sprd.rate$pb, burn=sprd.rate$burn))
      
    
      ## Mark that all these neighs have been visited (before breaking in case no burning)
      visit.cells <- c(visit.cells, sprd.rate$cell.id)
      burnt.cells <- c(burnt.cells, sprd.rate$cell.id[sprd.rate$burn])
      
      ## If at least there's a burning cell, continue, otherwise, stop
      if(!any(sprd.rate$burn))
        break
      
      ## Avoid fire overshooting at last iteration: Only burn cells with higher pb
      temp.burnt <- sprd.rate[sprd.rate$burn, c("cell.id", "pb")]
      if((aburnt.lowintens+aburnt.highintens+nrow(temp.burnt))>fire.size.target){
        max.burnt <- fire.size.target - (aburnt.lowintens+aburnt.highintens)
        temp.burnt <- temp.burnt[order(temp.burnt$pb, decreasing = TRUE),]
        def.burnt <- temp.burnt$cell.id[1:max.burnt]
        sprd.rate$burn <- (sprd.rate$cell.id %in% def.burnt)
      }
      
      
      ## Mark the burnt cells with fire.id, mark the fire.step in burnt areas
      map$id[map$cell.id %in% sprd.rate$cell.id[sprd.rate$burn]] <- i
      map$step[map$cell.id %in% sprd.rate$cell.id[sprd.rate$burn]] <- step
      
      ## Track the burnt cells, the suppressed, and the fire intensity for burnt cells
      track.burnt.cells <- rbind(track.burnt.cells, 
                                 data.frame(fire.id=i, cell.id=sprd.rate$cell.id[sprd.rate$burn], 
                                            igni=F, fintensity=sprd.rate$fi[sprd.rate$burn]))
      
      
      ## Increase area burnt 
      aburnt.lowintens <- aburnt.lowintens + sum(sprd.rate$burn & sprd.rate$fi<=fire.intens.th)
      aburnt.highintens <- aburnt.highintens + sum(sprd.rate$burn & sprd.rate$fi>fire.intens.th)
      
      ## Select the new fire front
      ## First, count the number of cells burnt in the current step
      nburn <- sum(sprd.rate$burn)
      ## If any cell has burnt in the current step, stop
      if(nburn==0)
        break
      ## Otherwise, select the new fire front 
      else if(nburn<=mn.ncell.ff){
        fire.front <- sprd.rate$cell.id[sprd.rate$burn]
        cumul.source <- sprd.rate$nsource[sprd.rate$burn]
      }
      else{
        ratio.burnt <- (aburnt.lowintens+aburnt.highintens)/fire.size.target
        z <- rdunif(1,mx.ncell.ff-5,mx.ncell.ff)
        ncell.ff <- min(nburn*runif(1,0.5,0.7), z, na.rm=T)
        # si el nombre cell del ff coincideix amb el màxim  
        # o bé aleatòriament cap al final de l'incendi, forço compacitat.
        if(ncell.ff==z | (ratio.burnt>=thruky & runif(1,0,1)>=0.75))
          fire.front <- sort(sample(sprd.rate$cell.id[sprd.rate$burn], round(ncell.ff), replace=F,
                                    prob=sprd.rate$nsource[sprd.rate$burn]/100 ) )  
        else
          fire.front <- sort(sample(sprd.rate$cell.id[sprd.rate$burn], round(ncell.ff), replace=F, 
                                    prob=wnsource^sprd.rate$pb[sprd.rate$burn]) )
        cumul.source <- sprd.rate$nsource[sprd.rate$cell.id %in% fire.front]
      }
      
      ## In the case, there are no cells in the fire front, stop trying to burn.
      ## This happens when no cells have burnt in the current spreading step
      if(length(fire.front)==0)
        break
      
      ## Track each spreading step
      track.step <- rbind(track.step, data.frame(fire.id=i, step, nneigh=nrow(neigh.id),
                                                 nneigh.in=length(i.land.in.neigh), nburn, ncell.ff=length(fire.front)))
      
    } # while 'fire.size.target'
    
    ## Write info about this fire
    track.fire <- rbind(track.fire, data.frame(fire.id=i, fst=fire.spread.type, wind=fire.wind, atarget=fire.size.target, 
                                               aburnt.highintens,  aburnt.lowintens))
    
    cat(paste("Year:", t+1988, "- Fire:", i, "- aTarget:", fire.size.target, "- aBurnt:", aburnt.highintens+aburnt.lowintens),
        "- pBurnt:", round(100*(aburnt.highintens+aburnt.lowintens)/fire.size.target), "\n")
    
  }  # for 'ids'
  
  if(print.maps){
    # fire.ids
    MAP[!is.na(MASK[])] <- map$id
    writeRaster(MAP, paste0(out.path, "/FireIds_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
    ## fire.step
    MAP[!is.na(MASK[])] <- map$step
    writeRaster(MAP, paste0(out.path, "/FireStep_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
  }
  save(map, file=paste0(out.path, "/Maps_r", irun, "t", t, ".rdata"))
  
  return(list(track.fire=track.fire[-1,], track.burnt.cells=track.burnt.cells[-1,], track.step=track.step[-1,],
              track.sprd=track.sprd[-1,]))
}

