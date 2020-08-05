######################################################################################
###  fire.regime()
###
######################################################################################

fire.regime <- function(land, ignis, coord, orography, t){
                        
  cat("Fires", "\n") 
  
  ## Function to select items not in a vector
  `%notin%` <- Negate(`%in%`)
  
  ## Read and load input data
  fire.ignis <- read.table(paste0("inputfiles/", file.fire.ignis, ".txt"), header=T)
  spp.flammability <- read.table("inputfiles/SppFlammability.txt", header=T)
  spp.flammability.large <- read.table("inputfiles/SppFlammabilityLarge.txt", header=T)
  fst.sprd.weight <- read.table(paste0("inputfiles/", file.sprd.weight, ".txt"), header=T)
  spp.ages <- read.table("inputfiles/SppAges.txt", header=T)
  shrub.fuel <- read.table("inputfiles/ShrubFuel.txt", header=T)
  
  ## To be sure that non-burnable covers do not burn (water, rock, urban)
  i <- land$spp<=13        # 3.043.317
  subland <- land[i,]
  suborography <- orography[i,]
  
  ## Reset TrackFires data frame each run 
  track.fire <- data.frame(year=NA, fire.id=NA, fst=NA, wind=NA, atarget=NA, aburnt=NA)
  track.sprd <- data.frame(year=NA, fire.id=NA, step=NA, cell.id=NA, slope=NA, wind=NA,
                           flam=NA, aspc=NA, fuel=NA, sr=NA, fi=NA, pb=NA, burn=NA)
  
  ## Wind direction between 12 neigbours
  ## Wind direction is coded as 0-N, 45-NE, 90-E, 135-SE, 180-S, 225-SW, 270-W, 315-NE
  default.neigh <- data.frame(x=c(-1,1,2900,-2900,2899,-2901,2901,-2899,-2,2,5800,-5800),
                              windir=c(270,90,180,0,225,315,135,45,270,90,180,0),
                              dist=c(100,100,100,100,141.421,141.421,141.421,141.421,200,200,200,200))
  default.nneigh <- nrow(default.neigh)
  
  ## Fuel
  aux <- left_join(subland, spp.ages, by="spp") %>% 
         mutate(fuel=ifelse(spp %in% c(11,12,13), 0.2,
                       ifelse(spp==10, 1,
                          ifelse(tsdist<=7, 0.2,
                              ifelse(tsdist<young, 0.4, 0.95)))))
  aux2 <- filter(aux, spp==10, tsdist<=24) %>% select(-fuel) %>%
          left_join(shrub.fuel, by="tsdist")
  aux$fuel[aux$cell.id %in% aux2$cell.id] <- aux2$fuel
  subland$fuel <- aux$fuel
  rm(aux); rm(aux2)
  # if(any(is.na(subland$fuel)))
  #   print(filter(subland, is.na(fuel)))
  
  ## Fire ids
  ids <- unlist(filter(fire.ignis, year==t) %>% select(fire.id))
  if(is_empty(ids))
    return(list(burnt.cells=numeric(), fire.ids=numeric(), track.fire=track.fire[-1,], track.sprd=track.sprd[-1,]))
  
  ## Start burning the fires of the current time step
  burnt.cells <- visit.cells <- fire.ids <- numeric()
  for(i in ids){
    
    ## cell.id of the ignition point
    igni.id <- unlist(filter(ignis, id==i) %>% select(cell.id))
    step=1
    
    ## Get the target area, the fire spread type and the wind direction
    fire.size.target <- unlist(filter(fire.ignis, fire.id==i) %>% select(area))
    fire.spread.type <- unlist(filter(fire.ignis, fire.id==i) %>% select(fst))
    fire.wind <- unlist(filter(fire.ignis, fire.id==i) %>% select(wind))
      
    ## According to the fire spread type, look at the weights of each factor on spread rate,
    ## and the species flammability
    rpb <- fst.sprd.weight[1,fire.spread.type+1]
    wwind <- fst.sprd.weight[2,fire.spread.type+1]
    wslope <- fst.sprd.weight[3,fire.spread.type+1]
    wflam <- fst.sprd.weight[4,fire.spread.type+1]
    waspc <- fst.sprd.weight[5,fire.spread.type+1]
    if(fire.size.target>2000)
      spp.flam <- filter(spp.flammability.large, fst==fire.spread.type) %>% select(-fst)
    else
      spp.flam <- filter(spp.flammability, fst==fire.spread.type) %>% select(-fst)
    fi.acc <- ifelse(fire.size.target>2000, fi.accelerate, 1)
    
    ## Initialize tracking variables
    ## Ignition always burnt, and it does in high intensity when no-PB
    fire.front <- igni.id
    aburnt <- 1
    burnt.cells <- c(burnt.cells, igni.id)
    visit.cells <- c(visit.cells, igni.id) # to account for visit (and burnt) cells 
    fire.ids <- c(fire.ids, i)
      
    ## Start speading from active cells (i.e. the fire front)
    while(aburnt<fire.size.target){
      
      ## Build a data frame with the theoretical 12 (=default.nneigh) neighbours of cells in fire.front, 
      ## and add the per definition wind direction and the distance.
      ## Filter cells thathave not been visited yet.
      neigh.id <- data.frame(cell.id=as.integer(rep(fire.front, each=default.nneigh)+
                                                rep(default.neigh$x, length(fire.front))),
                             source.id=rep(fire.front, each=default.nneigh),
                             dist=rep(default.neigh$dist,length(fire.front)),
                             windir=rep(default.neigh$windir,length(fire.front)) ) %>%
                  filter(cell.id %notin% visit.cells)
      
      ## Now find those neighbours that are currenty in Catalonia
      ## is_inCpp returns the position of neigh.id$cell.id in the 'land' data.frame (not the cell.id)!
      neigh.in.land <- is_inCpp(neigh.id$cell.id, subland$cell.id)
      i.land.in.neigh <- unique(neigh.in.land[which(neigh.in.land!=-1)])
      
      ## For all neighbours, compute fire intenstiy and flammability factors
      ## fire intenstiy and flam will be NA for non burnable covers
      neigh.land <- subland[i.land.in.neigh,] %>% 
                    left_join(spp.flam, by="spp") %>% mutate(flam=wflam*flam) 
      
      ## Now, add to i.land.in.neigh, the indexes (positions) of fire.front cells.
      ## Further on, we'll need to know the elevation of the fire.front cells.
      i.land.in.neigh <- c(i.land.in.neigh, is_inCpp(fire.front, subland$cell.id)) 
      
      ## Retrieve the orography variables for fire.front and neigbhour cells, 
      ## and already compute aspect factor
      neigh.orography <- suborography[i.land.in.neigh,] %>%
                         mutate(aspc=waspc*ifelse(aspect==1, 0.1, 
                                                  ifelse(aspect==3, 0.9, ifelse(aspect==4, 0.4, 0.3))))
      
      ## Get spread rate by:
      ## Joining to the neig.id data.frame the neigh.land and keep only burnable neighs 
      ## Joining to this df, the neigh.orography to get the elevation of the source cells
      ## Joining to this df, the neigh.orography to get the elevation of the neighbour cells
      ## Computing slope and wind factors
      sprd.rate <- left_join(neigh.land, neigh.id, by="cell.id") %>%
                   left_join(select(neigh.orography, cell.id, elev), by=c("source.id"="cell.id")) %>%
                   left_join(select(neigh.orography, cell.id, elev, aspc), by="cell.id") %>% 
                   mutate(dif.elev = elev.y-elev.x, 
                          slope = wslope * (pmax(pmin(dif.elev/dist,0.5),-0.5)+0.5), 
                          wind = wwind * (ifelse(abs(windir-fire.wind)>180, 
                                            360-abs(windir-fire.wind), abs(windir-fire.wind)))/180) %>% 
                   mutate(sr=slope+wind+flam+aspc, fi=sr*fuel*fi.acc, pb=1+rpb*log(sr*fuel*fi.acc)) %>% 
                   group_by(cell.id) %>% 
                   summarize(fire.id=i, slope=max(slope), wind=max(wind), flam=max(flam),
                             aspc=max(aspc), fuel=max(fuel), sr=max(sr), fi=max(fi), pb=max(pb)) 
      
      ## Now compute actual burning state (T or F) according to pb and suppress:
      sprd.rate$burn <- sprd.rate$pb >= runif(nrow(sprd.rate), pb.lower.th, pb.upper.th)
      if(nrow(sprd.rate)>0)
      track.sprd <- rbind(track.sprd,
                          data.frame(year=t, fire.id=i, step=step, cell.id=sprd.rate$cell.id,
                                     slope=sprd.rate$slope, wind=sprd.rate$wind, flam=sprd.rate$flam, 
                                     aspc=sprd.rate$aspc, fuel=sprd.rate$fuel, sr=sprd.rate$sr, 
                                     fi=sprd.rate$fi, pb=sprd.rate$pb, burn=sprd.rate$burn))
      
      
      ## Mark that all these neighs have been visited (before breaking in case no burning)
      visit.cells <- c(visit.cells, sprd.rate$cell.id)
      
      ## Avoid fire overshooting at last iteration: Only burn cells with higher pb
      temp.burnt <- sprd.rate[sprd.rate$burn, c("cell.id", "pb")]
      if(aburnt+nrow(temp.burnt)>fire.size.target){
        max.burnt <- fire.size.target - aburnt
        temp.burnt <- temp.burnt[order(temp.burnt$pb, decreasing = TRUE),]
        def.burnt <- temp.burnt$cell.id[1:max.burnt]
        sprd.rate$burn <- (sprd.rate$cell.id %in% def.burnt)
      }
      
      ## If at least there's a burning cell, continue, otherwise, stop
      if(nrow(sprd.rate)==0)
        break
      
      ## Mark the burnt cells and the fire intensity 
      burnt.cells <- c(burnt.cells, sprd.rate$cell.id[sprd.rate$burn])
      fire.ids <- c(fire.ids, sprd.rate$fire.id[sprd.rate$burn])
      
      # ## Select the new fire front
      # ## First, count the number of cells burnt in the current step
      nburn <- sum(sprd.rate$burn)
      ## If any cell has burnt in the current step, stop
      if(nburn==0)
        break
      ## Otherwise, select the new fire front, a random number from 1 to n.cell.burnt 
      ## according to fire.intensity
      if(nburn==1)
        fire.front <- sprd.rate$cell.id[sprd.rate$burn]
      if(nburn>1){
        if(fire.size.target>2000 & aburnt/fire.size.target<=0.75)
          fire.front <- sprd.rate$cell.id[sprd.rate$burn]
        else  
          fire.front <- base::sample(sprd.rate$cell.id[sprd.rate$burn], rdunif(1, 1, round(nburn*3/4)),
                              replace=F, prob=sprd.rate$fi[sprd.rate$burn]*runif(nburn, 0.65, 1))  
      }
            # exclude.th <- min(max(sprd.rate$sr)-0.005,   ## 'mad' -> median absolute deviation
            #                   rnorm(1,mean(sprd.rate$sr[sprd.rate$burn], na.rm=T)-mad(sprd.rate$sr[sprd.rate$burn]/2, na.rm=T),
            #                         mad(sprd.rate$sr[sprd.rate$burn], na.rm=T)))
            # fire.front <- sprd.rate$cell.id[sprd.rate$burn & sprd.rate$sr>=exclude.th]
      
      ## Increase area burnt in either high or low intensity (Prescribed burns always burnt in low intensity)
      aburnt <- aburnt + sum(sprd.rate$burn)
      
      ## In the case, there are no cells in the fire front, stop trying to burn.
      ## This happens when no cells have burnt in the current spreading step
      if(length(fire.front)==0)
        break
      
      step=step+1
    } # while 'fire.size.target'
    
    ## Write info about this fire
    track.fire <- rbind(track.fire, data.frame(year=t, fire.id=i, fst=fire.spread.type, 
                                               wind=fire.wind, atarget=fire.size.target, aburnt))
    
    cat(paste("Fire:", i, "- aTarget:", fire.size.target, "- aBurnt:", aburnt), "\n")
    
  }  # for 'ids'
  
  return(list(burnt.cells=burnt.cells, fire.ids=fire.ids, track.fire=track.fire[-1,],
         track.sprd=track.sprd[-1,]))
}

