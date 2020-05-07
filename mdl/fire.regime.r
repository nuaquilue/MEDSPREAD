######################################################################################
###  fire.regime()
###
######################################################################################

fire.regime <- function(land, ignis, coord, orography, t){
                        
  cat("Fires", "\n") 
  
  ## Function to select items not in a vector
  `%notin%` <- Negate(`%in%`)
  
  ## Read and load input data
  fire.ignis <- read.table("inputfiles/FireIgnitions.txt", header=T)
  spp.flammability <- read.table("inputfiles/SppSpreadRate.txt", header=T)
  fst.sprd.weight <- read.table(paste0("inputfiles/", file.sprd.weight, ".txt"), header=T)
  
  ## To be sure that non-burnable covers do not burn (water, rock, urban)
  i <- land$spp<=13        # 3.043.317
  subland <- land[i,]
  suborography <- orography[i,]
  
  ## Reset TrackFires data frame each run 
  track.fire <- data.frame(year=NA, fire.id=NA, fst=NA, wind=NA, atarget=NA, aburnt=NA)
  
  ## Wind direction between 12 neigbours
  ## Wind direction is coded as 0-N, 45-NE, 90-E, 135-SE, 180-S, 225-SW, 270-W, 315-NE
  default.neigh <- data.frame(x=c(-1,1,2900,-2900,2899,-2901,2901,-2899,-2,2,5800,-5800),
                              windir=c(270,90,180,0,225,315,135,45,270,90,180,0),
                              dist=c(100,100,100,100,141.421,141.421,141.421,141.421,200,200,200,200))
  default.nneigh <- nrow(default.neigh)
  
  ## Fire ids
  ids <- unlist(filter(fire.ignis, year==t) %>% select(fire.id))
  if(is_empty(ids))
    return(list(burnt.cells=numeric(), fire.ids=numeric(), track.fire=track.fire[-1,]))
  
  ## Start burning the fires of the current time step
  burnt.cells <- visit.cells <- fire.ids <- numeric()
  for(id in ids){
    
    ## cell.id of the ignition point
    igni.id <- unlist(filter(ignis, fire.id==id) %>% select(cell.id))
    
    ## Get the target area, the fire spread type and the wind direction
    fire.size.target <- unlist(filter(fire.ignis, fire.id==id) %>% select(area))
    fire.spread.type <- unlist(filter(fire.ignis, fire.id==id) %>% select(fst))
    fire.wind <- unlist(filter(fire.ignis, fire.id==id) %>% select(wind))
      
    ## According to the fire spread type, look at the weights of each factor on spread rate,
    ## and the species flammability
    wwind <- fst.sprd.weight[1,fire.spread.type+1]
    wslope <- fst.sprd.weight[2,fire.spread.type+1]
    wflam <- fst.sprd.weight[3,fire.spread.type+1]
    waspc <- fst.sprd.weight[4,fire.spread.type+1]
    spp.flam <- filter(spp.flammability, fst==fire.spread.type) %>% select(-fst)
    
    ## Initialize tracking variables
    ## Ignition always burnt, and it does in high intensity when no-PB
    fire.front <- igni.id
    aburnt <- 1
    burnt.cells <- c(burnt.cells, igni.id)
    visit.cells <- c(visit.cells, igni.id) # to account for visit (and burnt) cells 
    fire.ids <- c(fire.ids, id)
      
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
      neigh.land <- subland[i.land.in.neigh,] %>% mutate(fuel=1) %>%
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
                          slope = wslope * pmax(pmin(dif.elev/dist,0.5),-0.5)+0.5, 
                          wind = wwind * (ifelse(abs(windir-fire.wind)>180, 
                                            360-abs(windir-fire.wind), abs(windir-fire.wind)))/180) %>% 
                   mutate(sr=slope+wind+flam+aspc, pb=1+rpb*log(sr*fuel))  %>%
                   group_by(cell.id) %>% 
                   summarize(fire.id=id, spp=mean(spp), sr=max(sr), pb=max(pb)) 
      
      ## Now compute actual burning state (T or F) according to pb and suppress:
      sprd.rate$burning <- (runif(nrow(sprd.rate),0,pb.upper.th)<=sprd.rate$pb & sprd.rate$pb>pb.lower.th)
      
      ## Mark that all these neighs have been visited (before breaking in case no burning)
      visit.cells <- c(visit.cells, sprd.rate$cell.id)
      
      ## If at least there's a burning cell, continue, otherwise, stop
      if(nrow(sprd.rate)==0)
        break
      
      ## Mark the burnt cells and the fire intensity 
      burnt.cells <- c(burnt.cells, sprd.rate$cell.id[sprd.rate$burning])
      fire.ids <- c(fire.ids, sprd.rate$fire.id[sprd.rate$burning])
      
      ## Select the new fire front
      exclude.th <- min(max(sprd.rate$sr)-0.005,   ## 'mad' -> median absolute deviation
                        rnorm(1,mean(sprd.rate$sr[sprd.rate$burning], na.rm=T)-mad(sprd.rate$sr[sprd.rate$burning], na.rm=T)/2,
                              mad(sprd.rate$sr[sprd.rate$burning], na.rm=T)))
      fire.front <- sprd.rate$cell.id[sprd.rate$burning & sprd.rate$sr>=exclude.th]
      
      ## Increase area burnt in either high or low intensity (Prescribed burns always burnt in low intensity)
      aburnt <- aburnt + sum(sprd.rate$burning)
      
      ## In the case, there are no cells in the fire front, stop trying to burn.
      ## This happens when no cells have burnt in the current spreading step
      if(length(fire.front)==0)
        break
      
    } # while 'fire.size.target'
    
    ## Write info about this fire
    track.fire <- rbind(track.fire, data.frame(year=t, fire.id=id, fst=fire.spread.type, 
                                               wind=fire.wind, atarget=fire.size.target, aburnt))
    
    cat(paste("Fire:", id, "- aTarget:", fire.size.target, "- aBurnt:", aburnt), "\n")
    
  }  # for 'ids'
  
  return(list(burnt.cells=burnt.cells, fire.ids=fire.ids, track.fire=track.fire[-1,]))
}

