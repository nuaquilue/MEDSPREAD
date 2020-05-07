######################################################################################
##
######################################################################################

post.fire <- function(land, coord, orography, burnt.cells){
  
  ## Tracking
  cat("Post-fire regeneration", "\n") 
  
  ## Read matrix of secondary species according to species - sqi
  ## and fire response trait per species
  secondary.spp <- read.table("inputfiles/SecondarySpp.txt", header=T)
  response.trait <- read.table("inputfiles/FireResponseTrait.txt", header=T)
  spp.ages <- read.table("inputfiles/SppAges.txt", header=T)
  
  
  ## Num of neighbours in a circular neighbourhood according to radius (radius is in pixels)
  ## Assume that the neighbourhood is a star, with the maximum number of pixels in the
  ## east-west or north-south direction is 2*radius + 1 (1 is the center cell).
  ## The num of pixels is sequentially: 3+1*2, 5+3*2+1*2, 7+5*2+3*2+1*2, ...
  nneigh <- seq(3,41,2) + cumsum(seq(1,40,2)*2)

  ## Coordinates of high-intensity burnt forest cells the current time step that
  ## - it doesn't regenerate per se (fire functional trait), or
  ## - it is out of its climatic range, or
  ## - it is younger than the regeneration age
  burnt.land <- filter(land, cell.id %in% burnt.cells & spp<9) %>% left_join(response.trait, by="spp") %>%
                 left_join(spp.ages, by="spp") %>% filter(tsdist<=regener | trait==0) %>% 
                 left_join(coord, by = "cell.id") 
  
  ## Only continue if there's any cell with change of spp dominance
  if(nrow(burnt.land)>0){
    ## Coordinates of their closest neighbours (do not count for the cell itself)
    neigh.id <- nn2(coord[,-1], select(burnt.land,x,y),  searchtype="priority", k=nneigh[spp.distrib.rad])
    neigh.id <- neigh.id$nn.idx
    neigh.spp <- data.frame(cell.id=coord$cell.id[neigh.id[,1]],
                            matrix(land$spp[neigh.id[,-1]], nrow=nrow(neigh.id), ncol=ncol(neigh.id)-1) )
    
    ## Count number of neighbors per spp, assume that always there's a shrub cell in the neighbourhood
    neigh.spp <- data.frame(cell.id=coord$cell.id[neigh.id[,1]],
                            t(apply(neigh.spp[,-1], 1, count.spp))>=1 )
    neigh.spp$X10 <- T
    
    ## For those cells that a transition must be done:
    ## Look up sqi data and sencondary species  (according to dominant spp and sqi), 
    ## then add sdm of all tree species and finally
    ## add the number of forest spp in the neighbourhood
    burnt.land <- left_join(burnt.land, secondary.spp, by="spp") %>% 
                   left_join(neigh.spp, by = "cell.id")
  
    ## Select spp among available
    new.cohort <- data.frame(cell.id=burnt.land$cell.id,
                             spp=apply(select(burnt.land, phalepensis:shrub) * 
                                        select(burnt.land, X1:X10), 1, select.cohort) )
    
    return(new.cohort)  
  }
  
  else
    return(numeric())
}

