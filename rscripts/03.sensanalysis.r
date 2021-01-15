play <- function(){
  rm(list=ls())
  library(readxl)
  library(raster)
  library(tidyverse)
  
  source("rscripts/03.sensanalysis.r")
  
  ## List of scenarios
  list.scn <- list.dirs(path=paste0(getwd(), "/outputs"), full.name=F, recursive=F); list.scn
  
  ## Auto-extinction
  extinct <- self.extinguishing(list.scn[-1])
  write.table(extinct, "rscripts/outs/extinct_33scn.txt", quote=F, row.names=F, sep="\t")
  
  ## Percentage burnt per land-cover type
  lctburnt <- pct.lct.burnt(list.scn[-1]) 
  write.table(lctburnt[[2]], "rscripts/outs/error.pctburntlctfst_33scn.txt", quote=F, row.names=F, sep="\t")
  write.table(lctburnt[[1]], "rscripts/outs/error.pctburntlct_33scn.txt", quote=F, row.names=F, sep="\t")
  
  ## Match area
  report <- match.area(list.scn[-1])
  write.table(report, "rscripts/outs/matcharea_33scn.txt", quote=F, row.names=F, sep="\t")
  best <- best.scn(report)
  write.table(best, "rscripts/outs/bestmatcharea_33scn.txt", quote=F, row.names=F, sep="\t")
  
  ## Plot
  
  plot.fires(list.scn[23], 2, 24)
  
}



#################################################### AUTO EXTINCTION ####################################################
self.extinguishing <- function(list.scn){
  result <- data.frame(scn=NA, fst=NA, pctg.reach=NA, ab.at=NA)
  for(scn in list.scn){
    fires <- read.table(paste0("outputs/", scn, "/_Fires.txt"), header=T)
    # percentage of fires that have reach target area, per fst
    nfires.fst <- group_by(fires, fst) %>% summarise(n=length(fst))
    aux1 <- group_by(fires,fst) %>% summarise(pctg.reach=sum(pextra==0)) %>% 
      left_join(nfires.fst, by="fst") %>% mutate(pctg.reach=100*pctg.reach/n) %>% 
      select(-n)
    # percentage of area effectively burnt over the target area
    aux2 <- group_by(fires,fst) %>% 
      summarise(ab=sum(aburnt.highintens+aburnt.lowintens), at=sum(atarget)) %>% 
      mutate(ab.at=100*ab/at) %>% select(-ab, -at)
    aux <- left_join(aux1, aux2, by="fst")  
    result <- rbind(result, data.frame(scn=scn, aux))
  }
  return(result[-1,])
}


######################################## PCT LCT BURNT PER FST ########################################
pct.lct.burnt <- function(list.scn){
  
  load("rscripts/outs/obs.ab.lct.rdata")
  load("rscripts/outs/obs.ab.fst.lct.rdata")
  report <- data.frame(scn=NA, err=NA)
  report.fst <- data.frame(scn=NA, fst=NA, err=NA)
  for(scn in list.scn){
    # Diference with observed per fst
    ab.fst.lct.year <- read.table(paste0("outputs/", scn, "/_PctBurntLCT.FST.txt"), header=T) 
    ab.fst <- group_by(ab.fst.lct.year, scn, fst) %>% summarise(tot=sum(ab))
    ab.fst.lct <- group_by(ab.fst.lct.year, scn, fst, lct) %>% summarise(ab=sum(ab)) %>% 
      left_join(ab.fst, by=c("scn", "fst")) %>% mutate(pct=100*ab/tot)
    dif.fst <- left_join(ab.fst.lct, select(obs.ab.fst.lct, fst, lct, pct), by=c("fst", "lct")) %>% group_by(scn, fst) %>% 
      summarize(err = sqrt(sum((pct.x-pct.y)^2)))
    report.fst <- rbind(report.fst, as.data.frame(dif.fst))
    # Diference with observed 
    ab <- group_by(ab.fst.lct.year, scn) %>% summarise(tot=sum(ab))
    ab.lct <- group_by(ab.fst.lct.year, scn, lct) %>% summarise(ab=sum(ab)) %>% left_join(ab, by="scn") %>% 
      mutate(pct=100*ab/tot)
    dif <- left_join(ab.lct, select(obs.ab.lct, lct, pct), by="lct") %>% group_by(scn) %>% 
      summarize(err = sqrt(sum((pct.x-pct.y)^2)))
    report <- rbind(report, as.data.frame(dif))
  }
 
  report <- report[-1,]
  report.fst <- report.fst[-1,]
  return(list(report=report, report.fst=report.fst))
   
}

cook.pct.lct.burnt <- function(list.scn, nrun){
  ignitions <- read.table("inputfiles/FireIgnitions.txt", header=T) %>% select(fire.id, fst)
  lctype <- data.frame(spp=1:16,lct=c(rep("FC",4), rep("FD",3), "FC", "FD", "SH", rep("GC",3), rep("UR",3)))
  burnt.lct <- data.frame(fst=NA, lct=NA, ab=NA, year=NA, scn=NA)
  for(scn in list.scn){
    print(scn)
    for(r in 1:nrun){
      for(t in c(1,3:19,21:24)){
        load(paste0("outputs/", scn, "/Maps_r",r, "t", t, ".rdata"))
        if(t==1){
          load("inputlyrs/rdata/land.rdata")      
          land <- left_join(land, lctype, by="spp") %>% select(-tsdist, -spp)
        }
        if(t==12){
          LCF <- raster("inputlyrs/asc/ForestMapSpp00_31N-ETRS89.asc")
          land <- data.frame(cell.id=1:ncell(LCF), spp=LCF[]) %>% filter(!is.na(spp)) %>% 
            left_join( lctype, by="spp") %>% select(-spp)
        }
        land$fire.id <- map$id
        aux <- left_join(land, ignitions, by="fire.id") %>% group_by(fst, lct) %>% summarise(ab=length(lct)) %>% 
          filter(!is.na(fst)) %>% mutate(year=t, scn=scn)
        burnt.lct <- rbind(burnt.lct, as.data.frame(aux))
      }
    }
  }
  burnt.lct <- burnt.lct[-1,]
  
  # Burnt per FST
  ab.fst <- group_by(burnt.lct, scn, fst) %>% summarise(tot=sum(ab)) %>% filter(!is.na(fst))
  ab.fst.lct <- group_by(burnt.lct, scn, fst, lct) %>% summarise(ab=sum(ab)) %>% filter(!is.na(fst)) %>% 
                left_join(ab.fst, by=c("scn", "fst"))  %>% mutate(pct=round(100*ab/tot,1))
  # Total burnt
  ab <- group_by(burnt.lct, scn) %>% summarize(tot=sum(ab))
  ab.lct <- group_by(burnt.lct, scn, lct) %>% summarise(ab=sum(ab)) %>% filter(!is.na(lct)) %>% 
            left_join(ab, by="scn")  %>% mutate(pct=round(100*ab/tot,1))
    
  # Diference with observed per fst
  load("rscripts/outs/obs.ab.fst.lct.rdata")
  dif.fst <- left_join(ab.fst.lct, select(obs.ab.fst.lct, fst, lct, pct), by=c("fst", "lct")) %>% group_by(scn, fst) %>% 
          summarize(err = sqrt(sum((pct.x-pct.y)^2)))
  # Diference with observed 
  load("rscripts/outs/obs.ab.lct.rdata")
  dif <- left_join(ab.lct, select(obs.ab.lct, lct, pct), by="lct") %>% group_by(scn) %>% 
    summarize(err = sqrt(sum((pct.x-pct.y)^2)))
  
  
  return(list(ab.fst.lct=ab.fst.lct, ab.lct=ab.lct, dif.fst=dif.fst, dif=dif))
  
}


######################################## MATCH AREA BURNT ########################################
match.area <- function(list.scn){
  ## Info fires
  ignitions <- read.table("inputfiles/FireIgnitions.txt", header=T)
  ## Pctg of area match per fire and scenario
  report <- data.frame(scn=NA, year=NA, fire.id=NA, area=NA, fst=NA, pct=NA, am=NA)
  for(scn in list.scn){
    validation <- read.table(paste0("outputs/", scn, "/_AreaMatch.txt"), header=T)
    aux <- left_join(validation, ignitions, by=c("fire.id", "year")) %>%
           mutate(pctg=am/area*100) %>% group_by(scn, year, fire.id, area, fst) %>% 
          summarise(pct=round(mean(pctg),1), am=round(mean(am),1))
    report <- rbind(report, as.data.frame(aux))
  }
  report <- report[-1,]
  return(report)
}

best.scn <- function(report){
  ## Pctg of area match per scenario and fire.spread.type
  report.scn <- group_by(report, scn, fst) %>% summarise(min.match=min(pct), mn.match=mean(pct), 
                        md.match=median(pct), max.match=max(pct))
  best <- data.frame(scn=NA,	fst=NA,	min.match=NA,	mn.match=NA,	md.match=NA,	max.match=NA)
  ## Best scn per fire spread type
  for(type in 1:3){
    fires <- filter(report.scn, fst==type)
    md.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(median=max(md.match))
    mn.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(mean=max(mn.match))
    # filter(report.scn, fst==type, mn.match==mn.fst$mean)
    aux <- fires[order(fires$md.match, decreasing = T),]
    print(ifelse(type==1, "Wind", ifelse(type==2, "Topo", "Convectiu")))
    print(head(aux,5))
    best <- rbind(best, as.data.frame(aux))
  }
  return(best[-1,])
}

doubt.fires <- function(report){
  
  ## TOPOGRAPHIC FIRES
  type <- 2
  ## Let's find fires with very low rate of matching
  topo <- filter(report.fire, fst==type) %>% group_by(fire.id) %>% summarize(n=sum(pctg<5))
  topo[order(topo$n, decreasing = T),]
  ## Remove a few fires and repeat searching the best
  less.fires <- filter(report.fire, scn %notin% c(2229,689,26,184,4,266), fst==type)
  report.scn <- group_by(less.fires, scn, fst) %>% 
    summarise(min.match=min(pctg), mn.match=mean(pctg), 
              md.match=median(pctg), max.match=max(pctg))
  report.scn <- left_join(report.scn, all.scn, by="scn")
  mn.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(mean=max(mn.match))
  filter(report.scn, fst==type, mn.match==mn.fst$mean)
  
  
  ## CONVECTIVE FIRES
  type <- 3
  convec <- filter(report.fire, fst==type) %>% group_by(fire.id) %>% summarize(n=sum(pctg<5))
  convec[order(convec$n, decreasing = T),]
  less.fires <- filter(report.fire, scn %notin% c(213, 2146,2240,218,259), fst==type)
  report.scn <- group_by(less.fires, scn, fst) %>% 
    summarise(min.match=min(pctg), mn.match=mean(pctg), 
              md.match=median(pctg), max.match=max(pctg))
  report.scn <- left_join(report.scn, all.scn, by="scn")
  mn.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(mean=max(mn.match))
  filter(report.scn, fst==type, mn.match==mn.fst$mean)
  
  
  ## WIND FIRES
  type <- 1
  wind <- filter(report.fire, fst==type) %>% group_by(fire.id) %>% summarize(n=sum(pctg<5))
  wind[order(wind$n, decreasing = T),]
  less.fires <- filter(report.fire, scn %notin% c(2014, 2214, 2139, 2300, 2212), fst==type)
  report.scn <- group_by(less.fires, scn, fst) %>% 
    summarise(min.match=min(pctg), mn.match=mean(pctg), 
              md.match=median(pctg), max.match=max(pctg))
  report.scn <- left_join(report.scn, all.scn, by="scn")
  mn.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(mean=max(mn.match))
  filter(report.scn, fst==type, mn.match==mn.fst$mean)
  
  
  
  
  
  
}


######################################## PLOT ########################################
plot.fires <- function(scn, irun, t){
  # Mask
  if(t<=11)
    load("inputlyrs/rdata/mask.89-99.rdata")
  else  
    load("inputlyrs/rdata/mask.00-12.rdata")
  # Map dataframe  
  load(paste0("outputs/", scn, "/Maps_r", irun, "t", t, ".rdata"))
  # fire.ids
  MAP <- MASK
  MAP[!is.na(MAP[])] <- map$id
  writeRaster(MAP, paste0("outputs/", scn, "/FireIds_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
  ## fire.step
  MAP <- MASK
  MAP[!is.na(MAP[])] <- map$step
  writeRaster(MAP, paste0("outputs/", scn, "/FireStep_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
}  
