rm(list=ls())
library(raster)
library(tidyverse)
load("inputlyrs/rdata/mask.rdata")
load("inputlyrs/rdata/fireperim.89-99.rdata")
id.scn <- c(paste0("00", 1:9), paste0("0", 10:99), 100:286)
rpb <- 0.1
i <- 1
result <- data.frame(scn=NA, run=NA, year=NA, fire.id=NA, am=NA)
for(i in 1:3){  #286
  scn.name <- paste0("Test", rpb*10, id.scn[i])
  for(r in 1:3){
    for(y in 1:11){
      cat(paste)
      perim.y <- perim.id[,c(1,y+1)]; names(perim.y)[2] <- "perim"
      FIRE <- raster(paste0("outputs/", scn.name, "/lyr/FireID_r", r, "t", y, ".tif"))
      dta <- data.frame(cell.id=1:ncell(MASK[]), mask=MASK[], fire.id=FIRE[])
      dta$fire.id[is.na(dta$fire.id)] <- 0
      aux <- filter(dta, !is.na(mask)) %>% select(-mask) %>% left_join(perim.y, by="cell.id") %>%
             filter(fire.id!=0 & perim!=0) %>% group_by(fire.id) %>% summarize(am=length(fire.id))
      if(nrow(aux)>0)
        result <- rbind(result, data.frame(scn=scn.name, run=r, year=y, aux))    
    }
  }
}

## Pctg of area match per fire and scenario
fire.ignis <- read.table("inputfiles/FireIgnitions.txt", header=T)
result <- result[-1,]
report.fire <- left_join(result, fire.ignis, by=c("fire.id", "year")) %>%
               mutate(pctg=am/area*100) %>% group_by(scn, year, fire.id, area, fst, wind) %>% 
               summarise(pctg=round(mean(pctg),1))

## Pctg of area match per scenario and fire.spread.type
report.scn <- group_by(report.fire, scn, fst) %>% 
              summarise(min.match=min(pctg), mn.match=mean(pctg), max.match=max(pctg))

## Best scn per fire spread type
scn.names <- paste0("Test", rpb*10, id.scn)
group_by(report.scn, fst) %>% summarise(best=which.max(mn.match)) %>%
  mutate(scn=scn.names[best])
