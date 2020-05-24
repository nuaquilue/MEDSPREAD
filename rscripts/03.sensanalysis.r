rm(list=ls())
library(raster)
library(tidyverse)
id.scn <- c(paste0("00", 1:9), paste0("0", 10:99), 100:286)
rpb <- 0.9
nscn <- 205
result <- data.frame(scn=NA, run=NA, year=NA, fire.id=NA, am=NA)
for(i in 1:nscn){  
  scn.name <- paste0("Test", rpb*10, id.scn[i])
  for(r in 1:3){
    load("inputlyrs/rdata/mask.89-99.rdata")
    load("inputlyrs/rdata/fireperim.89-99.rdata")
    for(y in 1:11){
      cat(paste("Scn:", scn.name, "- run:", r, "- year:", y), "\n")
      perim.y <- perim.id[,c(1,y+1)]; names(perim.y)[2] <- "perim"
      FIRE <- raster(paste0("outputs/", scn.name, "/lyr/FireID_r", r, "t", y, ".tif"))
      dta <- data.frame(cell.id=1:ncell(MASK[]), mask=MASK[], fire.id=FIRE[])
      aux <- filter(dta, !is.na(mask)) %>% select(-mask) %>% left_join(perim.y, by="cell.id") %>%
             filter(!is.na(fire.id)) %>% mutate(dif=fire.id-perim) %>% filter(dif==0) %>% 
             group_by(fire.id) %>% summarize(am=length(fire.id))
      if(nrow(aux)>0)
        result <- rbind(result, data.frame(scn=scn.name, run=r, year=y, aux))    
    }
    load("inputlyrs/rdata/mask.00-12.rdata")
    load("inputlyrs/rdata/fireperim.00-12.rdata")
    for(y in 12:24){
      cat(paste("Scn:", scn.name, "- run:", r, "- year:", y), "\n")
      perim.y <- perim.id[,c(1,y-10)]; names(perim.y)[2] <- "perim"
      FIRE <- raster(paste0("outputs/", scn.name, "/lyr/FireID_r", r, "t", y, ".tif"))
      dta <- data.frame(cell.id=1:ncell(MASK[]), mask=MASK[], fire.id=FIRE[])
      aux <- filter(dta, !is.na(mask)) %>% select(-mask) %>% left_join(perim.y, by="cell.id") %>%
             filter(!is.na(fire.id)) %>% mutate(dif=fire.id-perim) %>% filter(dif==0) %>% 
             group_by(fire.id) %>% summarize(am=length(fire.id))
      if(nrow(aux)>0)
        result <- rbind(result, data.frame(scn=scn.name, run=r, year=y, aux))    
    }
  }
}
result <- result[-1,]


## Pctg of area match per fire and scenario
fire.ignis <- read.table("inputfiles/FireIgnitions.txt", header=T)
report.fire <- left_join(result, fire.ignis, by=c("fire.id", "year")) %>%
               mutate(pctg=am/area*100) %>% group_by(scn, year, fire.id, area, fst) %>% 
               summarise(pctg=round(mean(pctg),1), am=round(mean(am),1))
write.table(report.fire, paste0("rscripts/outs/ReportAreaMatchFire_0", rpb*10, ".txt"), quote=F, row.names=F, sep="\t")
            
## Pctg of area match per scenario and fire.spread.type
report.scn <- group_by(report.fire, scn, fst) %>% 
              summarise(min.match=min(pctg), mn.match=round(mean(pctg),1), max.match=max(pctg))
write.table(report.scn, paste0("rscripts/outs/ReportAreaMatchScn_0", rpb*10, ".txt"), quote=F, row.names=F, sep="\t")

## Best scn per fire spread type
scn.names <- paste0("Test", rpb*10, id.scn)
group_by(report.scn, fst) %>% summarise(best=which.max(mn.match)) %>%
    mutate(scn=scn.names[best])
