rm(list=ls())
library(readxl)
library(raster)
library(tidyverse)
id.scn <- c(paste0("00", 1:9), paste0("0", 10:99), 100:286)  
rpb <- 0.1
nscn <- 286
result <- data.frame(scn=NA, run=NA, year=NA, fire.id=NA, am=NA)
for(i in 1:nscn){  
  scn.name <- paste0("Test", rpb*10, id.scn[i])
  for(r in 1:3){
    load("inputlyrs/rdata/mask.89-99.rdata")
    load("inputlyrs/rdata/fireperim.89-99.rdata")
    for(y in c(1,3:11)){
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
    for(y in c(12:19,21:24)){
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
save(result, file="rscripts/outs/result_02_set1.rdata")

## Pctg of area match per fire and scenario
fire.ignis <- read.table("inputfiles/FireIgnitions.txt", header=T)
report.fire <- left_join(result, fire.ignis, by=c("fire.id", "year")) %>%
               mutate(pctg=am/area*100) %>% group_by(scn, year, fire.id, area, fst) %>% 
               summarise(pctg=round(mean(pctg),1), am=round(mean(am),1))
write.table(report.fire, paste0("rscripts/outs/ReportAreaMatchFire_0", rpb*10, "_set1.txt"), quote=F, row.names=F, sep="\t")



## SEE ALL SCN TOGETHER
rm(list=ls())
library(readxl)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
rpb <- 0.1
report.fire <- read.table(paste0("rscripts/outs/ReportAreaMatchFire_0", rpb*10, ".txt"), header=T)
for(rpb in seq(0.2,0.3,0.1)){
  a <- read.table(paste0("rscripts/outs/ReportAreaMatchFire_0", rpb*10, ".txt"), header=T)
  print(nrow(a))
  report.fire <- rbind(report.fire,a)
}
            

## Pctg of area match per scenario and fire.spread.type
report.scn <- group_by(report.fire, scn, fst) %>% 
              summarise(min.match=min(pctg), mn.match=mean(pctg), 
                        md.match=median(pctg), max.match=max(pctg))
# write.table(report.scn, "rscripts/outs/ReportAreaMatchScn.txt", quote=F, row.names=F, sep="\t")

## Read scenario parameters and join
scn.param <- read_xlsx("d:/MEDMOD/SpatialModelsR/MEDSPREAD/scenarios/WeightFactors.xlsx", sheet="Hoja1") %>% select(-tot)
scn.param$n <- scn.param$scn
all.scn <- scn.param %>% select(-n)
all.scn$scn <- paste0("Test", 1000+scn.param$n)
aux <- scn.param %>% select(-n)
aux$scn <- paste0("Test", 2000+scn.param$n)
all.scn <- rbind(all.scn, aux)
aux <- scn.param %>% select(-n)
aux$scn <- paste0("Test", 3000+scn.param$n)
all.scn <- rbind(all.scn, aux)
report.scn <- left_join(report.scn, all.scn, by="scn")
rm(aux)

## Best scn per fire spread type
# wind
type <- 1
fires <- filter(report.scn, fst==type)
md.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(median=max(md.match))
mn.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(mean=max(mn.match))
filter(report.scn, fst==type, mn.match==mn.fst$mean)
fires[order(fires$mn.match, decreasing = T),]
# topo
type <- 2
fires <- filter(report.scn, fst==type)
md.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(median=max(md.match))
mn.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(mean=max(mn.match))
filter(report.scn, fst==type, mn.match==mn.fst$mean)
# conv
type <- 3
fires <- filter(report.scn, fst==type)
md.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(median=max(md.match))
mn.fst <- filter(report.scn, fst==type) %>% group_by(fst) %>% summarise(mean=max(mn.match))
filter(report.scn, fst==type, mn.match==mn.fst$mean)


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





