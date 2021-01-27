library(raster)
library(viridis)
library(tidyverse)
rm(list=ls())

scn.name <- "Scn_pbEXPFI_windAslopeAfuelB_acc2_rpb3_up8_wDefault"
fires <- read.table(paste0("outputs/", scn.name, "/_Fires.txt"), header=T)
fires$pextra <- round(fires$extra/fires$atarget*100,1)
group_by(fires, fst) %>% summarize(p=round(sum(extra)/sum(atarget)*100,1))

##
FIRES <- raster(paste0("outputs/", scn.name, "/FireIds_r1t6.tif"), header=T)
table(FIRES[])
plot(FIRES, col=viridis(12))

fires <- read.table(paste0("outputs/", scn.name, "/_Fires.txt"), header=T)
sprd <- read.table(paste0("outputs/", scn.name, "/_Spread.txt"), header=T)
id <- 191
kk <- filter(sprd, fire.id==id)

xungos <- filter(fires, extra>200)
ids <- xungos$fire.id
for(id in ids){
  kk <- filter(sprd, fire.id==id)
  print(paste("fire", id, "- area", unlist(filter(fires, fire.id==id) %>% dplyr::select(atarget)),
              "- burnt", unlist(filter(fires, fire.id==id) %>% dplyr::select(aburnt.highintens, aburnt.lowintens))               ))
  print(table(kk$step))
}
  ## some numbers
round(100*sum(fires$extra>0)/nrow(fires),1)
round(100*sum(fires$pextra>10)/nrow(fires),1)
round(100*sum(fires$pextra>10 & fires$atarget>2000)/sum(fires$atarget>2000),1)
round(100*sum(fires$pextra>10 & fires$atarget>4000)/sum(fires$atarget>4000),1)
round(100*sum(fires$pextra>10 & fires$atarget>10000)/sum(fires$atarget>10000),1)


## PCT LCT BURNT PER FIRE and FST
rm(list=ls())
scn.name <- "TestMdl"; r <- 1
ignitions <- read.table("inputfiles/FireIgnitions.txt", header=T) %>% select(fire.id, fst)
load("inputlyrs/rdata/land.rdata")
lctype <- data.frame(spp=1:16,lct=c(rep("FC",4), rep("FD",3), "FC", "FD", "SH", rep("GC",3), rep("UR",3)))
land <- left_join(land, lctype, by="spp")
burnt.lct.fire <- data.frame(id=NA,  lct=NA, ab=NA, tot=NA, fst=NA, pct=NA)
for(t in c(1,3:6)){
  load(paste0("C:/WORK/MEDMOD/SpatialModelsR/MEDSPREAD/outputs/", scn.name, "/Maps_r",r, "t", t, ".rdata"))
  ab.fire <- group_by(map, id) %>% summarise(tot=length(id)) %>% filter(!is.na(id)) %>%
            left_join(ignitions, by=c("id"="fire.id"))
  aux <-  left_join(land, map, by="cell.id") %>% group_by(id, lct) %>% summarise(ab=length(spp)) %>% 
          left_join(ab.fire, by="id") %>% mutate(pct=100*ab/tot) %>% filter(!is.na(id)) %>% filter(lct!="UR")
  burnt.lct.fire <- rbind(burnt.lct.fire, as.data.frame(aux))
}
burnt.lct.fire <- burnt.lct.fire[-1,]
print(group_by(burnt.lct.fire, fst, lct) %>% summarise(mn=round(mean(pct)), sd=round(sd(pct))))
print(group_by(burnt.lct.fire, lct) %>% summarise(mn=round(mean(pct)), sd=round(sd(pct))))


