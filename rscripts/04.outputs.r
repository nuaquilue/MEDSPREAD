library(raster)
library(viridis)
library(tidyverse)
rm(list=ls())
scn.name <- "TestAll16"
fires <- read.table(paste0("outputs/", scn.name, "/Fires.txt"), header=T)
fires$pextra <- round(fires$extra/fires$atarget*100,1)
group_by(fires, fst) %>% summarize(p=round(sum(extra)/sum(atarget)*100,1))


## in medfire
scn.name <- "Scn_FireShape01sp"
fires <- read.table(paste0("C:/WORK/MEDMOD/SpatialModelsR/MEDFIRE/outputs/", scn.name, "/Fires.txt"), header=T)
fires$prem <- round(fires$rem/fires$atarget*100,1)
group_by(fires, fst) %>% summarize(p=round(sum(rem)/sum(atarget)*100,1))


## nfires per year
nfires <- group_by(fires, year) %>% summarize(n=length(year), a=sum(atarget)) %>% mutate(y=year+1988)

FIRES <- raster(paste0("outputs/", scn.name, "/lyr/FireID_r1t6.tif"), header=T)
table(FIRES[])
plot(FIRES, col=viridis(12))

fires <- read.table(paste0("outputs/", scn.name, "/Fires.txt"), header=T)
sprd <- read.table(paste0("outputs/", scn.name, "/FiresSprd.txt"), header=T)
id <- 191
kk <- filter(sprd, fire.id==id)

xungos <- filter(fires, extra>200)
ids <- xungos$fire.id
for(id in ids){
  kk <- filter(sprd, fire.id==id)
  print(paste("fire", id, "- area", unlist(filter(fires, fire.id==id) %>% dplyr::select(atarget)),
              "- burnt", unlist(filter(fires, fire.id==id) %>% dplyr::select(aburnt))               ))
  print(table(kk$step))
}
  
## some numbers
round(100*sum(fires$extra>0)/nrow(fires),1)
round(100*sum(fires$pextra>10)/nrow(fires),1)

round(100*sum(fires$pextra>10 & fires$atarget>2000)/sum(fires$atarget>2000),1)
round(100*sum(fires$pextra>10 & fires$atarget>4000)/sum(fires$atarget>4000),1)
round(100*sum(fires$pextra>10 & fires$atarget>10000)/sum(fires$atarget>10000),1)
