library(raster)
library(viridis)

rm(list=ls())
scn.name <- "TestAll02"
fires <- read.table(paste0("outputs/", scn.name, "/Fires.txt"), header=T)
  # sum(fires$extra)/sum(fires$atarget)*100
group_by(fires, fst) %>% summarize(p=sum(extra)/sum(atarget)*100)

a <- group_by(fires, fire.id, fst) %>% summarize(p=mean(extra/atarget)*100)

FIRES <- raster(paste0("outputs/", scn.name, "/lyr/FireID_r1t6.tif"), header=T)
table(FIRES[])
plot(FIRES, col=viridis(12))


scn.name <- "TestSprd01"
fires <- read.table(paste0("outputs/", scn.name, "/Fires.txt"), header=T)
sprd <- read.table(paste0("outputs/", scn.name, "/FiresSprd.txt"), header=T)
fire.id <- 170
kk <- filter(sprd, fire.id==170)
