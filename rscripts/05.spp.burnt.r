############################### % BURNT / LCT / FST ############################### 
library(raster)
library(tidyverse)

rm(list=ls())
ignitions <- read.table("inputfiles/FireIgnitions.txt", header=T) %>% select(fire.id, fst)
lctype <- data.frame(spp=1:16,lct=c(rep("FC",4), rep("FD",3), "FC", "FD", "SH", rep("GC",3), rep("UR",3)))
burnt.lct <- data.frame(fst=NA, lct=NA, ab=NA, year=NA)

## Fires 1989 - 1999
load("inputlyrs/rdata/land.rdata")
land <- left_join(land, lctype, by="spp") %>% select(-tsdist, -spp)
load("inputlyrs/rdata/fireperim.89-99.rdata")
for(y in 2:ncol(perim.id)){
  land$fire.id <- perim.id[,y]
  aux <- left_join(land, ignitions, by="fire.id") %>% group_by(fst, lct) %>% summarise(ab=length(lct)) %>% 
    filter(!is.na(fst)) %>% mutate(year=names(perim.id)[y]) %>% filter(lct!="UR")
  burnt.lct <- rbind(burnt.lct, as.data.frame(aux))
}


## Fires 2000 - 2012
LCF <- raster("inputlyrs/asc/ForestMapSpp00_31N-ETRS89.asc")
land <- data.frame(cell.id=1:ncell(LCF), spp=LCF[]) %>% filter(!is.na(spp)) %>% 
        left_join( lctype, by="spp") %>% select(-spp)
load("inputlyrs/rdata/fireperim.00-12.rdata")
for(y in 2:ncol(perim.id)){
  land$fire.id <- perim.id[,y]
  aux <- left_join(land, ignitions, by="fire.id") %>% group_by(fst, lct) %>% summarise(ab=length(lct)) %>% 
    filter(!is.na(fst)) %>% mutate(year=names(perim.id)[y]) %>% filter(lct!="UR")
  burnt.lct <- rbind(burnt.lct, as.data.frame(aux))
}

## Totals
burnt.lct <- burnt.lct[-1,]
ab.fst <- group_by(burnt.lct, fst) %>% summarise(tot=sum(ab)) %>% filter(!is.na(fst))
ab.fst.lct <- group_by(burnt.lct, fst, lct) %>% summarise(ab=sum(ab)) %>% filter(!is.na(fst)) %>% 
  left_join(ab.fst, by="fst")  %>% mutate(pct=round(100*ab/tot))
ab.fst.lct
tot <- sum(burnt.lct$ab)
ab.lct <- group_by(burnt.lct, lct) %>% summarise(ab=sum(ab)) %>% filter(!is.na(lct)) %>% 
  mutate(pct=round(100*ab/tot))
ab.lct







############################### % BURNT / LCT / FIRE and then / FST ############################### 
library(raster)
library(tidyverse)

rm(list=ls())
fires <- read.table("C:/WORK/MEDMOD/SpatialModelsR/MEDSPREAD/inputfiles/FireIgnitions.txt", header=T)
fires.fst <- select(fires, fire.id, fst)
lctype <- data.frame(spp=1:16,lc=c(rep("FC",4), rep("FD",3), "FC", "FD", "SH", rep("GC",3), rep("OT",3)))
dta <- data.frame(year=NA, fire.id=NA, lc=NA, ab=NA, tot=NA, pct=NA, fst=NA)

## Fires 1989 - 1999
lctype <- data.frame(spp=1:16,lc=c(rep("FC",4), rep("FD",3), "FC", "FD", "SH", rep("GC",3), rep("OT",3)))
load("inputlyrs/rdata/land.rdata")
load("inputlyrs/rdata/fireperim.89-99.rdata")
perim.id <- pivot_longer(perim.id, y1989:y1999, names_to="year", values_to = "fire.id") %>%  filter(fire.id!=0)
for(y in sort(unique(perim.id$year))){
  one.year.fires <- filter(perim.id, year==y)
  ab.fire <- group_by(one.year.fires, fire.id) %>% summarise(tot=length(fire.id))
  aux <-  left_join(land, one.year.fires, by="cell.id") %>% left_join(lctype, by="spp") %>% 
    group_by(year, fire.id, lc) %>% summarise(ab=length(spp)) %>% left_join(ab.fire, by="fire.id") %>% 
    mutate(pct=round(100*ab/tot,1)) %>% left_join(fires.fst, by="fire.id") %>% filter(!is.na(fire.id))
  dta <- rbind(dta, as.data.frame(aux))
}

## Fires 2000 - 2012
lctype <- data.frame(spp=1:20,lc=c(rep("FC",7), rep("FD",6), "SH", rep("GC",3), rep("OT",3)))
LCF <- raster("inputlyrs/asc/ForestMapSpp00_31N-ETRS89.asc")
land <- data.frame(cell.id=1:ncell(LCF), spp=LCF[])
land <- land[!is.na(land$spp),]
load("inputlyrs/rdata/fireperim.00-12.rdata")
perim.id <- pivot_longer(perim.id, y2000:y2012, names_to="year", values_to = "fire.id") %>%  filter(fire.id!=0)
for(y in sort(unique(perim.id$year))){
  one.year.fires <- filter(perim.id, year==y)
  ab.fire <- group_by(one.year.fires, fire.id) %>% summarise(tot=length(fire.id))
  aux <-  left_join(land, one.year.fires, by="cell.id") %>% left_join(lctype, by="spp") %>% 
    group_by(year, fire.id, lc) %>% summarise(ab=length(spp)) %>% left_join(ab.fire, by="fire.id") %>% 
    mutate(pct=round(100*ab/tot,1)) %>% left_join(fires.fst, by="fire.id") %>% filter(!is.na(fire.id))
  dta <- rbind(dta, as.data.frame(aux))
}
dta <- dta[-1,]
dta.all <- filter(dta, lc!="OT")


## Plots pct burnt
p1 <- ggplot(dta.all, aes(pct)) +  geom_histogram(position="identity", bins=10) + facet_wrap("lc")
dta.wind <- filter(dta.all, fst==1)
p2 <- ggplot(dta.wind, aes(pct)) +  geom_histogram(position="identity", bins=10) + facet_wrap("lc")
dta.conv <- filter(dta.all, fst==3)
p3 <- ggplot(dta.conv, aes(pct)) +  geom_histogram(position="identity", bins=10) + facet_wrap("lc")
dta.topo <- filter(dta.all, fst==2)
p4 <- ggplot(dta.topo, aes(pct)) +  geom_histogram(position="identity", bins=10) + facet_wrap("lc")
jpeg("C:/WORK/MEDMOD/SpatialModelsR/MEDSPREAD/outputs/BurntLCType.tif", width = 700, height = 1000)
gridExtra::grid.arrange(p1,p2,p3,p4)
dev.off()

## Stats
group_by(dta.all, lc) %>% summarise(mn=mean(pct), sd=sd(pct))
group_by(dta.all, lc, fst) %>% summarise(mn=mean(pct), sd=sd(pct))
