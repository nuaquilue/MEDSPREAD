library(readxl)
library(raster)
library(viridis)
library(tidyverse)

############################################ RUN SOME SCN ##################################################
rm(list=ls())
source("mdl/define.scenario.r"); source("mdl/land.dyn.mdl.r") 
scenarios <- read_xlsx("Scenarios.xlsx", sheet="Obj5")
for(i in 29:33){
  scn.name <- scenarios$scn.name[i]
  define.scenario(scn.name)
  nrun <- 5
  print.maps <- F
  validation <- T
  fuel.opt <- scenarios$fuel.opt[i]
  rpb <- scenarios$rpb[i]
  pb.upper.th <- scenarios$pb.upper.th[i]
  facc <- scenarios$facc[i]
  wwind <- scenarios$wwind[i]
  wslope <- scenarios$wslope[i]
  dump(c("nrun", "fuel.opt", "rpb", "pb.upper.th", "facc", "wwind", "wslope", "validation", "print.maps"), 
       paste0("outputs/", scn.name, "/scn.custom.def.r"))
  land.dyn.mdl(scn.name)
}



############################################## INITIALIZATION ########################################################
## Create .Rdata with static variables of the model, only run once for all scenarios!
## Create .Rdata with initial values of variables of the model, used at each replicate of any scn.
rm(list=ls())
source("mdl/read.static.vars.r")
source("mdl/read.state.vars.r")
read.state.vars(year=89)
read.static.vars()


############################################ SENSITIVITY ANALYSIS ##################################################
rm(list=ls())
library(raster)
library(tidyverse)
select <- dplyr::select  
# Load functions
source("mdl/define.scenario.r")
source("mdl/land.dyn.mdl.r")
wfactors <- read.table("inputfiles/wfactors.txt", header=T)
id.scn <- c(paste0("00", 1:9), paste0("0", 10:99), 100:286)
nrun <- 3; fi.accelerate <- 5; rpb <- 0.4
write.sp.outputs <- F; validation <- T
for(i in 1:286){
  scn.name <- paste0("Test", rpb*10, id.scn[i])
  ## Change weights of spread factors
  x <- unlist(filter(wfactors, scn==i) %>% dplyr::select(-scn))
  sprdw <- data.frame(factor=c("r", "wind", "slope", "flam", "aspc"), fst.w=c(rpb, x), fst.t=c(rpb, x), fst.c=c(rpb, x))
  write.table(sprdw, paste0("inputfiles/WeightSprdFactors_", rpb*10, id.scn[i], ".txt"), quote=F, sep="\t", row.names=F)
  define.scenario(scn.name)
  file.sprd.weight <- paste0("WeightSprdFactors_", rpb*10, id.scn[i])
  dump(c("nrun", "fi.accelerate", "file.sprd.weight", "write.sp.outputs", "validation"), 
       paste0("outputs/", scn.name, "/scn.custom.def.r"))
  land.dyn.mdl(scn.name)
}



