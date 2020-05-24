############################################ RUN A SCN ##################################################
rm(list=ls())
# Load functions
source("mdl/define.scenario.r")
source("mdl/land.dyn.mdl.r")  
scn.name <- "Test001"
define.scenario(scn.name)
# Change target parameters
nrun <- 10
rpb <- 0.4
file.sprd.weight <- "WeightSprdFactors"
# Write the name of the customized parameters in the dump function. 
# It copies these R objects into the file outputs/test/scn.custom.def.r
dump(c("nrun", "rpb", "file.sprd.weight"), paste0("outputs/", scn.name, "/scn.custom.def.r"))
# Run the model
land.dyn.mdl(scn.name)


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
# Load functions
source("mdl/define.scenario.r")
source("mdl/land.dyn.mdl.r")
wfactors <- read.table("scenarios/wfactors.txt", header=T)
id.scn <- c(paste0("00", 1:9), paste0("0", 10:99), 100:286)
nrun <- 3; rpb <- 0.9
for(i in 206:232){
  scn.name <- paste0("Test", rpb*10, id.scn[i])
  ## Change weights of spread factors
  x <- unlist(filter(wfactors, scn==i) %>% select(-scn))
  sprdw <- data.frame(factor=c("wind", "slope", "flam", "aspc"), fst.w=x, fst.c=x, fst.t=x)
  write.table(sprdw, paste0("inputfiles/WeightSprdFactors_", rpb*10, id.scn[i], ".txt"), quote=F, sep="\t", row.names=F)
  define.scenario(scn.name)
  file.sprd.weight <- paste0("WeightSprdFactors_", rpb*10, id.scn[i])
  dump(c("nrun", "rpb", "file.sprd.weight"), paste0("outputs/", scn.name, "/scn.custom.def.r"))
  land.dyn.mdl(scn.name)
}



