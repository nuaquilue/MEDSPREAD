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
file.sprd.weight <- "SprdRateWeights"
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

