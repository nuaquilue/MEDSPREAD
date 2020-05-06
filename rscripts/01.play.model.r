rm(list=ls())
 setwd("c:/work/MEDMOD/SpatialModelsR/MEDFIRE")  #Nú HP
# setwd("d:/MEDMOD/SpatialModelsR/MEDFIRE")   #CTFC

############################################ RUN A SCN ##################################################
# Load functions
source("mdl/define.scenario.r")
source("mdl/land.dyn.mdl.r")  
# Name the new scenario and call the define.scenario function to load default initialization of model’s parameters
scn.name <- "TestPerformanceFire"
define.scenario(scn.name)
# Change target parameters
file.clim.severity <- "ClimaticSeverity_test"
rpb <- 0.4
nrun <- 1
write.sp.outputs <- T
time.horizon <- 10
processes <- c(TRUE,   # 1. Climate change
               FALSE,  # 2. Land-cover changes
               FALSE,  # 3. Forest management
               TRUE,   # 4. Wildfires
               FALSE,   # 5. Prescribed burns
               F,  # 6. Drought
               TRUE,   # 7. Post-fire regeneration
               F,  # 8. Cohort establihsment
               F,  # 9. Afforestation
               TRUE)   # 10. Growth
# Write the name of the customized parameters in the dump function. 
# It copies these R objects into the file outputs/test/scn.custom.def.r
dump(c( "file.clim.severity", "rpb", "nrun", "write.sp.outputs", "time.horizon", "processes"), 
     paste0("outputs/", scn.name, "/scn.custom.def.r"))
# Run the model
system.time(land.dyn.mdl(scn.name))


###########################################################################################################
## Create .Rdata with static variables of the model, only run once for all scenarios!
work.path <- getwd()
source("mdl/read.static.vars.r")
read.static.vars(work.path)
## Create .Rdata with initial values of variables of the model, used at each replicate of any scn.
source("mdl/read.state.vars.r")
read.state.vars(work.path)
## Create 2 data frames per climatic scn and decade: SDMs of all spp, and climatic variables (temp & precip) for CAT
source("mdl/read.climatic.vars.r")
source("mdl/read.sdm.r")
read.climatic.vars(work.path)
read.sdm(work.path, "base")
## Save interfaces
source("mdl/update.interface.r")
load("inputlyrs/rdata/land.rdata")
interface <- update.interface(land)
save(interface, file="inputlyrs/rdata/interface.rdata")
# write as raster
load("inputlyrs/rdata/mask.rdata")
INTERFACE <- MASK
INTERFACE[!is.na(MASK[])] <- interface$x
plot(INTERFACE, col=rainbow(7))
writeRaster(INTERFACE, "C:/WORK/MEDMOD/SpatialModelsR/MEDFIRE/inputlyrs/asc/Interface_100m_31N-ETRS89.asc",
            format="ascii", NAflag=-1, overwrite=T)



