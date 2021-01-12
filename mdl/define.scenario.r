######################################################################################
###  define.scenario()
###
###  Description >  Initialize the scenario parameters and the global variables of 
###   MEDFIRE model
###
###  Arguments >  
###   scn.name : identificative name of the scenario (string)
###
###  Details > By default, the output directory is ..\outputs\scn.name\ and all the
###   objects are saved in the file scn.def.r.
###
###  Value >  An R script.
######################################################################################

define.scenario <- function(scn.name){

  cat("Initializing parameters", "\n")

  ## Output directory (do not never change that, please!)
  out.path <- paste0("outputs/", scn.name)
  
  ## Time lenght (in years) of a model simulation, from 1989 to 2012 both inclusive
  time.horizon <- 24
  
  ## Number of runs (i.e. replicas)
  nrun <- 10
  
  ## Flags to write spatial and tabular output data
  print.maps <- TRUE
  validation <- FALSE
  
  ## File of fire ignitions
  file.fire.ignis <- "FireIgnitions"
  
  ## Radius of the neighborhood (in pixels) to find out if a species is present in a region
  spp.distrib.rad <- 20 	# i.e. 2 km
  
  ## Fire parameters (should not change to much): Spread rate, burn probability, prescribed burns
  pb.upper.th <- 0.85
  pb.lower.th <- -1
  facc <- 1
  file.sprd.weight <- "WeightSprdFactors"
  rpb <- 0
  fire.intens.th <- 0.35
  fuel.opt <- "A"
  
  ## Save all the variables in .r file to be further loaded by landscape.dyn.r
  if(!file.exists(out.path))
    dir.create(file.path(getwd(), out.path), showWarnings = T) 
  dump(ls(), paste0(out.path, "/scn.def.r"))
  
}