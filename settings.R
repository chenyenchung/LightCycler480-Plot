###################### Read Me ####################
# The sample name is expected to be sample name and replicate number
# Those after spliting character (set in the variable splitChr) will be read as the replicate number
# For example: This script will read "HGSN_1" as sample: "HGSN " and replicate 1
# When splitChr <- "_"
###################################################

intCtrl <- "Gapdh" # Internal ctrl gene
ctrlGroup <- "WTB" # Name of ctrl sample (Check readme for details)
ignore <- "ignore"
bioRep <- T # If there's multiple biological sample in this file, set as TRUE (NO QUOTE!)
logScale <- F # If you want to use log scale, set as TRUE. (Right now only for bioRep == TRUE)
fontSize <- 16 # Set font size for figure output
splitChr <- "_" #Setting the character spliting sample name and replicate number