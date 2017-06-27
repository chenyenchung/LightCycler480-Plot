###################### Read Me ####################
# The sample name is expected to be sample name and replicate number
# The last character will be read as  replicate number
# For example: This script will read "HGSN 1" as sample: "HGSN " and replicate 1
###################################################

intCtrl <- "Gapdh" # Internal ctrl gene
ctrlGroup <- "HGSN" # Name of ctrl sample (Check readme for details)
ignore <- "ignore"
bioRep <- F # If there's multiple biological sample in this file, set as TRUE (NO QUOTE!)
logScale <- F # If you want to use log scale, set as TRUE. (Right now only for bioRep == TRUE)
fontSize <- 18 # Set font size for figure output
splitChr <- " " #Setting the character spliting sample name and replicate number

#Load required packages. (Auto install if necessary)
dplyrEx <- require("dplyr")
if(!dplyrEx){
  install.pacakages("dplyr")
  library("dplyr")
}

ggplotEx <- require("ggplot2")
if(!ggplotEx){
  install.packages("ggplot2")
  library("ggplot2")
}

hmiscEx <- require("Hmisc")
if (!hmiscEx) {
  install.packages("Hmisc")
  library("Hmisc")
}


#Read raw file and primer layout

rawPath <- file.choose()
setwd(dirname(rawPath))
raw <- read.table(rawPath, header = TRUE, stringsAsFactors = FALSE, skip = 1, sep = "\t", fill = TRUE)
layout <- read.csv("primerlayout.csv", header = FALSE, stringsAsFactors = FALSE)
raw$gene <- as.character(unlist(as.data.frame(t(layout))))
workTbl <- tbl_df(raw[,c(4,5,9)])
workTbl$Cp[workTbl$Cp == 0] <- 40 
workTbl$Cp[is.na(workTbl$Cp)] <- 40 
workTbl <- workTbl[-which(workTbl$gene == ignore),]


if(!bioRep){
  sumTbl <- summarise(group_by(workTbl, Name, gene), mean = mean(Cp), sd = sd(Cp))
  
  #Check if Cp sd is high. If so, send a warning message.
  if (any(sumTbl$sd > 0.2)) {
    warnings("Some Cp values are not consistent between replicates.")
    write.table(sumTbl[which(sumTbl$sd > 0.2),], "CpSd.txt")
  }
  
  sumTbl$ref <- rep(filter(sumTbl, gene == intCtrl)$mean, each = nlevels(as.factor(sumTbl$gene)))
  sumTbl <- mutate(sumTbl, relCp = mean - ref, relExp = 2^-relCp)
  sumTbl$refExp <- rep(filter(sumTbl, Name == ctrlGroup)$relExp, nlevels(as.factor(sumTbl$Name)))
  sumTbl <- mutate(sumTbl, fc = relExp/refExp)
  
  #Plots for relative expression and fold change
  relBar <- ggplot(data = filter(sumTbl, !(gene == intCtrl)),
                  aes(x = gene, y = relExp, fill = Name))+
                  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "black") +
                  labs(x = "Gene", y = paste0("Relative Expression to ", intCtrl), fill = "Sample")+
                  theme_classic(base_size = fontSize)
  
  print(relBar)
  
  fcBar <- ggplot(data = filter(sumTbl, !(gene == intCtrl)),
                  aes(x = gene, y = fc, fill = Name))+
                  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "black") +
                  labs(x = "Gene", y = paste0("Fold Change to ", ctrlGroup), fill = "Sample")+
                  geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed")+
                  theme_classic(base_size = fontSize)
  
  print(fcBar)
                
  
}else{
  
#  workTbl$rep <- substr(workTbl$Name, nchar(workTbl$Name),nchar(workTbl$Name))
#  workTbl$Name <- substr(workTbl$Name, 1, nchar(workTbl$Name) - 1)
  workTbl$rep <- sapply(strsplit(workTbl$Name, splitChr),"[",2)
  workTbl$Name <- sapply(strsplit(workTbl$Name, splitChr),"[",1)
  sumTbl <- summarise(group_by(workTbl, rep, Name, gene), mean = mean(Cp), sd = sd(Cp))
  
  #Check if Cp sd is high. If so, send a warning message.
  if (any(sumTbl$sd > 0.2)) {
    warnings("Some Cp values are not consistent between replicates.")
    write.table(sumTbl[which(sumTbl$sd > 0.2),], "CpSd.txt")
  }
  
  sumTbl$ref <- rep(filter(sumTbl, gene == intCtrl)$mean, each = nlevels(as.factor(sumTbl$gene)))
  sumTbl <- mutate(sumTbl, relCp = mean - ref, relExp = 2^-relCp)
  refTbl <- filter(sumTbl, Name == ctrlGroup)
  expTbl <- filter(sumTbl, !(Name == ctrlGroup))
  expTbl$fc <- expTbl$relExp/refTbl$relExp
  
  #Plots for relative expression and fold change
  relPlot <- ggplot(data = filter(sumTbl, !(gene == intCtrl)), aes(x = gene, y = relExp, colour = Name, group = Name)) +
                    geom_jitter(position = position_dodge(width = 0.5), size = 3, alpha = 0.3)+
                    stat_summary(fun.y = mean,
                                 geom = "point",
                                 pch = "\U2013",
                                 position = position_dodge(width = 0.5),
                                 size = 8,
                                 show.legend = FALSE)+
                    stat_summary(fun.data = mean_sdl,
                                 geom = "errorbar",
                                 position = position_dodge(width = 0.5),
                                 fun.args = list(mult = 1),
                                 width = 0.5,
                                 alpha = 0.5,
                                 show.legend = FALSE)+
                    labs(x = "Gene", y = paste0("Relative Expression to ", intCtrl), colour = "Sample")+
                    theme_classic(base_size = fontSize)
  
  if (logScale) {print(relPlot + scale_y_continuous(trans = "log10"))}else{print(relPlot)}
  
  fcplot <- ggplot(data = filter(expTbl, !(gene == intCtrl)),
                   aes(x = gene, y = fc, colour = Name, group = Name))+
    geom_jitter(position = position_dodge(width = 0.5), size = 3, alpha = 0.3)+
    stat_summary(fun.y = mean,
                 geom = "point",
                 pch = "\U2013",
                 position = position_dodge(width = 0.5),
                 size = 8,
                 show.legend = FALSE)+
    stat_summary(fun.data = mean_sdl,
                 geom = "errorbar",
                 position = position_dodge(width = 0.5),
                 fun.args = list(mult = 1),
                 width = 0.5,
                 alpha = 0.5,
                 show.legend = FALSE)+
    labs(x = "Gene", y = paste0("Fold Change to ", ctrlGroup), colour = "Sample")+
    geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed")+
    theme_classic(base_size = fontSize)
  
  if(logScale){print(fcplot + scale_y_continuous(trans = "log10"))}else{print(fcplot)}
}