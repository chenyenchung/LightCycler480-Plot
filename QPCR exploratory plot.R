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

# Append primerlayout on rawdata
primerArr <- function(file) {
  raw <- read.table(paste0(file, ".txt"), header = TRUE, stringsAsFactors = FALSE, skip = 1, sep = "\t", fill = TRUE)
  layout <- read.csv(paste0(file, ".csv"), header = FALSE, stringsAsFactors = FALSE)
  raw$gene <- as.character(unlist(as.data.frame(t(layout))))
  return(raw)
}

#Read raw file and primer layout
rawPath <- file.choose()
setwd(dirname(rawPath))
source("settings.R")
fllist <- list.files(path = dirname(rawPath))[grep(".txt", list.files(path = dirname(rawPath)))]
flname <- gsub(".txt", "", fllist)
flname <- flname[!flname == "CpSd"]
raw <- do.call("rbind", lapply(flname, primerArr))
workTbl <- tbl_df(raw[,c(4,5,9)])
workTbl$Cp[workTbl$Cp == 0] <- 40 
workTbl$Cp[is.na(workTbl$Cp)] <- 40 
workTbl <- workTbl[!(workTbl$gene == ignore),]


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
  
####  Order the sample and genes here (No biorep) ####  
  #sumTbl$Name <- factor(sumTbl$Name, levels = c("H3Gold", "H3Gnew", "KO14", "KO22"))
  #sumTbl$gene <- factor(sumTbl$gene, levels = c("Meg3v1", "Meg3v5","Ezh2"))#"Hb9", "Hoxa5","Hoxc8"))
###
  
  #Plots for relative expression and fold change
  relBar <- ggplot(data = filter(sumTbl, !(gene == intCtrl)),
                  aes(x = gene, y = relExp, fill = Name))+
                  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "black") +
                  labs(x = "Gene", y = paste0("Relative Expression to ", intCtrl), fill = "Sample")+
                  theme_classic(base_size = fontSize) #+ scale_fill_manual(values = c("blue","lightblue","blue","lightblue","red","pink"))
  
  print(relBar)
  
  fcBar <- ggplot(data = filter(sumTbl, !(gene == intCtrl)),
                  aes(x = gene, y = fc, fill = Name))+
                  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "black") +
                  labs(x = "Gene", y = paste0("Fold Change to ", ctrlGroup), fill = "Sample")+
                  geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed")+
                  theme_classic(base_size = fontSize) #+ scale_fill_manual(values = c("blue","lightblue","blue","lightblue","red","pink"))
  
  
  print(fcBar)
                
  
}else{
  

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

####  Order the sample and genes here (Biorep) ####
  
#  sumTbl$Name <- factor(sumTbl$Name, levels = c("WTB","IgDMR", "V1OE1","V1OE7", "V1OE8","V1OE15"))
#  sumTbl$gene <- factor(sumTbl$gene, levels = c("Meg3v1", "Meg3v5","Ezh2"))#"Hb9", "Hoxa5","Hoxc8"))
  
####
    
  refTbl <- filter(sumTbl, Name == ctrlGroup)
  expTbl <- arrange(filter(sumTbl, !(Name == ctrlGroup)), Name, rep)
  expTbl$fc <- expTbl$relExp/refTbl$relExp
  refTbl$fc <- 1
  fcTbl <- rbind(expTbl, refTbl)
  
  #Plots for relative expression and fold change
  relplot <- ggplot(data = filter(sumTbl, !(gene == intCtrl)), aes(x = gene, y = relExp, colour = Name, group = Name)) +
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
                    theme_classic(base_size = fontSize)+
                    scale_color_manual(values = c("black","red",rep("green",4)))
                    
  
  if (logScale) {print(relplot + scale_y_continuous(trans = "log10"))}else{print(relplot)}
  
  fcplot <- ggplot(data = filter(fcTbl, !(gene == intCtrl)),
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
    theme_classic(base_size = fontSize)+
    scale_color_manual(values = c("black","red",rep("green",4)))
  
  if(logScale){print(fcplot + scale_y_continuous(trans = "log10"))}else{print(fcplot)}
}
