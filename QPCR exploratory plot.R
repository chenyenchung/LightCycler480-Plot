intCtrl <- "Gapdh" # Internal ctrl gene
ctrlGroup <- "1-1 " # Name of ctrl sample (Check readme for details)
bioRep <- TRUE # If there's multiple biological sample in this file, set as TRUE (NO QUOTE!)

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


rawPath <- file.choose()
setwd(dirname(rawPath))
raw <- read.table(rawPath, header = TRUE, stringsAsFactors = FALSE, skip = 1, sep = "\t")
layout <- read.csv("primerlayout.csv", header = FALSE, stringsAsFactors = FALSE)
raw$gene <- as.character(unlist(as.data.frame(t(layout))))
workTbl <- tbl_df(raw[,c(4,5,9)])

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
  
  print(ggplot(data = filter(sumTbl, !(gene == intCtrl)),
               aes(x = gene, y = relExp, fill = Name))+
          geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
          labs(x = "Gene", y = paste0("Relative Expression to ", intCtrl), fill = "Sample")+
          theme_classic(base_size = 15))
  
  print(ggplot(data = filter(sumTbl, !(gene == intCtrl)),
               aes(x = gene, y = fc, fill = Name))+
          geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
          labs(x = "Gene", y = paste0("Fold Change to ", ctrlGroup), fill = "Sample")+
          geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed")+
          theme_classic(base_size = 15))
}else{
  workTbl$rep <- substr(workTbl$Name, nchar(workTbl$Name),nchar(workTbl$Name))
  workTbl$Name <- substr(workTbl$Name, 1, nchar(workTbl$Name) - 1)
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
  
  print(ggplot(data = filter(sumTbl, !(gene == intCtrl)), aes(x = gene, y = relExp, colour = Name, group = Name)) +
    geom_jitter(position = position_dodge(width = 0.5), size = 3)+
    ylim(0, max(sumTbl$relExp)*1.2)+
    theme_classic(base_size = 15))
  
  print(ggplot(data = filter(expTbl, !(gene == intCtrl)),
         aes(x = gene, y = fc, colour = Name, group = Name))+
          geom_jitter(position = position_dodge(width = 0.5), size = 3, alpha = 0.3)+
          stat_summary(fun.data = mean_sdl,
                       position = position_dodge(width = 0.5),
                       fun.args = list(mult = 1))+
          labs(x = "Gene", y = paste0("Fold Change to ", ctrlGroup), colour = "Sample")+
          geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed")+
          theme_classic(base_size = 15))
}
