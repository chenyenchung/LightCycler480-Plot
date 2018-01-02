# Load required packages. (Auto install if necessary)
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

cowEx <- require("cowplot")
if (!cowEx) {
  install.packages("cowplot")
  library("cowplot")
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
#  sumTbl$Name <- factor(sumTbl$Name, levels = c("H3G@ES","H3G@D2","H3G@D4", "H3G@D7"))
#  sumTbl$gene <- factor(sumTbl$gene, levels = c("Gapdh","Oct4","Sox1","Pax6","Hba-a1","Hbb-bh1"))
###
  
  #Plots for relative expression and fold change
  relBar <- ggplot(data = filter(sumTbl, !(gene == intCtrl)),
                  aes(x = gene, y = relExp, fill = Name))+
                  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "black") +
                  labs(x = "Gene", y = paste0("Relative Expression to ", intCtrl), fill = "Sample") +
                  #geom_text(aes(label = Name),
                  #          position = position_dodge(width = 0.8), size = 3,
                  #          angle = -90, hjust = -0.5) +
                  theme_classic(base_size = fontSize)# +
                  #scale_fill_manual(values = c("Dodgerblue4","Dodgerblue3","Dodgerblue2","Dodgerblue1"))#+
                  #guides(fill = FALSE)
  
  if(logScale){print(relBar + scale_y_continuous(trans = "log10"))}else{print(relBar)}
  
  fcBar <- ggplot(data = filter(sumTbl, !(gene == intCtrl)),
                  aes(x = gene, y = fc, fill = Name))+
                  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "black") +
                  labs(x = "Gene", y = paste0("Fold Change to ", ctrlGroup), fill = "Sample")+
                  geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed")+
                  #geom_text(aes(label = Name),
                  #position = position_dodge(width = 0.8), size = 3,
                  #angle = -90, hjust = 0.5) +
                  theme_classic(base_size = fontSize)# +
                  #scale_fill_manual(values = c("Dodgerblue4","Dodgerblue3","Dodgerblue2","Dodgerblue1")) #+
                  #guides(fill = FALSE)
  
  
  if(logScale){print(fcBar + scale_y_continuous(trans = "log10"))}else{print(fcBar)}
                
  
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
  
#  sumTbl$Name <- factor(sumTbl$Name, levels = c("WTB","IgDMR","V1OE1","V1OE7","V1OE8","V1OE15"))
# #                                               "V5OE13","V5OE24","DOE13"))
# # levels(sumTbl$Name)[levels(sumTbl$Name) == "V1OE15"] <- "Meg3v1 OE"
#  levels(sumTbl$Name)[levels(sumTbl$Name) == "WTB"] <- "Wildtype Ctrl (p10-13)"
#  levels(sumTbl$Name)[levels(sumTbl$Name) == "IgDMR"] <- "IgDMR MatKO (p12-15)"
#  levels(sumTbl$Name)[levels(sumTbl$Name) == "V1OE1"] <- "Meg3v1 OE #1 (p11-14)"
#  levels(sumTbl$Name)[levels(sumTbl$Name) == "V1OE7"] <- "Meg3v1 OE #7 (p11-14)"
#  levels(sumTbl$Name)[levels(sumTbl$Name) == "V1OE8"] <- "Meg3v1 OE #8 (p11-14)"
#  levels(sumTbl$Name)[levels(sumTbl$Name) == "V1OE15"] <- "Meg3v1 OE #15 (p11-14)"
# # levels(sumTbl$Name)[levels(sumTbl$Name) == "V5OE13"] <- "Meg3v5 OE #13 (p26-28)"
# # levels(sumTbl$Name)[levels(sumTbl$Name) == "V5OE24"] <- "Meg3v5 OE #24 (p26-28)"
#  sumTbl$gene <- factor(sumTbl$gene, levels = c("Meg3v1","Meg3v5","Mirg","Rian"))#,"Jarid2", "Dlk1","Dio3","Rtl1","Mirg","Rian"))
#   
####
    
  refTbl <- filter(sumTbl, Name == ctrlGroup)
  expTbl <- arrange(filter(sumTbl, !(Name == ctrlGroup)), Name, rep)
  expTbl$fc <- expTbl$relExp/refTbl$relExp
  refTbl$fc <- 1
  fcTbl <- rbind(expTbl, refTbl)
  
  #Plots for relative expression and fold change
  relplot <- ggplot(data = filter(sumTbl, !(gene == intCtrl)), aes(x = gene, y = relExp, colour = Name, group = Name)) +
                    geom_jitter(#aes(pch = rep),
                                position = position_dodge(width = 0.5), size = 3, alpha = 0.5)+
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
                    labs(x = "Gene", y = paste0("Relative Expression to ", intCtrl),
                         colour = "Sample", pch = "Replicate")+
                    theme_classic(base_size = fontSize)#+
                    #scale_color_manual(values = c("navyblue","maroon",rep("springgreen3",4),
                    #              rep("dodgerblue",2),"purple"))
                    
  
  if (logScale) {print(relplot + scale_y_continuous(trans = "log10"))}else{print(relplot)}
  
  fcplot <- ggplot(data = filter(fcTbl, !(gene == intCtrl)),
                   aes(x = gene, y = fc, colour = Name, group = Name))+
    geom_jitter(#aes(pch = rep),
                position = position_dodge(width = 0.5), size = 3, alpha = 0.5)+
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
    labs(x = "Gene", y = paste0("Fold Change to ", ctrlGroup),
         colour = "Sample", pch = "Replicate")+
    geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed")+
    theme_classic(base_size = fontSize) #+
    #scale_color_manual(values = c("navyblue","maroon",rep("springgreen3",4),
    #                              rep("dodgerblue",2)))# +
    #scale_y_continuous(limits = c(0,5))
  
  if(logScale){print(fcplot + scale_y_continuous(trans = "log10"))}else{print(fcplot)}
}

output <- group_by(sumTbl, Name, gene)
exp_sum <- summarise(output, mean_cp = mean(mean), rel_exp = mean(relExp))
exp_tbl <- dcast(exp_sum, Name ~ gene, value.var = "rel_exp")
cp_tbl <- dcast(exp_sum, Name ~ gene, value.var = "mean_cp")
write.csv(x = exp_tbl, file = "summary_table_exp.csv")
write.csv(x = cp_tbl, file = "summary_table_cp.csv")
