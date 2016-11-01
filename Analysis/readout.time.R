# Readout times for mLI-RADS Paper
# Author: Anton Becker (@ASBecker)

source('lirads.R')

library('reshape2')

f.mlirads.time <- function(x){
  x <- as.character(as.numeric(x))
  capswo <- sum(as.numeric(x[4:5])) #-4
  if (is.na(capswo)) {
    return(NA)
  } else {
    enh <- x[3]#,
    thre <- x[6]#,
    diam <- as.numeric(as.character(x[1]))#,
    timefactor <- 0
    
    if (diam >= 20){
      if (capswo>=1) {score <- 5; timefactor <- 3 }
      else if (thre ==1) {score <- 5; timefactor <- 4}
      else {score <- 4; timefactor <- 4}
    } else if (capswo>=1){
      if (is.na(thre)) {
        return(NA)  #"Missing Threshold" 
      } else if (thre == 1) {
        score <- 5
        timefactor <- 4
      } else {
        score <- 4
        timefactor <- 3
      }
    } else if (enh == 1){
      if (is.na(thre)) {
        return(NA)  #"Missing Threshold" 
      } else if (thre == 1) {
        score <- 4
        timefactor <- 5
      } else { 
        score <- 3
        timefactor <- 4
      } 
    }
    else {
      score <- 2
      timefactor <- 5
    }
    
    if (exists('score')){
      return(timefactor)
    } else {
      return(NA)
    }
  }
}

f.readouttime <- function (x) {
  timefactor <- x[1]
  time <- x[2]
  timefragment <- time/6
  readouttime <- timefragment*timefactor
  return(readouttime)
}

f.mmtoss <- function  (string) {
  mmss <- strsplit (string, ":", T)
  mm <- as.numeric (mmss[[1]][1])
  ss <- as.numeric (mmss[[1]][2])
  return (mm * 60 + ss)
}

###### Cleanup #######
d.re <- d.re[1:104,]

####### Convert MM:SS to SS ########
d.re['Timestring1'] <- apply(d.re['TIME_LR'], 1, as.character)
d.re['Timestring1'] <- apply(d.re['Timestring1'], 1, f.mmtoss)
d.re['Timestring2'] <- apply(d.re['TIME_LR2'], 1, as.character)
d.re['Timestring2'] <- apply(d.re['Timestring2'], 1, f.mmtoss)
d.re['Timestring3'] <- apply(d.re['TIME_LR3'], 1, as.character)
d.re['Timestring3'] <- apply(d.re['Timestring3'], 1, f.mmtoss)
d.re['Timestring4'] <- apply(d.re['TIME_LR4'], 1, as.character)
d.re['Timestring4'] <- apply(d.re['Timestring4'], 1, f.mmtoss)


######## Calculate Time Factors #########

d.re["m.lirads1.tf"] <- apply(d.re[c("MAX_DM_LR", "MASS", "ART_HYP", "WASHOUT",  "CPS", "TRESHOLD")], 1,  f.mlirads.time)
d.re["m.lirads2.tf"] <- apply(d.re[c("MAX_DM_LR2", "MASS2", "ART_HYP2", "WASHOUT2",  "CPS2", "TRESHOLD2")], 1,  f.mlirads.time)
d.re["m.lirads3.tf"] <- apply(d.re[c("MAX_DM_LR3", "MASS3", "ART_HYP3", "WASHOUT3",  "CPS3", "TRESHOLD3")], 1,  f.mlirads.time)
d.re["m.lirads4.tf"] <- apply(d.re[c("MAX_DM_LR4", "MASS4", "ART_HYP4", "WASHOUT4",  "CPS4", "TRESHOLD4")], 1,  f.mlirads.time)
d.re["lirads.tf"]  <- 5


######## Calculate readout times ########

d.re['liradstime1']  <- apply(d.re[c('Timestring1','lirads.tf')], 1, f.readouttime)
d.re['mliradstime1']  <- apply(d.re[c('Timestring1','m.lirads1.tf')], 1, f.readouttime)
d.re['liradstime2']  <- apply(d.re[c('Timestring2','lirads.tf')], 1, f.readouttime)
d.re['mliradstime2']  <- apply(d.re[c('Timestring2','m.lirads2.tf')], 1, f.readouttime)
d.re['liradstime3']  <- apply(d.re[c('Timestring3','lirads.tf')], 1, f.readouttime)
d.re['mliradstime3']  <- apply(d.re[c('Timestring3','m.lirads3.tf')], 1, f.readouttime)
d.re['liradstime4']  <- apply(d.re[c('Timestring4','lirads.tf')], 1, f.readouttime)
d.re['mliradstime4']  <- apply(d.re[c('Timestring4','m.lirads4.tf')], 1, f.readouttime)

###### Create Boxplot #########

avgtime.lirads <- rowMeans(d.re[c("liradstime1", "liradstime2", 
                                  "liradstime3", "liradstime4")], na.rm=TRUE)
avgtime.mlirads <- rowMeans(d.re[c("mliradstime1", "mliradstime2", 
                                   "mliradstime3", "mliradstime4")], na.rm=TRUE)
d.readtime <- data.frame(avgtime.lirads, avgtime.mlirads)
colnames(d.readtime) <- c('LiRads','mLiRads')

plottime.avg <- ggplot(data = melt(d.readtime), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme(legend.position='none') +
  ggtitle('Average') +
  labs(x=' ', y='Readout Time [s]') +
  scale_fill_brewer(palette="Set1")+
  scale_y_continuous(limits=c(-2,170))

d.indtimes <- data.frame(d.re[c("liradstime1", "mliradstime1",
                                "liradstime2", "mliradstime2",
                                "liradstime3", "mliradstime3",
                                "liradstime4", "mliradstime4")])
colnames(d.indtimes) <- c('R1-LR', 'R1-mLR', 'R2-LR', 'R2-mLR', 
                          'R3-LR', 'R3-mLR', 'R4-LR', 'R4-mLR')

plottime.indiv <- ggplot(data = melt(d.indtimes), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme(legend.position='none') +
  ggtitle('Individual') +
  labs(x=' ', y=' ') +
  scale_y_continuous(limits=c(-2,170))

qplot(d.readtime$LiRads)
qplot(d.readtime$mLiRads)

ks.test(d.readtime$LiRads, 'pnorm')
ks.test(d.readtime$mLiRads, 'pnorm')

paste(median(d.readtime$LiRads), range(d.readtime$LiRads), sep='; ')
paste(median(d.readtime$mLiRads, na.rm=TRUE), range(d.readtime$mLiRads, na.rm=TRUE), sep='; ')

wilcox.test(d.readtime$LiRads, d.readtime$mLiRads, paired=TRUE)

#tiff('Presentation/Figure_05.tiff', res=300, width=3000, height=1750, units='px')
multiplot(plottime.avg, plottime.indiv, cols=2)
#dev.off()
