# Main Analysis for mLI-RADS Paper
# Author: Anton Becker (@ASBecker)

source("../Data/csv.to.dataframe.R")

library(irr)
library(DescTools)
library(ROCR)
library(pROC)
library(ggplot2)

# Calculate LiRADS and mLiRADS

## Custom function for calculating LiRADS v2014-C. 
### Accepts list of 6 observations when using apply() on the dataframe
### Order: MAX_DM_LR MASS ART_HYP WASHOUT  CPS TRESHOLD

f.lirads <- function(x){
  # If-Else Blocks work in a row-wise fashion accroding to the ACR LI-RADS criteria
  x <- as.character(as.numeric(x))
  properties <- x[4:6]
  properties[is.na(properties)] <- 0 # Threshold if often NA
  sumrf <- sum(as.numeric(properties))

  thre <- x[6]#,
  diam <- as.numeric(as.character(x[1])) #,
  if (is.na(diam)) {return(NA)} else {
    enh <- x[3]
    if (is.na(enh)) {return(NA)} else {
    
    if (sumrf == 0) {
      # 1st row
      if (diam >= 20 && enh == 1) {score <- 4} else {score <- 3}} else if (sumrf == 1) {
      # 2nd row
      if (enh == 0) {
        if (diam < 20) {
          score <- 3} else if (diam >= 20) {score <- 4}
      } else if (enh == 1) {
        if (diam < 10) {
          score <- 4} else if (diam >= 20) {
          score <- 5} else {
            if (thre == 1) {score <- 5}
            else {score <- 4}
            }
        }
    } else if (sumrf >= 2) {
      #3rd row
      if (enh == 0 || diam < 10) {score <- 4} else {score <- 5}
    }
    }
  return(score)}
}

##### mLI-RADS algorithm:
f.mlirads <- function(x){
  x <- as.character(as.numeric(x))
  capswo <- sum(as.numeric(x[4:5])) #-4
  if (is.na(capswo)) {
    return(NA)
  } else {
    enh <- x[3]#,
    thre <- x[6]#,
    diam <- as.numeric(as.character(x[1]))#,

    if (capswo >= 1){
      if (diam >= 20){
        score <- 5
      } else if (is.na(thre)){
        return(NA) #"Missing Threshold"
      } else if (thre == 1){
        score <- 5
      } else {
        score <- 4
      }
    } else if (enh == 1){
      if (diam >= 20) {
        if (is.na(thre)) {
          return(NA)  #"Missing Threshold"
        } else if (thre == 1) {
          score <- 5} else {
            score <- 4}
      }
      else if (diam >= 15) {
        if (is.na(thre)) {
          return(NA)  #"Missing Threshold"
        } else if (thre == 1) {
          score <- 4} else {
            score <- 3}
      }
      else if (diam < 15) {
        if (is.na(thre)) {
          return(NA)  #"Missing Threshold"
        } else if (thre == 1) {
          score <- 3} else {
            score <- 2}
      }
    } else {
      score <- 2
    }
    return (score)
  }
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Calculate LiRADS and mLiRADS for every observed lesion

d.re["lirads1"] <- apply(d.re[c("MAX_DM_LR", "MASS", "ART_HYP", "WASHOUT",  "CPS", "TRESHOLD")], 1,  f.lirads)
d.re["m.lirads1"] <- apply(d.re[c("MAX_DM_LR", "MASS", "ART_HYP", "WASHOUT",  "CPS", "TRESHOLD")], 1,  f.mlirads)
d.re["lirads2"] <- apply(d.re[c("MAX_DM_LR2", "MASS2", "ART_HYP2", "WASHOUT2",  "CPS2", "TRESHOLD2")], 1,  f.lirads)
d.re["m.lirads2"] <- apply(d.re[c("MAX_DM_LR2", "MASS2", "ART_HYP2", "WASHOUT2",  "CPS2", "TRESHOLD2")], 1,  f.mlirads)
d.re["lirads3"] <- apply(d.re[c("MAX_DM_LR3", "MASS3", "ART_HYP3", "WASHOUT3",  "CPS3", "TRESHOLD3")], 1,  f.lirads)
d.re["m.lirads3"] <- apply(d.re[c("MAX_DM_LR3", "MASS3", "ART_HYP3", "WASHOUT3",  "CPS3", "TRESHOLD3")], 1,  f.mlirads)
d.re["lirads4"] <- apply(d.re[c("MAX_DM_LR4", "MASS4", "ART_HYP4", "WASHOUT4",  "CPS4", "TRESHOLD4")], 1,  f.lirads)
d.re["m.lirads4"] <- apply(d.re[c("MAX_DM_LR4", "MASS4", "ART_HYP4", "WASHOUT4",  "CPS4", "TRESHOLD4")], 1,  f.mlirads)

####### For d.complete -> Fill missing with reference LIRADS (cat 1 and 2 lesions) #########

d.re$lirads1[is.na(d.re$lirads1)] <- d.re$LIRADS_Category[is.na(d.re$lirads1)]
d.re$lirads2[is.na(d.re$lirads2)] <- d.re$LIRADS_Category2[is.na(d.re$lirads2)]
d.re$lirads3[is.na(d.re$lirads3)] <- d.re$LIRADS_Category3[is.na(d.re$lirads3)]
d.re$lirads4[is.na(d.re$lirads4)] <- d.re$LIRADS_Category4[is.na(d.re$lirads4)]

d.re$m.lirads1[is.na(d.re$m.lirads1)] <- d.re$LIRADS_Category[is.na(d.re$m.lirads1)]
d.re$m.lirads2[is.na(d.re$m.lirads2)] <- d.re$LIRADS_Category2[is.na(d.re$m.lirads2)]
d.re$m.lirads3[is.na(d.re$m.lirads3)] <- d.re$LIRADS_Category3[is.na(d.re$m.lirads3)]
d.re$m.lirads4[is.na(d.re$m.lirads4)] <- d.re$LIRADS_Category4[is.na(d.re$m.lirads4)]

# Calculate Fleiss' Kappa for LiRADS/mLiRADS
flk.lirads <- kappam.fleiss(d.re[c("lirads1", "lirads2", "lirads3", "lirads4")])
flk.ci.lirads <- KappaM(d.re[c("lirads1", "lirads2", "lirads3", "lirads4")], method = "Fleiss", conf.level=0.95)
flk.mlirads <- kappam.fleiss(d.re[c("m.lirads1", "m.lirads2", "m.lirads3", "m.lirads4")])
flk.ci.mlirads <- KappaM(d.re[c("m.lirads1", "m.lirads2", "m.lirads3", "m.lirads4")], method = "Fleiss", conf.level=0.95)

# Calculate Fleiss' Kappa for Imaging Features
icc.diam <- icc(d.re[c("MAX_DM_LR", "MAX_DM_LR2", "MAX_DM_LR3", "MAX_DM_LR4")], model="twoway", type="agreement")
flk.enh <- KappaM(d.re[c("ART_HYP", "ART_HYP2", "ART_HYP3", "ART_HYP4")], method = "Fleiss", conf.level=0.95)
flk.wo <- KappaM(d.re[c("WASHOUT", "WASHOUT2", "WASHOUT3", "WASHOUT4")], method = "Fleiss", conf.level=0.95)
flk.caps <- KappaM(d.re[c("CPS", "CPS2", "CPS3", "CPS4")], method = "Fleiss", conf.level=0.95)
flk.thre <- KappaM(d.re[c("TRESHOLD", "TRESHOLD2", "TRESHOLD3", "TRESHOLD4")], method = "Fleiss", conf.level=0.95)

####### Kappa - Size Threshold ########

diam.thres <- seq(7,45)

f.diamkappas <- function (feat.diam, diam.name, thresholds) {
  kappas <- matrix(ncol=2)
  feat.diam <- feat.diam[complete.cases(feat.diam),]
  for (x in thresholds) {
    feat <- feat.diam[which(feat.diam[,diam.name]>x),]
    feat <- feat[ , -which(names(feat) %in% c(diam.name))]
    k <- KappaM(feat, method="Fleiss")
    kappas <- rbind(kappas,matrix(c(x,k),ncol=2))
  }
  qplot(x=kappas[,1], y=kappas[,2], xlim=range(thresholds), ylim=c(0,0.7), 
       geom='line', xlab='Diameter [> x mm]', ylab='Fleiss Kappa')
  return(kappas)
}

par(mfrow=c(2,2))

f.diamkappas(d.re[c("ART_HYP", "ART_HYP2", "ART_HYP3", "ART_HYP4", "MAX_DM_LR")], 'MAX_DM_LR', diam.thres)

f.diamkappas(d.re[c("WASHOUT", "WASHOUT2", "WASHOUT3", "WASHOUT4", "MAX_DM_LR")], 'MAX_DM_LR', seq(7,34))

f.diamkappas(d.re[c("CPS", "CPS2", "CPS3", "CPS4", "MAX_DM_LR")], 'MAX_DM_LR', diam.thres)

f.diamkappas(d.re[c("TRESHOLD", "TRESHOLD2", "TRESHOLD3", "TRESHOLD4", "MAX_DM_LR")], 'MAX_DM_LR', seq(7,30))


v1 <- KappaM(d.re[which(d.re$MAX_DM_LR>20),c("ART_HYP", "ART_HYP2", "ART_HYP3", "ART_HYP4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))
v1 <- KappaM(d.re[which(d.re$MAX_DM_LR<=20),c("ART_HYP", "ART_HYP2", "ART_HYP3", "ART_HYP4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))


v1 <- KappaM(d.re[which(d.re$MAX_DM_LR>20),c("WASHOUT", "WASHOUT2", "WASHOUT3", "WASHOUT4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))
v1 <- KappaM(d.re[which(d.re$MAX_DM_LR<=20),c("WASHOUT", "WASHOUT2", "WASHOUT3", "WASHOUT4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))

v1 <- KappaM(d.re[which(d.re$MAX_DM_LR>20),c("CPS", "CPS2", "CPS3", "CPS4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))
v1 <- KappaM(d.re[which(d.re$MAX_DM_LR<=20),c("CPS", "CPS2", "CPS3", "CPS4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))

v1 <- KappaM(d.re[which(d.re$MAX_DM_LR>20),c("TRESHOLD", "TRESHOLD2", "TRESHOLD3", "TRESHOLD4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))
v1 <- KappaM(d.re[which(d.re$MAX_DM_LR<=20),c("TRESHOLD", "TRESHOLD2", "TRESHOLD3", "TRESHOLD4")], method="Fleiss", conf.level=0.95)
print(c(v1[1], v1[[1]]-v1[[2]]))