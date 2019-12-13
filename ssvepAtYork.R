##########################################
#######          Paths        ############
##########################################

## analyse the behavioural data

#ggsave(pl, filename = "figures/condition curves base.tiff", width = sc*8, height = sc*5, units = "cm", dpi = 600, scale= 1)

zscore = 0
analyseOneSub = 0
showAllTrials = 0 
subi = 23
analyseEEG = 0
runAnova = 0

# main folder
# setwd("C:/Users/Richard Naar/Documents/dok/ssvep/Visit to York/DATA") 
setwd("C:/Users/Richard Naar/Documents/dok/ssvep/Visit to York/EEG data/behavioural") 

# behavioural data
behPath = "/Users/Richard Naar/Documents/dok/ssvep/Visit to York/EEG data/behavioural"
# silmaandmete alamkaust
eegPath = ""

plotSave = "C:/Users/Richard Naar/Documents/dok/ssvep/Visit to York/Plots"

# extract file names for begavioural data
filesB <- list.files(path = behPath, pattern = "*.csv", full.names = TRUE)

if (analyseOneSub == 0) {

# subject loop
allDat <- data.frame()
for (subi in 1:length(filesB)){
# reading data in
behDat <- read.csv(filesB[subi], stringsAsFactors=FALSE)

# .6 = 52.5%, 1 = 87.5% 
behDat$trials_2.intensity <- (behDat$trials_2.intensity * 60)/0.6
behDat$trials_2.intensity.z <- scale(behDat$trials_2.intensity , center=TRUE, scale = TRUE)

allDat = rbind(allDat, behDat)

}


NAN <- as.data.frame( is.na(allDat) )
allDat <- subset(allDat, NAN$trials_2.intensity == FALSE)  #!= 'NA'

for (ti in 1:length(allDat$trials_2.label)){
  if (allDat$trials_2.thisRepN[ti] < 6) {
    allDat$label_2[ti] <- '1...5'
  } else if (allDat$trials_2.thisRepN[ti] < 21) { 
    allDat$label_2[ti] <- '6...20'
  } else if (allDat$trials_2.thisRepN[ti] < 36) {
    allDat$label_2[ti] <- '21...35'
  } else if (allDat$trials_2.thisRepN[ti] < 51) {
    allDat$label_2[ti] <- '36...50'
  } else if (allDat$trials_2.thisRepN[ti] < 66) {
    allDat$label_2[ti] <- '51...65'
  } else {
    allDat$label_2[ti] <- 'NA' 
  }

}

allDat$label_2 <- ordered(allDat$label_2, levels = c("1...5", "6...10", "11...15","21...26","26...30", "35...65"))

#allDat$label_2 <- ordered(allDat$label_2, levels = c("1...5", "6...20", "21...35","36...50", "51...65"))


allDat$trials_2.label <- as.factor(allDat$trials_2.label)

levels(allDat$trials_2.label)[levels(allDat$trials_2.label)=="high"] <- "Cued high"
levels(allDat$trials_2.label)[levels(allDat$trials_2.label)=="low"] <- "Cued Low"
levels(allDat$trials_2.label)[levels(allDat$trials_2.label)=="high50"] <- "Non-cued (high)"
levels(allDat$trials_2.label)[levels(allDat$trials_2.label)=="low50"] <- "Non-cued (low)"

levels(allDat$trials_2.label) 

allDat$cued[allDat$trials_2.cueText == 'fast'] <- 'Cued'   
allDat$cued[allDat$trials_2.cueText == 'slow'] <- 'Cued'   
allDat$cued[allDat$trials_2.cueText == '?'] <- 'Non-Cued'   


# means
require(Rmisc)
# meanAllDat <- summarySE(allDat, measurevar="trials_2.intensity.z", groupvars=c("trials_2.label", "label_2")) #
meanAllDat <- summarySE(allDat, measurevar="trials_2.intensity", groupvars=c("trials_2.label", "label_2")) #


# head(meanAllDat)
# meanAllDat1 <- subset(meanAllDat, trials_2.label == "Non-cued (high)" | trials_2.label == "Non-cued (low)" | trials_2.label == "Cued high")
# meanAllDat1 <- subset(meanAllDat, trials_2.label == "Non-cued (high)" | trials_2.label == "Non-cued (low)" | trials_2.label == "Cued high")


library(wesanderson)
require(ggplot2)
if (showAllTrials == 0) {

# plot windows 

ggplot(meanAllDat, aes(x=label_2, y= trials_2.intensity,color = trials_2.label,  group = trials_2.label)) + # , linetype=trials_2.label 
  geom_line(size = 1.3) +
  geom_errorbar(aes(ymin=trials_2.intensity-ci, ymax=trials_2.intensity+ci),
                width=.1, size= 1)+ # , color='red'
  theme_bw()+
  #  facet_wrap(~ KI) +
  theme(text = element_text(size=12, face = 'bold'))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.ticks = element_line(size = 1)) +
 #   theme(axis.line = element_line(color = 'black', size = 1)) +
  xlab("Trials") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size="12", angle = 0, hjust = 0.5, face = 'bold',color = 'black'))+
  ylab("Contrast thresholds (lower == higher acuity) ") +
  theme(axis.text.y = element_text(size="12", angle = 0, hjust = 0, face = 'bold',color = 'black')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  expand_limits(y=c(15, 35)) +                        # Expand y range
  theme(legend.title = element_blank(), legend.text=element_text(size=12)) +
  theme(legend.position= c(0.18,0.18)) +
  scale_colour_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00")) 

  # ratios 930 520 (582 423)

  #levels(meanAllDat$trials_2.label)  
  #"Cued high"  "Non-cued (high)" "Cued Low"  "Non-cued (low)" 
    
  # scale_colour_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00")) 
  # #999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" # A colorblind-friendly palette

  #scale_color_manual(values = wes_palette("Rushmore1",4,"discrete")) 
  # scale_color_manual(values = wes_palette("Zissou1",4,"discrete")) 

dataAll <- subset(allDat, label_2 == '35...65') 
allDat$label_2 <- as.factor(allDat$label_2)
#dataAll <- subset(allDat, label_2 == '61...65')  

meanLast <- summarySE(dataAll, measurevar="trials_2.intensity.z", groupvars=c("trials_2.label")) #

tsize= 12 # text size

ggplot(meanLast, aes(x=trials_2.label, y=trials_2.intensity.z)) + 
  geom_bar(fill = c( "#D55E00", "#009E73", "#F0E442", "#0072B2"),# 'white',
           stat="identity",
           width=0.6,
           colour="black",
           size = 1.3) +
  geom_point(size=2, position=position_dodge(1)) +
  theme_bw()+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.ticks = element_line(size = 1)
  ) +
  theme(text = element_text(size=tsize)) +
  ylab("Contrast thresholds (lower == higher acuity) ") +
  xlab(" ") +
#     coord_cartesian(ylim = c(15, 30)) +
  theme(axis.text.y = element_text(size=tsize,face="bold", angle = 0, hjust = 0), 
        axis.title=element_text(size=tsize,face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size=tsize,face="bold", angle = 0, hjust = 0.5)) +  # element_text(size="12", angle = 0, hjust = 0)
  geom_errorbar(aes(ymin=trials_2.intensity.z-ci, ymax=trials_2.intensity.z+ci), color = 'black', width=.25, size= 1.25, position=position_dodge(.1)) +
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid"))

    #if (runAnova == 1){
#  
    
dataAll$cued[dataAll$trials_2.cueText == 'fast'] <- 'Cued'   
dataAll$cued[dataAll$trials_2.cueText == 'slow'] <- 'Cued'   
dataAll$cued[dataAll$trials_2.cueText == '?'] <- 'Non-Cued'   

dataAll$cued <- as.factor(dataAll$cued)


Last <- subset(allDat, trials_2.thisN > 30)  

library(ez)
options(contrasts=c("contr.sum", "contr.poly"))

require(Rmisc)

meanLast <- summarySE(Last, measurevar="trials_2.intensity.z", groupvars=c("cued", "trials_2.frex", "participant")) #


ezANOVA(data=meanLast, dv=.(trials_2.intensity.z), wid=.(participant), within=.(trials_2.frex, cued), detailed = TRUE, type=3)
#ezANOVA(data=allDat, dv=.(trials_2.intensity.z), wid=.(participant), within=.(trials_2.label, label_2), detailed = TRUE, type=3)


meanLast$trials_2.frex <- as.factor(meanLast$trials_2.frex)
meanLast$participant <- as.factor(meanLast$participant)

#  allDat1 <- subset(allDat, label_2 == '51...65')

meanLast <- summarySE(Last, measurevar="trials_2.intensity.z", groupvars=c("trials_2.label","participant")) #
  
with(meanLast, pairwise.t.test(trials_2.intensity.z, trials_2.label, p.adjust.method="BH", paired=T))


## non-parametric tests
# kruskal.test(dataAll$trials_2.intensity.z ~ dataAll$trials_2.label)

x = subset(dataAll, trials_2.label == 'Cued Low')
y = subset(dataAll, trials_2.label == 'Non-cued (low)')
wilcox.test(x$trials_2.intensity.z, y$trials_2.intensity.z,paired=TRUE)

## violin plot
#install.packages('ggpubr')

# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list(c("Cued Low","Non-cued (low)" ), c("Cued high", "Non-cued (high)"))

require(ggpubr)
ggviolin(dataAll, x = 'trials_2.label', y = 'trials_2.intensity', fill = 'trials_2.label',
         palette = c("#0072B2","#E7B800","#00AFBB", "#FC4E07"), # "#0072B2","#E7B800","#00AFBB", "#FC4E07"
         add = "boxplot", add.params = list(fill = "white"), xlab = FALSE, 
         ylab = 'Contrast thresholds (lower == higher acuity)')+
  theme(legend.title = element_blank(),
        axis.text.x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = 't.test', 
                     tip.length = 0.03, paired = TRUE) # wilcox.test t.test Add significance levels

## ggpaired
meanLast2 <- summarySE(Last, measurevar="trials_2.intensity.z", groupvars=c("cued", "trials_2.frex", "trials_2.label", "participant")) #

meanLast2$trials_2.frex <- as.factor(meanLast2$trials_2.frex)
levels(meanLast2$trials_2.frex) 
levels(meanLast2$trials_2.frex) <- c('High', 'Low')
  
ggpaired(meanLast2, x = 'cued', y = 'trials_2.intensity.z', 
         color = 'trials_2.label', line.color = "gray", line.size = 0.1,
         facet.by = "trials_2.frex",
         palette = "jco",
         ylab = 'Contrast thresholds (lower == higher acuity)',
         xlab = FALSE) +
         theme(legend.title = element_blank(),
         axis.text.x = element_blank(),
         panel.border = element_blank(),
         axis.ticks.y = element_line(size = 1),
         axis.ticks.x = element_blank()) 
#  stat_compare_means(paired = TRUE, method = 't.test') 

#}  
  
## other  
  
#describe(allDat$trueRT)

#allDatClean <- subset(allDat, trueRT != 'no response')
#allDatClean$trueRT <- as.numeric(allDatClean$trueRT) 

#hist(allDatClean$trueRT, 50)  
#describe(allDatClean$trueRT)
#names(allDatClean)


# Change color by groups 
#dp <- ggplot(allDatClean, aes(x=trials_2.label, y=trials_2.intensity, fill=trials_2.label)) + 
#  geom_violin(trim=FALSE)+
#  geom_boxplot(width=0.1, fill="white")+
#  facet_wrap(~ label_2) +
#  labs(title="Plot of length  by dose",x="Dose (mg)", y = "Length")
#dp + scale_color_manual(values = wes_palette("Rushmore1",4,"discrete")) 


} else {
  
meanAllDat <- summarySE(allDat, measurevar="trials_2.intensity.z", groupvars=c("trials_2.label", "trials_2.thisRepN")) #

# Make the plot
ggplot(data=meanAllDat, aes(x=trials_2.thisRepN, y=trials_2.intensity.z, ymin=trials_2.intensity.z-ci, ymax=trials_2.intensity.z+ci, fill=trials_2.label, linetype=trials_2.label)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5) + 
#   scale_x_log10() + 
#    scale_y_log10() + 
  xlab("") + # as.expression(expression( paste("Radius (", R[500], ")") ))
  ylab("") +
  scale_fill_manual(values = wes_palette("Royal1",4,"discrete")) 

}

} else {

behDat <- read.csv(filesB[subi], stringsAsFactors=FALSE)
  
  # .6 = 52.5%, 1 = 87.5% 
behDat$trials_2.intensity <- (behDat$trials_2.intensity * 60)/0.6
behDat$trials_2.intensity.z <- scale(behDat$trials_2.intensity , center=TRUE, scale = TRUE)

lablesh <- subset(behDat, trials_2.label == 'high')
lablesl <- subset(behDat, trials_2.label == 'low')
lablesh5 <- subset(behDat, trials_2.label == 'high50')
lablesl5 <- subset(behDat, trials_2.label == 'low50')

if (zscore) {
  
  matplot(lablesh$trials_2.thisN, cbind(lablesh$trials_2.intensity.z,lablesl$trials_2.intensity.z,  
                                        lablesh5$trials_2.intensity.z, lablesl5$trials_2.intensity.z),
          xlab = "Trial", ylab = "Change in contrast (%)",xlim = c(-15,270), ylim = c(-5, 3),
          col=c(1,6,1,6), pch = c(19, 18, 21, 22),  cex=1) # pch=20 , , 
  
  
  legend("bottomright", legend = c('High','Low','Random (High)','Random (Low)'), 
         col=c(1,6,1,6), pch = c(19, 18, 21, 22)) # c("darkgreen","darkred", "green", "red")
  
} else {
  matplot(lablesh$trials_2.thisN, cbind(lablesh$trials_2.intensity,lablesl$trials_2.intensity,  
                                        lablesh5$trials_2.intensity, lablesl5$trials_2.intensity),
          xlab = "Trial", ylab = "Change in contrast (%)",xlim = c(-15,270), ylim = c(-15, 60),
          col=c(1,6,1,6), pch = c(19, 18, 21, 22),  cex=1) # pch=20 , , 
  
  
  legend("bottomright", legend = c('High','Low','Random (High)','Random (Low)'), 
         col=c(1,6,1,6), pch = c(19, 18, 21, 22)) # c("darkgreen","darkred", "green", "red")
}

  
}


## EEG

if (analyseEEG == 1){
  
imf <- c('F2-F1', 'F1+F2', '(2*F2)-(2*F1)','3*F1+F2', '(F1+F2)*2','3*F2 - F1')

eegPath = "/Users/Richard Naar/Documents/dok/ssvep/Visit to York/EEG data"
filesEEG <- list.files(path = eegPath, pattern = "*.csv", full.names = TRUE)

eegDat <- read.csv(filesEEG[15], stringsAsFactors=FALSE)
filesEEG[15]

#subDat <- data.frame()
#eegDat.z <- data.frame()
#for (subi in 1:length(unique(eegDat$SubId))){
#  subDat <- subset(eegDat, SubId == subi)
#  subDat$low8.z <- scale(subDat$low8 , center=TRUE, scale = TRUE)
#  subDat$F19.z <- scale(subDat$F19 , center=TRUE, scale = TRUE)
#  subDat$high30.z <- scale(subDat$high30 , center=TRUE, scale = TRUE)
#  eegDat.z <- rbind(eegDat.z, subDat)
#}

# install.packages('reshape')
library(reshape)
require(Rmisc)
require(ggplot2)
require(psych)

exclude = which(as.vector(unlist(tapply(eegDat$high30,eegDat$SubId, function(x) x %in% boxplot.stats(x)$out))))
eegCleanDat <- subset(eegDat[-exclude,])
eegCleanDat <- eegDat

#names(eegCleanDat)
meanEEG <- summarySE(eegDat, measurevar="high15", groupvars=c("cond")) # eegDatCleanDat
#describe(eegDat$low12) #describe(subset(rtm, cueInf == '50%')$reactionTime)

# Kernel Density Plot
# plot(density(subset(rtm, cueInf == '70%')$reactionTime), col='green')

names(meanEEG)[3] <- 'yVal'
#meanEEG$yVal <- meanEEG$yVal - mean(meanEEG$yVal)
#meanF <- summarySE(mdata, measurevar="value", groupvars=c("cond"))


ggplot(meanEEG, aes(x=cond, y= yVal )) + 
  geom_bar(width=0.8, position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=yVal-ci, ymax=yVal+ci),
                width=.2, size= 0.5,                    # Width of the error bars
                position=position_dodge(0.4), color='red')+
  theme_bw()+
  #  facet_wrap(~ KI) +
  theme(text = element_text(size=12))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.ticks = element_line(size = 1)) +
  theme(axis.line = element_line(color = 'black', size = 1)) +
  xlab("") +
  theme(axis.text.x = element_text(size="12", angle = 0, hjust = 0))+
  ylab("") +
  theme(axis.text.y = element_text(size="12", angle = 0, hjust = 0))
  

# eegDat$trial <- as.factor(eegDat$trial)

eegCleanDat$cued[eegCleanDat$cond == "Cued high"] <- 'Cued'   
eegCleanDat$cued[eegCleanDat$cond == "Cued low"] <- 'Cued'   
eegCleanDat$cued[eegCleanDat$cond == "Non-cued high"] <- 'Non-Cued'   
eegCleanDat$cued[eegCleanDat$cond == "Non-cued low"] <- 'Non-Cued'   

eegCleanDat$frex[eegCleanDat$cond == "Cued high"] <- 'High'   
eegCleanDat$frex[eegCleanDat$cond == "Cued low"] <- 'Low'   
eegCleanDat$frex[eegCleanDat$cond == "Non-cued high"] <- 'High'   
eegCleanDat$frex[eegCleanDat$cond == "Non-cued low"] <- 'Low'   


eegCleanDat$cued <- as.factor(eegCleanDat$cued)
eegCleanDat$frex <- as.factor(eegCleanDat$frex)


# install.packages('ez')
library(ez)
options(contrasts=c("contr.sum", "contr.poly"))
eegCleanDat$trial <- as.factor(eegCleanDat$trial)

meanClean <- summarySE(eegCleanDat, measurevar="high15", groupvars=c("frex", "cued", "cond", "SubId")) # eegDatCleanDat
#meanClean <- subset(meanClean, SubId != 3 & SubId != 22)
ezANOVA(data=meanClean, dv=.(high15), wid=.(SubId), within=.(cued, frex), type=3) 
#(data=eegCleanDat, dv=.(F11), wid=.(SubId), within=.(cued, frex), type=3) 


#ezANOVA(data=eegDat, dv=.(high30), wid=.(SubId), within=.(cond), type=3) 


with(meanClean, pairwise.t.test(high30, cond, p.adjust.method="BH", paired=T))


names(eegCleanDat)

N <- summarySE(meanClean, measurevar="N", groupvars=c("cond"))

}

## other

require(dplyr)

ratioDat <- data.frame()

N <- length(unique(allDat$participant))
trials <- length(unique(allDat$trials_2.thisRepN))
nRows <- N*(trials-24)

ratioDat <- setNames(data.frame(matrix(ncol = 4, nrow = nRows)), c("Subject", "Trial", "LowRatio", "HighRatio") )
counter = 0
for (subi in 1:N) {
  subDat <- subset(allDat, participant == subi)
  for (ti in 25:trials) {
    counter = counter + 1
    
    trialDat <- subset(subDat, trials_2.thisRepN == ti)
    trialOrder <- (order(trialDat$trials_2.label) )
    trialDat <- trialDat[trialOrder,]
    
    ch <- trialDat$trials_2.intensity.z[1] #"Cued high"       
    nch <- trialDat$trials_2.intensity.z[2] #"Non-cued (high)" 
    cl <- trialDat$trials_2.intensity.z[3] #"Cued Low"        
    nch <- trialDat$trials_2.intensity.z[4] #"Non-cued (low)" 
    
    #ch nch 
    ratioDat$HighRatio[counter] <- ch - nch
    ratioDat$LowRatio[counter] <- cl - nch
    ratioDat$Subject[counter] <- subi
    ratioDat$Trial[counter] <- ti
    
  }
}

#ratioDat$HighRatio <- exp(ratioDat$HighRatio)

x <- ratioDat$LowRatio
y <- ratioDat$HighRatio 

plot( x,type="l",col="red", lwd = 2)

lines(ratioDat$HighRatio,col="green", lwd = 2)

cor(x,y)

#acf(ratioDat[,3:4], lag.max = 20,
#    type = c("correlation"), # "correlation", "covariance", "partial"
#    plot = TRUE, demean = TRUE)

# pm
dataAll <- data.frame()
counter = 0
for (subi in 1:length(unique(allDat$participant))) {
  
  subDat <- subset(allDat, participant == subi)
  
  conds = c("Cued high", "Non-cued (high)", "Cued Low", "Non-cued (low)")
  
  for (condi in 1:length(conds)) {
#    counter = counter + 1
    
    data <- subset(subDat, trials_2.label == conds[condi])
    #data$localz <- scale(data$trials_2.intensity , center=TRUE, scale = TRUE)
   data <- subset(data, trials_2.thisRepN > 5)
    
    ord = order(data$trials_2.intensity)
    data = data[ord, ]
    
    #p <- c(20, 40, 60, 80, 100)/100
    #Qs = quantile(data$trials_2.intensity, probs = p)
    
#    for (tri in 1:length(data$localz)){
#      zi = data$localz[tri]
#      
#      
#      if (zi <= Qs[1]) { 
#        data$quantile[tri] <- '20'
    #   } else if (zi <= Qs[2]) {
    #     data$quantile[tri] <- '40'
    #   } else if (zi <= Qs[3]) {
    #     data$quantile[tri] <- '60'
    #   } else if (zi <= Qs[4]) {
    #     data$quantile[tri] <- '80'
    #   } else if (zi <= Qs[5]) {
    #     data$quantile[tri] <- '100'
    #   } else {
    #     data$quantile[tri] <- 'NA'
    #   print(subi)
    #   print(tri)
    #   print(zi)
    #   
    #     
    #   }
    # }
    
    
    for (tri in 1:length(data$trials_2.intensity)){
      zi = data$trials_2.intensity[tri]
      
      
      if (zi <= 10) { 
        data$delta[tri] <- '5'  #'0-10'
      } else if (zi <= 17) {
        data$delta[tri] <- '15' # '8-23'
      } else if (zi <= 23) {
        data$delta[tri] <- '20' # '17-23'
      } else if (zi <= 26) {
        data$delta[tri] <- '24' # '23-26'
      } else if (zi <= 29) {
        data$delta[tri] <- '28' # '26-29'
      } else if (zi <= 45) {
        data$delta[tri] <- '37' # '29-45'
      } else if (zi <= 60) {
        data$delta[tri] <- '53' # '45-60'
      } else {
        data$delta[tri] <- 'NA'
        print(subi)
        print(tri)
        print(zi)
        
        
      }
    }
    
    
    dataAll = rbind(dataAll, data)
  }
}

require(Rmisc)


dataAll$delta <- as.factor(dataAll$delta)
levels(dataAll$delta)
#dataAll$quantile <- ordered(dataAll$quantile, levels = c("12.5", "25", "37.5", "50", "62.5", "75","87.5" ,"100"  ))
#dataAll$delta <- ordered(dataAll$delta, levels = c("10", "20","30","40", "50", "60"  ))
#dataAll$delta <- ordered(dataAll$delta, levels = c("10", "15","20","25","30","35","40","45","50","55", "60"  ))
#dataAll$delta <- ordered(dataAll$delta, levels = c("15", "25","30","35","40","50","60"  ))
#dataAll$delta <- ordered(dataAll$delta, levels = c("1-17", "17-23","23-29","29-60"))
# dataAll$delta <- ordered(dataAll$delta, levels = c("8", "20","26","44"))



# meanContrast <- summarySE(dataAll, measurevar='trials_2.response', groupvars=c('delta', 'trials_2.label','participant'))

require(plyr)


meanContrast = ddply(dataAll, .( delta, trials_2.label, participant, trials_2.thisRepN, trials_2.intensity), summarize,
                     NumPos = sum(trials_2.response),
                     N      = length(trials_2.response))

meanContrast = ddply(dataAll, .( delta, trials_2.label), summarize,
               NumPos = sum(trials_2.response),
               N      = length(trials_2.response))


meanContrast$trials_2.intensity <- as.numeric(meanContrast$trials_2.intensity )


write.table(meanContrast, paste0(getwd(),"/contrastId3.txt"), sep="\t", dec = ",") # meanContrast



# plot(juku$quantile, juku$trials_2.response)

library(ggplot2)

#pData <- ggplot(juku, aes(x = quantile, y = trials_2.response)) + 
 # geom_point()


ggplot(juku, aes(x=delta, y= trials_2.response, color = trials_2.label,  group = trials_2.label)) + # , linetype=trials_2.label 
  geom_line(size = 1) +
  geom_point(size = 2) +
#  geom_errorbar(aes(ymin=trials_2.response-ci, ymax=trials_2.response+ci),
#                width=.1, size= 1)+ # , color='red'
  theme_bw()+
  #  facet_wrap(~ KI) +
  theme(text = element_text(size=12, face = 'bold'))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.ticks = element_line(size = 1)) +
  #   theme(axis.line = element_line(color = 'black', size = 1)) +
  xlab("Contrast") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size="12", angle = 0, hjust = 0.5, face = 'bold',color = 'black'))+
  ylab("Percent correct") +
  theme(axis.text.y = element_text(size="12", angle = 0, hjust = 0, face = 'bold',color = 'black')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  #expand_limits(y=c(15, 35)) +                        # Expand y range
  theme(legend.title = element_blank(), legend.text=element_text(size=12)) +
  #theme(legend.position= c(0.18,0.18)) +
  scale_colour_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00")) 


parme

juku$nYes <- juku$trials_2.response * juku$N   

N <- sum(juku$N)

fit <- quickpsy(data, quantile, nYes, N)




juku <- summarySE(dataAll, measurevar='trials_2.response', groupvars=c('quantile'))


library(quickpsy)

fit <- quickpsy(data, quantile, trials_2.response)
plot(fit)


dataAll$quantile <- as.numeric(dataAll$quantile)
dataAll$trials_2.response <-  as.numeric(dataAll$trials_2.response) 


fit <- quickpsy(dataAll, quantile, trials_2.response,prob = 0.80) 
plot(fit)

plotthresholds(fit)

# end
