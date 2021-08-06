#############################################################
# 3.SC_JAGS_viewoutput.R
# created by Joanna Burgar, 8-May-2018; modified 09-Apr-2020
# updated by Joanna Burgar, 25-Jul-2021 for Leroy Gonzales (streamlined)
# script for viewing JAGS output
#############################################################


######################################
#Load Packages
list.of.packages <- c("xtable", "plyr","tidyverse", "Cairo", "coda")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

######################################
#load files

info <- list.files(recursive = T, pattern="*.RData") # pulls in all mcmc files, from subfolders
info <- info[!grepl("dat|map", info)] # removes the dat and density map (large) file
for(i in 1:length(info)) load(info[[i]]) 

StudyAreas <- substr(info, 5, str_locate(info, "_")[1]-1)

######################################
#create function to load files and write output to table

#function to create output table for JAGS output
get_JAGS_output <- function(filename){
  out <- filename
  s <- summary(window(out, start = 41000))
  gd <- gelman.diag(window(out, start = 41000),multivariate = FALSE)
  output_table <- rbind(as.data.frame(t(s$statistics)),
                        as.data.frame(t(s$quantiles)),
                        as.data.frame(t(gd$psrf)))
  return(output_table)	
}

##################################################################################################
###---format and aggregate model data
info.mc <- str_replace(info.mc, pattern = "_2020_", replacement = "_")
info.mc <- str_replace(info.mc, pattern = "out/", replacement = "")
JAGS.output <- lapply(info.mc, function(i) get_JAGS_output(get(i)))
JAGS.output <- do.call(rbind, JAGS.output)
glimpse(JAGS.output)

JAGS.output.df <- as.data.frame(JAGS.output)
summary(JAGS.output.df)
JAGS.model <- rep(info.mc, times=1, each=11)
JAGS.output.df$model <- JAGS.model
JAGS.StudyArea <- rep(StudyAreas, times=1, each=11)
JAGS.output.df$Area <- JAGS.StudyArea
head(JAGS.output.df,20)

estimates <- c("Mean", "SD", "Naive SE","Time-series SE", "CI_2.5","CI_25", "CI_50", "CI_75", "CI_97.5", "gd.point","gd.upper")
JAGS.output.df$estimates <- as.factor(rep(estimates, nrow(JAGS.output.df)/11))

write.csv(JAGS.output.df, file=paste("Koala_modeloutput.csv"))

###--- manipulate data to be in correct format for ggplot graphs
data <- read.csv("Koala_modeloutput.csv") 

glimpse(data)
data <- data[c(-1)]
head(data)

data2 <- gather(data, key=parameter, value=value, 1:5)
# data checks
glimpse(data2)
data2$parameter <- as.factor(data2$parameter)
str(data2)
summary(data2)
table(data2$estimate, data2$parameter)

data.plot <- spread(data2, estimates, value)
levels(data.plot$parameter)
data.plot$parameter <- factor(data.plot$parameter, levels = c("sigma", "lam0", "psi", "D", "N"))

data2$parameter <- mapvalues(data2$parameter, 
                             from = c("sigma", "lam0", "psi", "D", "N"),
                             to = c("Sigma", "Lam0","Psi","Density","Abundance"))
data2$parameter <- factor(data2$parameter, levels = c("Density","Abundance", "Lam0","Sigma", "Psi"))

density <- filter(data2, parameter=="Density")
density.plot <- spread(density, estimates, value)
density.plot$model <- as.factor(substr(density.plot$model,5,9))

abundance <- filter(data2, parameter=="Abundance")
abundance.plot <- spread(abundance, estimates, value)
abundance.plot$model <- as.factor(substr(abundance.plot$model,5,9))

lam0 <- filter(data2, parameter=="Lam0")
lam0.plot <- spread(lam0, estimates, value)
lam0.plot$model <- as.factor(substr(lam0.plot$model,5,9))

sigma <- filter(data2, parameter=="Sigma")
sigma.plot <- spread(sigma, estimates, value)
sigma.plot$model <- as.factor(substr(sigma.plot$model,5,9))

psi <- filter(data2, parameter=="Psi")
psi.plot <- spread(psi, estimates, value)
psi.plot$model <- as.factor(substr(psi.plot$model,5,9))

data3 <- rbind(density.plot, abundance.plot,lam0.plot, sigma.plot, psi.plot)
data3$CV <- data3$SD/data3$Mean
data3$Prior <- ifelse(data3$model=="KUI1" |data3$model=="K_KUI1","Uniformed","Informed")
summary(data3$CV)
summary(data3$gd.point) # Gelman-Rubin diagnostic <1.1 suggests model convergence
summary(data3$gd.upper) # Gelman-Rubin diagnostic <1.1 suggests model convergence
write.csv(data3, file="Koala_JAGS_modeloutput.csv")

##########################################################################
##########################################################################
# Plot the Mean and standard deviation for each model by species:  

data3 <- read.csv("Koala_JAGS_modeloutput.csv")
data3 <- data3[,-1]
head(data3)
summary(data3)
data3$Prior <- factor(data3$Prior, levels = c("Uniformed","Informed"))
levels(data3$model) # this may have to be manually changed depending on priors and number of models
data3$model2 <- mapvalues(data3$model, 
                          from = c("KSI1","KUI1","KWI1", "KWI2"),
                          to = c("SI1", "UI", "WI1", "WI2"))

p = ggplot(data3[data3$parameter!="Abundance" & data3$parameter!="Psi" ,], aes(model2, CI_50))+
  geom_point(colour="white", shape=21, size=10,
             aes(fill=Prior))+
  scale_fill_manual(values=c("black","grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name="Median ± 50% Credible Intervals") +
  geom_linerange(aes(model2, ymin = CI_25, ymax = CI_75)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14))+
  facet_wrap(~parameter+Area, scales="free_y", ncol=2)


###--- density
p.density <- ggplot(data3[data3$parameter=="Density",], aes(x = as.factor(model2), y = CI_50)) +
  geom_point(aes(shape=Prior, colour=model2),size=10)+
  scale_fill_manual(values=c("grey","black","grey","grey")) +
  geom_linerange(aes(ymin = CI_2.5, ymax = CI_97.5)) +
  theme_bw() +
  ylab(expression(paste("Estimated Density (ha); Median ± 50% BCI"))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14))+
  facet_wrap(~Area)


Cairo(file="Figure_Density.PNG", 
      type="png",
      width=2500, 
      height=1600, 
      pointsize=15,
      bg="white",
      dpi=300)
p.density
dev.off()

density = ggplot(data3[data3$parameter=="Density",], 
              aes(x = as.factor(model2), y = CI_50))+
  geom_point(colour="white", shape=21, size=10,aes(fill=Prior))+
  scale_fill_manual(values=c("black","grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(expression(paste("Median Density (individuals/ha) ± 50 BCI"))) +
  geom_linerange(aes(model2, ymin = CI_25, ymax = CI_75)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=14))+
  theme(legend.position = "none", strip.text = element_text(size=16))+
  facet_wrap(~Area)


ggsave("Density_50BCI.png",plot=density, device="png", dpi=300,
       width=11, height=8)


###--- lam0
p.lam0 <- ggplot(data3[data3$parameter=="Lam0",], aes(x = as.factor(model2), y = CI_50)) +
  geom_point(colour="white", shape=21, size=10,aes(fill=Prior))+
  scale_fill_manual(values=c("black","grey")) +
  geom_linerange(aes(ymin = CI_2.5, ymax = CI_97.5)) +
  theme_bw() +
  ylab(expression(paste("Estimated ",lambda[0],"; Median ± 95% Credible Intervals"))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14))+
  facet_wrap(~Area)

Cairo(file="Figure_Lam0.PNG", 
      type="png",
      width=2500, 
      height=1600, 
      pointsize=15,
      bg="white",
      dpi=300)
p.lam0
dev.off()


lam0 = ggplot(data3[data3$parameter=="Lam0",], 
              aes(x = as.factor(model2), y = CI_50))+
  geom_point(colour="white", shape=21, size=10,aes(fill=Prior))+
  scale_fill_manual(values=c("black","grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(expression(paste("Median ", lambda[0]," ± 50 BCI"))) +
  geom_linerange(aes(model2, ymin = CI_25, ymax = CI_75)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=14))+
  theme(legend.position = "none", strip.text = element_text(size=16))+
  facet_wrap(~Area)

ggsave("Lam0_50BCI.png",plot=lam0, device="png", dpi=300,
       width=11, height=8)


###--- sigma
p.sigma <- ggplot(data3[data3$parameter=="Sigma",], aes(x = as.factor(model2), y = CI_50*100)) +
  geom_point(colour="white", shape=21, size=10,aes(fill=Prior))+
  scale_fill_manual(values=c("black","grey")) +
  geom_linerange(aes(ymin = CI_2.5*100, ymax = CI_97.5*100)) +
  theme_bw() +
  ylab(expression(paste("Estimated ",sigma," (m); Median ± 95% Credible Intervals"))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14))+
  facet_wrap(~Area)

  
Cairo(file="Figure_Sigma.PNG", 
      type="png",
      width=2500, 
      height=1600, 
      pointsize=15,
      bg="white",
      dpi=300)
p.sigma
dev.off()


sigma = ggplot(data3[data3$parameter=="Sigma",], 
              aes(x = as.factor(model2), y = CI_50*100))+
  geom_point(colour="white", shape=21, size=10,aes(fill=Prior))+
  scale_fill_manual(values=c("black","grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(expression(paste("Median ", sigma," (m) ± 50 BCI"))) +
  geom_linerange(aes(model2, ymin = CI_25*100, ymax = CI_75*100)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=14))+
  theme(legend.position = "none", strip.text = element_text(size=16))+
  facet_wrap(~Area)

ggsave("Sigma_50BCI.png",plot=sigma, device="png", dpi=300,
       width=11, height=8)


report = ggplot(data3[data3$parameter!="Abundance" & data3$parameter!="Psi",], 
               aes(x = as.factor(model2), y = CI_50))+
  geom_point(colour="white", shape=21, size=10,aes(fill=Prior))+
  scale_fill_manual(values=c("black","grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(expression(paste("Median ± 50 BCI"))) +
  geom_linerange(aes(model2, ymin = CI_25, ymax = CI_75)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=14))+
  theme(legend.position = "none", strip.text = element_text(size=16))+
  facet_wrap(~parameter+Area, scales="free_y",ncol=2)


ggsave("DLS_50BCI.png",plot=report, device="png", dpi=300,
       width=11, height=20)
