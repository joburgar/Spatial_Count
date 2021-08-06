#############################################################
# 5.SC_Informed Sigma.R
# adapated from Chandler & Royle 2013 (appendix)
# created by Joanna Burgar, 08-May-2018
# updated by Joanna Burgar, 04-Aug-2021 for Leroy Gonzales (streamlined)
# script for estimating calculating vaguely informed priors
#############################################################

#########################################################################
# Informative sigma priors (from Chandler & Royle)

###--- home range and call attentuation to sigma prior
# home range size in ha
# detection range in metres

sigma_prior <- function (hrsize = hrszie, detect_min = detect_min, detect_max = detect_max){
  hrlength <- sqrt(hrsize/pi*10000)
  
  hr_detect_min_area <- (hrlength+detect_min)^2*pi / 10000 # 65.62 ha, and
  hr_detect_max_area <- (hrlength+detect_max)^2*pi / 10000 # 135.61 ha
  
  # Following Royle et. al (2011), and assuming a
  # chi-squared distribution with 2 degrees of freedom,
  # the range of sigma is given by
  
  sigma_lower <- sqrt(hr_detect_min_area*10000/pi)/sqrt(5.99)   # 187 m
  sigma_upper <- sqrt(hr_detect_max_area*10000/pi)/sqrt(5.99)  # 268 m
  
  # Assuming a grid spacing of 1 unit = 100 m, aim to have a prior with most
  # of the density between:
  
  pd_lower <- sigma_lower/100
  pd_upper <- sigma_upper/100 
  
  return(list(pd_lower, pd_upper))
}


###---
# now use the output to determine the range of dgamma by adjusting the last shape and scale
# for example, for a home range of 40 ha with a call detection range between 100-300 m
# hrsize <- 40
# detect_min <- 100
# detect_max <- 300

sigma_prior(hrsize = 40, detect_min = 100, detect_max = 300)
# [1] 1.866536
# [1] 2.683713
# so in this example, want the density of the simga prior between 1.866536 and 2.683713
# and the peak ~2.28

qgamma(c(0.001,0.5,0.999),60,27) #  1.439910 2.209889 3.215138 - tight home ranges centred on 40
curve(dgamma(x,60,27), col='black',xlim=c(0,5), ylim=c(0,2))

# for the model, will want to specify sigma ~ dgamma(60, 27) 
