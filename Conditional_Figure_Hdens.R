###################################################################
################ Plotting Conditional Probabilities ################
################ along an environmental gradient    ################
################ from MCMC output                   ################
###################################################################

#Set your working directory
setwd("C:/Users/Arielle/Desktop/NRC/Biodiversity Lab/Daily Planet Talks/TWS/2018/Multi-species Workshop/Mod_graphing")

#Load the necessary packages
library(ggplot2)
library(rstan)

####################### Step 1 #################################
#Load the MCMC object
#The object we are about to load is called "fit"
load("C:/Users/Arielle/Desktop/NRC/Biodiversity Lab/Daily Planet Talks/TWS/2018/Multi-species Workshop/Mod_graphing/Workshop_Mod.rda")

####################### Step 2 #################################
#Extract the beta coefficient estimates from MCMC
#The dimensions of the betas are iterations by natural parameter
#by covariate
b <- extract(fit, 'beta')[[1]]
dim(b)

####################### Step 3 ##################################
#Read in the original covariates used to fit the model
psi.cov <- read.csv('psi_cov3wk.csv')

####################### Step 4 ################################
#Extract and create a regular gradient of the covariate of interest
#Center and scale since that's what scale the betas are in from
#the model
hden <- scale(psi.cov$HDens_5km)[, 1]

#Create a regular gradient from the lowest to the highest value
#of your covariate
hd <- seq(min(hden), max(hden), length.out = 100)
summary(hd)  

####################### Step 5 ###############################
#Calculate the f parameters as a linear function of covariates 
#from the model

#Create an empty matrix with rows equal to the # iterations
#columns equal to the length of hd for each of the 
#8 unique capture histories possible (3 interacting species)
#This will hold the overall occupancy estimates
psi <- array(dim = c(8, nrow(b), length(hd)))

#Create an empty matrix with rows equal to the # iterations
#columns equal to the length of hd for each pairwise species
#interaction

#These will hold the resulting conditional occupancy estimates
#b is bobcat, c is coyote, r is red fox, g is gray fox
#1 and 0 represent presence/absence of the second species
#in the interaction
bc0 <- bc1 <- br0 <- br1 <-
  cb0 <- cb1 <- cr0 <- cr1 <-
  rb0 <- rb1 <- rc0 <- rc1 <-
  matrix(nrow = nrow(b), ncol = length(hd))

#Loop through each iteration and each value of hd and calculate the
#f parameters (there are 6, 3 for the single species and 
#3 parameters for the pairwise interactions)
for(i in 1:nrow(b)){
  for(j in 1:length(hd)){
    #This model had all f parameters modeled as:
    #~1+HDens5+Dist5km
    
    #Since the natural parameters (f) can be defined as linear functions
    #of covariates, we can calculate the f parameters
    #As the product of the relevant columns of the covariates we used 
    #in our model and the beta coefficients estimated by the model
    
    #We define the values of the covariate of interest using the continuous
    #vector of values defined above
    #and hold everything else constant at their mean 
    #(i.e. at 0 for centered and scaled continuous covariates)
    f <- numeric(6)
    f[1] <- crossprod(c(1, 0, hd[j]), b[i, 1, 1:3])[1, 1]  # f1
    f[2] <- crossprod(c(1, 0, hd[j]), b[i, 2, 1:3])[1, 1]  # f2
    f[3] <- crossprod(c(1, 0, hd[j]), b[i, 3, 1:3])[1, 1]  # f3
    f[4] <- crossprod(c(1, 0, hd[j]), b[i, 4, 1:3])[1, 1]  # f12
    f[5] <- crossprod(c(1, 0, hd[j]), b[i, 5, 1:3])[1, 1]  # f13
    f[6] <- crossprod(c(1, 0, hd[j]), b[i, 6, 1:3])[1, 1]  # f23
    
    ######################## Step 6 ###############################
    #Calculate the probabilities of each unique capture history
    #This is the numerator necessary to calculate occupancy 
    #probabilities
    n <- numeric(6)
    n[1] <- exp(sum(f))             # 111
    n[2] <- exp(sum(f[c(1, 2, 4)])) # 110
    n[3] <- exp(sum(f[c(1, 3, 5)])) # 101
    n[4] <- exp(f[1])               # 100
    n[5] <- exp(sum(f[c(2, 3, 6)])) # 011
    n[6] <- exp(f[2])               # 010
    n[7] <- exp(f[3])               # 001
    n[8] <- 1                       # 000
    
######################### Step 7 ###############################
    #Calculate Psi: Divide the probability of each capture 
    #history by the sum of capture history probabilities
    psi[, i, j] <- n / sum(n)

######################### Step 8 ###############################
    #Calculate the conditional probabilities for each species
    #Conditional on the presence/absence of an interacting species
    #marginalizing over all other species.
    
    #Conditional on absence: We sum the Psi's for only those 
    #capture histories where a given species was present in the absence
    #of the interacting species and divide by the sum of the 
    #Psi's where that given species was present without the interacting
    #species and where both species were absent
    
    #Conditional on presence: We sum the Psi's for only those 
    #capture histories where a given species was present in the presence
    #of the interacting species and divide by the sum of the 
    #Psi's where that given species was present WITH the interacting
    #species and where that species was absent but the interacing 
    #species was present 
    bc0[i, j] <- sum(psi[c(3, 4), i, j]) /
      sum(psi[c(3, 4, 7, 8), i, j])
    bc1[i, j] <- sum(psi[c(1, 2), i, j]) /
      sum(psi[c(1, 2, 5, 6), i, j])
    
    br0[i, j] <- sum(psi[c(2, 4), i, j]) /
      sum(psi[c(2, 4,6, 8), i, j])
    br1[i, j] <- sum(psi[c(1, 3), i, j]) /
      sum(psi[c(1, 3, 5, 7), i, j])
    
    cb0[i, j] <- sum(psi[c(5, 6), i, j]) /
      sum(psi[c(5, 6, 7, 8), i, j])
    cb1[i, j] <- sum(psi[c(1, 2), i, j]) /
      sum(psi[c(1, 2, 3, 4), i, j])
    
    cr0[i, j] <- sum(psi[c(2, 6), i, j]) /
      sum(psi[c(2, 6, 4, 8), i, j])
    cr1[i, j] <- sum(psi[c(1, 5), i, j]) /
      sum(psi[c(1, 5, 3, 7), i, j])
    
    rb0[i, j] <- sum(psi[c(5, 7), i, j]) /
      sum(psi[c(5, 7, 6, 8), i, j])
    rb1[i, j] <- sum(psi[c(1, 3), i, j]) /
      sum(psi[c(1, 3, 2, 4), i, j])
    
    rc0[i, j] <- sum(psi[c(3, 7), i, j]) /
      sum(psi[c(3, 7, 4, 8), i, j])
    rc1[i, j] <- sum(psi[c(1, 5), i, j]) /
      sum(psi[c(1, 5, 2, 6), i, j])
    
  }
}

########################## Step 9 #############################
#Make a dataframe of results for graphing

#x is unscaled/unceneterd housing density
#y is the conditional occupancy of each species in the presence
#and absence of each interacting species
#averaged over the MCMC iterations
#ymin is the lower bound of the 95%CI of each pairwise interaction averaged over the MCMC iterations
#ymax is the upper bound of the 95%CI  of each pairwise interaction averaged over the MCMC iterations
#We than add a column labeling the focal species, the interacting
#species and whether or not the interacting species is present
#or absent
cond_hdens <- data.frame(
  x = rep(seq(min(psi.cov$HDens_5km),max(psi.cov$HDens_5km), length.out = 100), 12),
  y = c(apply(bc0, 2, mean), apply(bc1, 2, mean),
        apply(br0, 2, mean), apply(br1, 2, mean),
        apply(cb0, 2, mean), apply(cb1, 2, mean),
        apply(cr0, 2, mean), apply(cr1, 2, mean),
        apply(rb0, 2, mean), apply(rb1, 2, mean), 
        apply(rc0, 2, mean), apply(rc1, 2, mean)), 
  ymin = c(apply(bc0, 2, quantile, probs = 0.025),
           apply(bc1, 2, quantile, probs = 0.025),
           apply(br0, 2, quantile, probs = 0.025),
           apply(br1, 2, quantile, probs = 0.025),
           apply(cb0, 2, quantile, probs = 0.025),
           apply(cb1, 2, quantile, probs = 0.025),
           apply(cr0, 2, quantile, probs = 0.025),
           apply(cr1, 2, quantile, probs = 0.025),
           apply(rb0, 2, quantile, probs = 0.025),
           apply(rb1, 2, quantile, probs = 0.025),
           apply(rc0, 2, quantile, probs = 0.025),
           apply(rc1, 2, quantile, probs = 0.025)),
  ymax = c(apply(bc0, 2, quantile, probs = 0.975),
           apply(bc1, 2, quantile, probs = 0.975),
           apply(br0, 2, quantile, probs = 0.975),
           apply(br1, 2, quantile, probs = 0.975),
           apply(cb0, 2, quantile, probs = 0.975),
           apply(cb1, 2, quantile, probs = 0.975),
           apply(cr0, 2, quantile, probs = 0.975),
           apply(cr1, 2, quantile, probs = 0.975),
           apply(rb0, 2, quantile, probs = 0.975),
           apply(rb1, 2, quantile, probs = 0.975),
           apply(rc0, 2, quantile, probs = 0.975),
           apply(rc1, 2, quantile, probs = 0.975)),
  species = rep(c('Bobcat\noccupancy',
                  'Coyote\noccupancy',
                  'Red fox\noccupancy'), each = 400),
  con = rep(c('Conditional on\ncoyote',
              'Conditional on\nred fox', 'Conditional on\nbobcat',
              'Conditional on\nred fox',
              'Conditional on\nbobcat', 'Conditional on\ncoyote'),
            each = 200),
  pa = rep(c('Absent', 'Present'), each = 100, times = 6)
)

########################## Step 10 ############################################
#Graph

#ggplot 2 takes 2 main arguments, the dataset and the aesthetics
#aesthetics define the x and y components, and a min and max for the credible intervals
condplot <- ggplot(cond_hdens, aes(x, y, ymin = ymin, ymax = ymax)) +
  #Then we splot the plot by focal species (columns) AND by the
  #interacting species (rows)
  facet_grid(con ~ species) +
  #and add credible interval shading and choose a pretty color
  geom_ribbon(aes(fill = pa), alpha = 0.3) +
  #and add a mean line which is dotted when the interacting species
  #is present and solid when the interacting species is absent
  geom_line(aes(linetype = pa)) +
  #Choose pretty colors to further highlight when the interacting
  #species is present and absent
  scale_fill_manual(values=c("green", "orange"))+
  #Label our axes
  ylab('Occupancy probability') +
  xlab('Average housing density in 5km radius') +
  #Define our "theme" to be fairly minimal, i.e. the overall look of the plot
  #This includes things like text size and color
  theme(strip.text = element_text(size = 8),
        legend.title = element_blank()) +
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  theme(text=element_text(size=8, colour="black"))+
  theme(axis.text.x = element_text(colour="black",size=8),
        axis.text.y = element_text(colour="black",size=8))

condplot

#Save our plot
ggsave('ConditionalHD.jpg', condplot, width = 14, height = 12, units = 'cm')
