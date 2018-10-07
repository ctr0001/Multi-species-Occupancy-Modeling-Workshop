###################################################################
################## Plotting Marginal Probabilities ################
################## along an environmental gradient ################
################## from MCMC output                ################
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
#columns equal to the length of hd for each species
#These will hold the resulting marginal occupancy estimates
bob <- coy <-  red <- matrix(nrow = nrow(b), 
                                   ncol = length(hd))

#Loop through each iteration and each value of hd and calculate the
#f parameters (there are 6, three for the single species and 
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
    #Calculate the marginal Psi probabilities for 
    #each species by summing the Psi's for only those capture 
    #histories where that species was present
    bob[i, j] <- sum(psi[1:4, i, j])
    coy[i, j] <- sum(psi[c(1:2, 5:6), i, j])
    red[i, j] <- sum(psi[c(1, 3, 5, 7), i, j])
    
  }
}

head(bob)

########################## Step 9 #############################
#Make a dataframe of results for graphing

#x is unscaled/unceneterd housing density
#y is the marginal occupancy of each species averaged over the MCMC iterations
#ymin is the lower bound of the 95%CI of each species averaged over the MCMC iterations
#ymax is the upper bound of the 95%CI  of each species averaged over the MCMC iterations
#We than add a column labeling the species
hdens <- data.frame(
  x = rep(seq(min(psi.cov$HDens_5km),max(psi.cov$HDens_5km), 
              length.out = 100), 3),
  y = c(apply(bob, 2, mean), apply(coy, 2, mean),
        apply(red, 2, mean)),
  ymin = c(apply(bob, 2, quantile, probs = 0.025),
           apply(coy, 2, quantile, probs = 0.025),
           apply(red, 2, quantile, probs = 0.025)),
  ymax = c(apply(bob, 2, quantile, probs = 0.975),
           apply(coy, 2, quantile, probs = 0.975),
           apply(red, 2, quantile, probs = 0.975)),
  species = rep(c('Bobcat', 'Coyote', 'Red fox'), each = 100)
)

head(hdens)

########################## Step 10 ############################################
#Graph

#ggplot 2 takes 2 main arguments, the dataset and the aesthetics
#aesthetics define the x and y components, and a min and max for the credible intervals
marplot <- ggplot(hdens, aes(x=x, y=y, ymin = ymin, ymax = ymax)) +
  #Then we split the plot by the species
  facet_grid(species~.)+
  #and add credible interval shading and choose a pretty color
  geom_ribbon(fill="purple", alpha = 0.3) +
  #and add a mean line
  geom_line() + 
  #Label our axes
  ylab('Marginal occupancy probability') + 
  xlab('Housing density\n5km radius') +
  #Define our "theme" to be fairly minimal, i.e. the overall look of the plot
  #This includes things like text size and color
  theme(strip.background = element_rect(fill = 'white'))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  theme(text=element_text(size=12, colour="black"))+
  theme(axis.text.x = element_text(colour="black",size=12),
        axis.text.y = element_text(colour="black",size=12))

#Draw the plot in the viewer
marplot

ggsave('Marginal Plot.jpg', marplot, width = 14, height = 14, units = 'cm')
