########################################################################
######################### Make a MARK unput file #######################
################# from raw detection/nondetection data #################
########################################################################

#Set your working directory
setwd("C:/Users/Arielle/Desktop/NRC/Biodiversity Lab/Daily Planet Talks/TWS/2018/Multi-species Workshop/Inp")

#Read in the detection/non detection data files for each species
#In this example there are 3 species, two jackals and a genet
#The files are one row for each site, one column for each occasion
#and a 1 or 0 for whether that species was detected at that site on 
#that occasion
sp1 <- read.csv('Bobcat_3wk.csv', header = F)
sp2 <- read.csv('Coyote_3wk.csv', header = F)
sp3 <- read.csv('RedFox_3wk.csv', header = F)
head(sp3)

#Read in the p and psi covariates
p <- read.csv('detection data.csv')
psi <- read.csv('psi_cov3wk.csv')
head(p)
head(psi)

#Create an empty matrix to hold the capture histories
ch <- matrix(nrow = nrow(sp1), ncol = ncol(sp1))

#Loop through each site and occasion, adding ".." to the
#capture history file if the site was not
#sampled on that occassion, '00' if none of the species were detected
#at that site on that occasion, '01' if only species1, '02' if only
#species2, '03' if species1 and species2, '04' if only species3, '05'
#if species1 and species3, '06' if species2 and species3, '07' if all
#species were detected at that site on that occasion.
for(i in 1:nrow(sp1)){
  for(j in 1:ncol(sp1)){
    if(is.na(sp1[i, j]) | is.na(sp2[i, j]) | is.na(sp3[i, j])){
      ch[i, j] <- '..'
    } else {
      if(sp1[i, j] == 0 & sp2[i, j] == 0 & sp3[i, j] == 0){
        ch[i, j] <- '00'
      } 
      if(sp1[i, j] == 1 & sp2[i, j] == 0 & sp3[i, j] == 0){
        ch[i, j] <- '01'
      }
      if(sp1[i, j] == 0 & sp2[i, j] == 1 & sp3[i, j] == 0){
        ch[i, j] <- '02'
      }
      if(sp1[i, j] == 1 & sp2[i, j] == 1 & sp3[i, j] == 0){
        ch[i, j] <- '03'
      }
      if(sp1[i, j] == 0 & sp2[i, j] == 0 & sp3[i, j] == 1){
        ch[i, j] <- '04'
      }
      if(sp1[i, j] == 1 & sp2[i, j] == 0 & sp3[i, j] == 1){
        ch[i, j] <- '05'
      }
      if(sp1[i, j] == 0 & sp2[i, j] == 1 & sp3[i, j] == 1){
        ch[i, j] <- '06'
      }
      if(sp1[i, j] == 1 & sp2[i, j] == 1 & sp3[i, j] == 1){
        ch[i, j] <- '07'
      }
    }
  }
}

head(ch)

#Concatenate across the occasions for each site 
ch_full <- paste0(ch[, 1], ch[, 2], ch[, 3])
head(ch_full)

#Make a dataframe for the output, adding the covariates after the 
#concatenated capture histories for each site
out <- data.frame(
  ch = ch_full,
  group = rep(1, nrow(sp1)),
  Temp1 = scale(p$Temp1)[, 1],
  Temp2 = scale(p$Temp2)[, 1],
  Temp3 = scale(p$Temp3)[, 1],
  Dist_5km = scale(psi$Dist_5km)[, 1],
  HDens_5km = scale(psi$HDens_5km)[, 1]
  )

head(out)

#Save the fild as a .inp which can be read directly into MARK
write.table(out, 'multi_species.inp', eol = ";\n",
            row.names = F, col.names = F, quote = F)
