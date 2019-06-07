
setwd("/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/Scripts/2_Meta_RT_Mem")
# OPTIONAL: Clear R's memory and graphics:
rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off() # This closes all of R's graphics windows.
#Note:This is what it looks like on my MacBook 
library(metafor)# Get the functions loaded into R's working memory:


dat <- read.csv("1_Meta_RT_Mem.csv",header=TRUE);

res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="FE") ### use r-to-z transformed correlations
funnel(res, yaxis="sei")
forest(res, slab=paste("NIH",dat$SubjID,"(",dat$Trials,")"),xlim=c(-1.8,1.8),alim=c(-1,1), ylim=c(-1.8, 22))  # ylim=c(-4.5, 12) for smaller table; ylim=c(-2, 21.8) for large table
res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML") ### use r-to-z transformed correlations
addpoly(res,cex=1)

text(-0.6,  20.8, "StudyID (#of Intrusion Trials)",  pos=2)
text(1.6,  20.8, "Zr [95% CI]", pos=2)

dat <- read.csv("Meta_Intrusionwithnoise.csv",header=TRUE);

resf <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="FE") ### use r-to-z transformed correlations
funnel(resf, yaxis="sei")
forest(resf, slab=paste("NIH",dat$SubjID,"(",dat$Trials,")"),xlim=c(-1.8,1.8),alim=c(-1,1), ylim=c(-1.5, 28))  # ylim=c(-4.5, 12) for smaller table; ylim=c(-2, 21.8) for large table
res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="REML") ### use r-to-z transformed correlations
addpoly(res,cex=1)

text(-0.55,  27, "StudyID (#of Response Trials)",  pos=2)
text(1.6,  27, "Zr [95% CI]", pos=2)

## ==============
dat <- read.csv("Meta_RT_Mem.csv",header=TRUE);

res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="FE") ### use r-to-z transformed correlations
funnel(res, yaxis="sei")
forest(res, slab=paste("NIH",dat$SubjID,"(",dat$Trials,")"),xlim=c(-1.8,1.8),alim=c(-1,1), ylim=c(-1.5, 28))  # ylim=c(-4.5, 12) for smaller table; ylim=c(-2, 21.8) for large table
res2 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML") ### use r-to-z transformed correlations
addpoly(res2,cex=1)

text(-0.55,  27, "StudyID (#of Response Trials)",  pos=2)
text(1.6,  27, "Zr [95% CI]", pos=2)

taf <- trimfill(res)
funnel(taf)
funnel(res, main="Standard Error")

# with noise
dat <- read.csv("Meta_RT_Memwithnoise.csv",header=TRUE);
resf <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="FE") ### use r-to-z transformed correlations
funnel(resf, yaxis="sei")
forest(resf, slab=paste("NIH",dat$SubjID,"(",dat$Trials,")"),xlim=c(-1.8,1.8),alim=c(-1,1), ylim=c(-1.5, 31))  # ylim=c(-4.5, 12) for smaller table; ylim=c(-2, 21.8) for large table
res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="REML") ### use r-to-z transformed correlations
addpoly(res,cex=1)

text(-0.55,  30, "StudyID (#of Response Trials)",  pos=2)
text(1.6,  30, "Zr [95% CI]", pos=2)

taf <- trimfill(res)
funnel(taf)
funnel(res, main="Standard Error")


## test of moderation
res0 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="REML", data=dat)
res1 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="REML", mods = ~ ni + ACC, data=dat)
anova(res0,res1)


#=======================



## test of moderation
res0 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="REML", data=dat)
res1 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="REML", mods = ~ ACC, data=dat)



# For Pos_Neu task by stimuli
predict(res01,newmods=cbind(c(0, 0, 0 ,0),c(1,-1,1,-1),c(1 ,1, -1, -1),c(0.5,0.5,0.5,0.5),c(1,-1,-1,1)) )# gender by stim, keep gender constant


# For Neg_Neu domain by stimuli
predict(res06,newmods=cbind(c(1,1,-1,-1),c(0,0,0,0),c(1,-1,1,-1),c(0.5,0.5,0.5,0.5), c(0,0,0,0), c(1, -1, -1,1)*0.5, c(0,0,0,0),c(0,0,0,0), c(1,-1,-1,1),  c(1,1,-1,-1)*0.5)) 
                            # domain task stimuli gender # taskbystimuli + genderbystimuli + taskbygender # domainbytask + domainbystimuli +domainbygender
                             

# For Neg_Neu domain by task
predict(res06,newmods=cbind(c(1,1,-1,-1),c(1,-1,1,-1),c(0,0,0,0),c(0.5,0.5,0.5,0.5), c(0,0,0,0), c(0,0,0,0),c(1,-1,1,-1)*0.5, c(1,-1,-1,1),c(0,0,0,0),  c(1,1,-1,-1)*0.5)) 
# domain task stimuli gender # taskbystimuli + genderbystimuli + taskbygender # domainbytask + domainbystimuli +domainbygender



# For overall effect domain by gender by stimuli
predict(res06,newmods=cbind(c(0,0,0,0),c(0,0,0,0),c(1,-1,1,-1),c(1,1,0,0),
                            c(0,0,0,0), c(1,-1,0,0),c(0,0,0,0),
                            c(0,0,0,0), c(0,0,0,0), c(0,0,0,0)))# gender by stimuli, 



resall <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML", mods = ~ domain+ task + stimuli + gender + 
                domainbytask + domainbystimuli+ domainbygender + taskbystimuli + taskbygender + genderbystimuli   ,
              slab=paste(author, year, sep=", "), data=dat)


x <- predict(res3,newmods=cbind(c(1,1,-1,-1),c(1,1,1,1),c(1,0,1,0),c(1,1,-1,-1),c(1,0,-1,0)))# task by gender, keep stimu constant
predict(res4,newmods=cbind(c(1,1,-1,-1),c(1,-1,1,-1),c(1 ,1, 1, 1),c(1,-1,-1,1))) # task by stim, keep gender constant


# for effect 4

res0 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML", mods = ~ matcharousal+domain + task + stimuli + gender ,
            slab=paste(author, year, sep=", "), data=dat)

res01 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML", mods = ~ matcharousal+domain + task + stimuli + gender + 
               taskbystimuli,
             slab=paste(author, year, sep=", "), data=dat)

res02 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML", mods = ~matcharousal+ domain + task + stimuli + gender + 
               taskbystimuli+genderbystimuli ,
             slab=paste(author, year, sep=", "), data=dat)

res03 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML", mods = ~ matcharousal+domain + task + stimuli + gender + 
               taskbystimuli + genderbystimuli  + taskbygender ,
             slab=paste(author, year, sep=", "), data=dat)

res06 <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="ML", mods = ~ matcharousal + domain + task + stimuli + gender + 
               domainbytask + domainbystimuli+ domainbygender+ taskbystimuli + taskbygender + genderbystimuli,
             slab=paste(author, year, sep=", "), data=dat)


### Conservative effects

dat <- read.csv("negneu_conservative.csv",sep=",",header=TRUE);
dat <- read.csv("posneu_conservative.csv",sep=",",header=TRUE);
dat <- read.csv("negpos_conservative.csv",sep=",",header=TRUE);

res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="FE") ### use r-to-z transformed correlations
funnel(res, yaxis="sei")
forest(res, slab=paste(dat$author, dat$year, sep=", "),
       xlim=c(-3,2),alim=c(-1,1), ylim=c(-2, 22) )  # ylim=c(-4.5, 12) for smaller table; ylim=c(-2, 21.8) for large table
res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="HE") ### use r-to-z transformed correlations
addpoly(res)
text(-3,                21, "Study ID. Author(s) and Year",  pos=4)
text(1.8,                21, "Zr[95% CI]", pos=2)





