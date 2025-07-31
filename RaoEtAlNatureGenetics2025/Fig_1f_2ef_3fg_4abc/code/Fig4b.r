 
########################################################################################
#
# Couturier, D.-L., Insolia, L. and Guerrier, S., 20250730
# Code related to FIGURE 4b of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig4b.rdata")


########################################################################################
library(gamlss)


########################################################################################
# fit
fit = gamlss(n.met~group+log(n.day),data=Fig4b,family=NBI)
ar.hat.ggp = array(NA, dim=c(n.group,n.group,2),
                   dimnames=list(id.group$id,id.group$id,c("estimates","p-values")))
for(gw in 1:n.group){# gw=1
    data = Fig4b
    # levels
    levelw = c(id.group$id[gw],id.group$id[-gw])
    data$group = factor(as.character(data$group),levels=levelw)
    # fit
    fitw = summary(gamlss(n.met~group+log(n.day),data=data,family=NBI),save=TRUE)$coef.table

    ar.hat.ggp[levelw,gw,1] = fitw[1:n.group,1]
    ar.hat.ggp[levelw,gw,2] = fitw[1:n.group,4] 
    ar.hat.ggp[gw,gw,2]     = NA
}



########################################################################################
# plot
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(2.5,3.5,0,.25))
ylimw = c(0,max(Fig4b$n.met))
xlimw = c(0,120)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
legend("topleft",ncol=1,pch=NA,lty=19,pt.lwd=4,
       col=id.group$col,
       legend=c(id.group$id),lwd=3,
       cex=1,box.lwd=NA)   
axis(1)
axis(1,at=mean(xlimw),padj=1.5,"Number of days",tick=FALSE)
axis(2)
axis(2,at=mean(ylimw),padj=-1.5,"Number of metastases",tick=FALSE)
med_group = rep(NA,n.group)
set.seed(1)
points(Fig4b$n.day+runif(nrow(Fig4b),-.05,.05),
       Fig4b$n.met+runif(nrow(Fig4b),-.05,.05),
       col=paste0(id.group[Fig4b$group,"col"],50),pch=16)
for(gw in 1:n.group){
    axe.x = seq(0,120,length=250)
    f.x   = exp(ar.hat.ggp[gw,gw,1]+log(axe.x)*fit$mu.coefficients["log(n.day)"])
    lines(axe.x,f.x,col=id.group$col[gw])
}


########################################################################################
# p-values
print(ar.hat.ggp[,,"p-values"])


