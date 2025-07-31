 
########################################################################################
#
# Couturier, D.-L., Insolia, L. and Guerrier, S., 20250730
# Code related to FIGURE 4c of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig4c.rdata")


########################################################################################
library(gamlss)
library(gamlss.cens)
gen.cens(LOGNO,type="right")


########################################################################################
# fit
ar.hat.ggp = array(NA, dim=c(n.group,n.group,2),
                   dimnames=list(id.group$id,id.group$id,c("estimates","p-values")))
for(gw in 1:n.group){# gw=1
    data = Fig4c
    # levels
    levelw = c(id.group$id[gw],id.group$id[-gw])
    data$group = factor(as.character(data$group),levels=levelw)
    # fit
    fitw  = summary(gamlss(n.day_surv~group,sigma.formula=~group,
                    data=data,family=LOGNOrc),save=TRUE)$coef.table
    ar.hat.ggp[levelw,gw,1] = fitw[1:n.group,1]
    ar.hat.ggp[levelw,gw,2] = fitw[1:n.group,4] 
    ar.hat.ggp[gw,gw,2]     = NA
}


########################################################################################
# plot
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(3.5,3.5,0,.25))
ylimw = c(0,1)
xlimw = c(0,100)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
legend("bottomleft",ncol=1,pch=NA,lty=19,pt.lwd=4,
       col=id.group$col,
       legend=c(id.group$id),lwd=3,
       cex=1,box.lwd=NA)   
axis(1)
axis(1,at=mean(xlimw),padj=1.5,"Number of days",tick=FALSE)
axis(2)
axis(2,at=mean(ylimw),padj=-1.5,"Survival probability",tick=FALSE)
med_group = rep(NA,n.group)
set.seed(1)
sfitw  = survfit(n.day_surv~group,data=Fig4c)
groupw = rep(id.group$id,sfitw[[10]])
for(gw in 1:n.group){
    sprobw = c(1,sfitw[[6]][groupw==id.group$id[gw]])
    timew  = c(0,sfitw[[2]][groupw==id.group$id[gw]])
    lines(timew,sprobw,type="s",col=paste0(id.group$col[gw],99),lwd=2)
    censw  = sfitw[[5]][groupw==id.group$id[gw]]
    x.cens = rep(timew[-1][censw>0],censw[censw>0])
    y.cens = rep(sprobw[-1][censw>0],censw[censw>0])
    x.cens = x.cens+runif(length(x.cens),-.5,.5)
    points(x.cens,y.cens,pch=1,col=paste0(id.group$col[gw],99))
    }
abline(h=0.5,lty=2)


########################################################################################
# p-values
print(ar.hat.ggp[,,"p-values"])


