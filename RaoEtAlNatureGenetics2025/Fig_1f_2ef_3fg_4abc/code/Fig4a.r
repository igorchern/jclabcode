  
########################################################################################
#
# Couturier, D.-L., Insolia, L. and Guerrier, S., 20250730
# Code related to FIGURE 4a of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig4a.rdata")


########################################################################################
library(nlme)
library(multcomp)


########################################################################################
# fit
dataw     = Fig4a
for(gw in 1:n.group){
    dataw = cbind(dataw,as.numeric(dataw$group==id.group$id[gw])*as.numeric(dataw$day))
    colnames(dataw)[ncol(dataw)] = paste0("x",gw)
}
formula.full = as.formula(paste0("log1p(volume)~0+",paste0("x",c(1:n.group),collapse="+")))
fit.lme      = lme(log1p(volume) ~ 0 + x1 + x2 + x3 + x4,data=dataw,random=~day+0|mouse)
coef.lme     = coef(summary(fit.lme))    


########################################################################################
# plot
set.seed(1)
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(3.5,3.5,0,0))
ylimw = c(0,2000)
xlimw = c(0,50)+c(-.5,.5)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
axis(1) 
axis(1,mean(xlimw),"Time (in days)",tick=FALSE,padj=1.2)
axis(2)     
axis(2,mean(c(0,2000)),"Volume (mm3)",tick=FALSE,padj=-2)
legend("topleft",ncol=1,pch=NA,lty=19,pt.lwd=4,
       col=id.group$col,
       legend=c(id.group$id),lwd=3,
       cex=1,box.lwd=NA)   
for(mw in 1:n.mouse){
    dataw = Fig4a[Fig4a$mouse==id.mouse$id[mw],]
    dataw = dataw[!is.na(dataw$volume),]
    if(nrow(dataw)>0){
        colw  = paste0(id.group[dataw$group[1],"col"],50)
        lines(dataw$day,I(dataw$volume),col=colw)
        points(dataw$day,I(dataw$volume),col=colw,pch=16,cex=.5)
        }
    }
for(gw in 1:n.group){
    axe.x = seq(0,50,length=250)
    betaw = coef.lme[gw,"Value"]
    f.x   = betaw*axe.x
    lines(axe.x,exp(f.x)+1,col=id.group$col[gw],lwd=2)
    }


########################################################################################
# p-values     
K = cbind(-1,diag(n.group)[-1,-1])
glhtw<-summary(glht(fit.lme,linfct=K))
print(glhtw)


