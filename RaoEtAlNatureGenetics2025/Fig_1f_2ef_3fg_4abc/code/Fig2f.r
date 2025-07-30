 
########################################################################################
#
# Couturier, D.-L., Insolia, I. and Guerrier S. [20250730]
# Code related to FIGURE 2f of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig2f.rdata")


########################################################################################
library(nlme)



########################################################################################
# fit
fit1 = lme(I(volume^(1/3))~group*day,data=Fig2f,random=~1|mouse)


########################################################################################
# plot
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(2.5,3.5,0,.25))
ylimw = c(0,(max(Fig2f$volume))*1.1)
xlimw = c(0,max(Fig2f$day))+c(-.5,.5)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
legend(1,1500,ncol=1,pch=NA,lty=19,pt.lwd=4,
       col=id.group$col,
       legend=c(id.group$name),lwd=3,
       cex=1,box.lwd=NA)    
axis(1) 
axis(1,mean(xlimw),"Time (in days)",tick=FALSE,padj=1.2)
axis(2)     
axis(2,mean(ylimw),"Tumour Volume (mm3)",tick=FALSE,padj=-3.2)
for(mw in 1:nrow(Fig2e)){
    dataw = Fig2f[Fig2f$mouse==Fig2e$id[mw],]
    dataw = dataw[!is.na(dataw$volume),]
    if(nrow(dataw)>0){
        colw  = ifelse(dataw$group[1]=="Cas-9 control",
                       gray(.25),id.group[dataw$group[1],"col"])
        id.group[dataw$group[1],"col"]#.p(id.group[dataw$group[1],"col"],50)
        lines(dataw$day,I(dataw$volume),col=colw,lwd=.75)
    }
}
axe.x = seq(0,100,length=2500)
n.x   = length(axe.x)
fake  = data.frame(group=factor(rep(id.group$id,each=n.x),levels=id.group$id),
                   day=rep(axe.x,n.group),
                   mouse=Fig2e$id[1])
X1 = model.matrix(~group*day,data=fake)    
pred1 = X1%*%coef(summary(fit1))[,1]
which.group = rep(NA,n.group)
for(gw in 1:n.group){#gw=1
    which.group[gw] = min(which(((pred1[fake$group==id.group$id[gw]])^3)>1500))
    predw = ((pred1[fake$group==id.group$id[gw]])^3)
    predw[predw<0] = 0
    lines(axe.x[1:which.group[gw]],predw[1:which.group[gw]],
          col=id.group$col[gw],lwd=5)
}



########################################################################################
# p-value
fit0 = lme(I(volume^(1/3))~day,data=Fig2f,random=~1|mouse, )
anova(update(fit0,method="ML"),update(fit1,method="ML"))
