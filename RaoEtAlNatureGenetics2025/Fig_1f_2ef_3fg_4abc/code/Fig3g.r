  
########################################################################################
#
# Couturier, D.-L., Insolia, L. and Guerrier, S., 20250730
# Code related to FIGURE 3g of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig3g.rdata")


########################################################################################
library(nlme)
library(multcomp)



########################################################################################
# fit
fit = lme(I(volume^(1/3))~0+x1+x2+x3+x4+x5,data=Fig3g,random=~day|mouse)


########################################################################################
# plot
set.seed(1)
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(3.5,3.5,0,0))
ylimw = c(0,2000)
xlimw = c(10,80)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
axis(1) 
axis(1,mean(xlimw),"Time (in days)",tick=FALSE,padj=1.2)
axis(2)     
axis(2,mean(c(0,2000)),"Volume (mm3)",tick=FALSE,padj=-2)
legend(60,500,ncol=1,pch=NA,lty=19,pt.lwd=4,# 12.5,2050
       col=id.group$col,
       legend=c(id.group$id),lwd=3,
       cex=1,box.lwd=NA)    
for(mw in 1:nrow(Fig3f)){
    dataw = Fig3g[Fig3g$mouse==Fig3f$id[mw],]
    dataw = dataw[!is.na(dataw$volume),]
    if(nrow(dataw)>0){
        colw  = paste0(id.group[dataw$group[1],"col"],50)
        lines(dataw$day,I(dataw$volume),col=colw)
        points(dataw$day,I(dataw$volume),col=colw,pch=16,cex=.5)
        }
    }
kw    = id.knot$pos
axe.x = seq(10,80,length=1000)
n.x   = length(axe.x)
fake  = data.frame(group=factor(rep(id.group$id,each=n.x),levels=id.group$id),
                   type=NA,
                   day=rep(axe.x,n.group),
                   mouse=Fig3f$id[1],
                   volume=NA)
fake$type = factor(as.numeric(as.numeric(fake$group)>2),levels=0:1,labels=c("Cas-9","HNF4G-KO"))
fake$day.shited = fake$day-id.knot$knot1[kw] 
X1 = fake$day.shited*model.matrix(~group,data=fake)
X2 = fake$day.shited*model.matrix(~type,data=fake)[,2]
X3 = X1[,3:4]-id.knot$knot2[kw]
X3[X3<0] = 0
X = cbind(X1[,1:2],X2,X3)#,X1[,1]^2)
colnames(X) = paste0("x",1:ncol(X))
fake$volume = X%*%coef(summary(fit))[,1]
fake$volume[fake$day<id.knot$knot1[kw]] = 0
fake = fake[fake$volume<(1750^(1/3)),]
for(gw in 1:n.group){# gw=1
        posw = fake$group==id.group$id[gw]
        lines(fake$day[posw],fake$volume[posw]^3,
              col=id.group$col[gw],lwd=5)
        }


########################################################################################
# p-values      
glhtw<-summary(glht(fit,linfct=K[,1:5]))
print(glhtw)


