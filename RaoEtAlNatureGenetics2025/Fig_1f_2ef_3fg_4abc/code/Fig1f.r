 
########################################################################################
#
# Couturier, D.-L., Insolia, L. and Guerrier, S., 20250730
# Code related to FIGURE 1f of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig1f.rdata")


########################################################################################
library(survival)


########################################################################################
# plot
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(2.5,2.75,1.25,.25))
ylimw = c(0,100)
xlimw = c(0,90)+c(-.5,.5)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
legend(60,100,ncol=1,pch=c(15,15,4),pt.cex=c(2,2,1),lty=NA,
       col=c(id.group$col,"1"),
       legend=c(id.group$id,"Censored participants"),lwd=3,
       cex=1,box.lwd=NA)    
med_group = rep(NA,n.group)
for(gw in 1:n.group){
    Fig1fw = Fig1f[Fig1f$group==id.group$id[gw],]
    Fig1fw = Fig1fw[Fig1fw$status=="deceased",]
    Fig1fw = Fig1fw[order(Fig1fw$pi,decreasing = TRUE),]
    whichw = min(which(Fig1fw$pi<=0.5))
    med_group[gw] = Fig1fw[whichw,"survival"]
    }
axis(1,c(seq(0,90,10),100)) 
axis(1,mean(xlimw),"Time (in months)",tick=FALSE,padj=1.2)
axis(2,seq(0,100,10),cex.axis=.75)     
axis(2,mean(c(0,100)),"% Progression free",tick=FALSE,padj=-2)
for(gw in 1:n.group){
    axis(3,med_group[gw],labels=round(med_group[gw],1),font=2,font.axis=2,
         pos=97.5,cex.axis=.75,lwd=2,tick=FALSE,col.axis=paste0(id.group$col[gw]))         
}
for(gw in 1:n.group){
    Fig1fw = Fig1f[Fig1f$group==id.group$id[gw],]
    Fig1fw1 = Fig1fw
    Fig1fw1 = Fig1fw1[order(Fig1fw1$survival,decreasing = FALSE),]
    yw.new    = 100
    xw.new    = 0
    medw      = FALSE
    for(dw in 1:nrow(Fig1fw1)){
        yw.old    = yw.new
        xw.old    = xw.new
        xw.new    = Fig1fw1[dw,"survival"]
        segments(xw.old,yw.old,xw.new,yw.old,col=paste0(id.group$col[gw]),lwd=3.5)
        yw.new = Fig1fw1[dw,"pi"]*100
        segments(xw.new,yw.old,xw.new,yw.new,col=paste0(id.group$col[gw]),lwd=3.5)
        if(yw.new<50&medw==FALSE){
            medw = TRUE
            segments(xw.new,0,xw.new,yw.new,col=paste0(id.group$col[gw]),lwd=0.75,lty=3)
            segments(xw.new,100,xw.new,yw.new,col=paste0(id.group$col[gw]),lwd=0.75,lty=3)
            med_group[gw] = xw.new
            }
        }
    Fig1fw2 = Fig1fw[Fig1fw$status!="deceased",]
    points(Fig1fw2$survival,Fig1fw2$pi*100,col="black",pch=4)
    }
abline(h=50,lty=2)


########################################################################################
# p-value
surv_diff <- survdiff(Surv(survival, I(status=="deceased")) ~ group, data = Fig1f)
print(surv_diff)
