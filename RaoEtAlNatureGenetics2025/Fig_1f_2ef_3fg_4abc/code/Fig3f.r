 
########################################################################################
#
# Couturier, D.-L., Insolia, L. and Guerrier, S., 20250730
# Code related to FIGURE 3f of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig3f.rdata")


########################################################################################
library(survival)


########################################################################################
# plot
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(2.5,3.5,0,.25))
ylimw = c(0,125)
xlimw = c(0,80)+c(-.5,0)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
legend(2,30,ncol=1,pch=NA,lty=19,pt.lwd=4,
       col=id.group$col,
       legend=c(id.group$id),lwd=3,
       cex=1,box.lwd=NA)
med_group = rep(NA,n.group)
for(gw in 1:n.group){
    dataw = Fig3f[Fig3f$group==id.group$id[gw],]
    dataw = dataw[order(dataw$survival,decreasing = FALSE),]
    dataw_day = split(dataw,dataw$survival)
    dataw_day = dataw_day[order(as.numeric(names(dataw_day)))]
    vect.dayw = as.numeric(names(dataw_day))
    n         = nrow(dataw)
    survivalw    = length(vect.dayw)
    yw.new    = 100
    xw.new    = 0
    medw      = FALSE
    for(dw in 1:survivalw){
        yw.old    = yw.new
        xw.old    = xw.new
        xw.new    = vect.dayw[dw]
        yw.new = yw.old - (nrow(dataw_day[[dw]])/n)*100
        if(yw.new<50&medw==FALSE){
            medw = TRUE
            med_group[gw] = xw.new
            }
        }    
    }
axis(1) 
axis(1,mean(xlimw),"Time (in days)",tick=FALSE,padj=1.2)
axis(2)  
axis(2,mean(c(0,100)),"% Survival",tick=FALSE,padj=-2)
for(gw in 1:n.group){
    axis(3,med_group[gw],font=2,font.axis=2,
         pos=97.5,cex.axis=1,lwd=2,tick=FALSE,col.axis=paste0(id.group$col[gw]))         
}
for(gw in 1:n.group){
    dataw = Fig3f[Fig3f$group==id.group$id[gw],]
    dataw = dataw[order(dataw$survival,decreasing = FALSE),]
    dataw_day = split(dataw,dataw$survival)
    dataw_day = dataw_day[order(as.numeric(names(dataw_day)))]
    vect.dayw = as.numeric(names(dataw_day))
    n         = nrow(dataw)
    survivalw    = length(vect.dayw)
    yw.new    = 100
    xw.new    = 0
    medw      = FALSE
    for(dw in 1:survivalw){
        yw.old    = yw.new
        xw.old    = xw.new
        xw.new    = vect.dayw[dw]
        segments(xw.old,yw.old,xw.new,yw.old,col=paste0(id.group$col[gw]),lwd=5)
        yw.new = yw.old - (nrow(dataw_day[[dw]])/n)*100
        segments(xw.new,yw.old,xw.new,yw.new,col=paste0(id.group$col[gw]),lwd=5)
        if(yw.new<50&medw==FALSE){
            medw = TRUE
            segments(xw.new,0,xw.new,yw.new,col=paste0(id.group$col[gw]),lwd=0.75,lty=3)
            segments(xw.new,100,xw.new,yw.new,col=paste0(id.group$col[gw]),lwd=0.75,lty=3)
            med_group[gw] = xw.new
            }
        }    
    }
abline(h=50,lty=2)


########################################################################################
# p-values
id.comp$pval.welch = NA
for(cw in 1:n.comp){# cw=1
    xw = Fig3f[Fig3f$group==id.group$id[id.comp$grp1[cw]],"survival"]
    yw = Fig3f[Fig3f$group==id.group$id[id.comp$grp2[cw]],"survival"]
    dataw = data.frame(y=c(xw,yw),group=rep(c(0,1),c(length(xw),length(yw))))
    id.comp$pval.welch[cw] = t.test(log(xw),log(yw),alternative="less")$p.value
    }
id.comp$adjpval.welch = p.adjust(id.comp$pval.welch,method="bonferroni")
print(id.comp)



