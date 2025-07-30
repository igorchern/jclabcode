 
########################################################################################
#
# Couturier, D.-L., Insolia, I. and Guerrier S. [20250730]
# Code related to FIGURE 2e of Rao S. et al (2025), 
# "Transcription Factor Switching Drives Subtype Specific Pancreatic Cancer"
# 
########################################################################################


rm(list=ls())


########################################################################################
load("data/Fig2e.rdata")


########################################################################################
library(survival)


########################################################################################
# plot
par(mfrow=c(1,1),omi=c(.1,.1,0,0),mar=c(2.5,3.5,0,.25))
ylimw = c(0,100)
xlimw = c(0,71)+c(-.5,.5)
plot(1,1,pch="",axes=FALSE,
     xlab="",ylab="",
     main="",col.main="red",cex.main=1.5,
     ylim=ylimw,xlim=xlimw)   
legend(1,18.5,ncol=1,pch=NA,lty=19,pt.lwd=4,
       col=id.group$col,
       legend=c(id.group$name),lwd=3,
       cex=1,box.lwd=NA)   
axis(1) 
axis(1,mean(xlimw),"Time (in days)",tick=FALSE,padj=1.2)
axis(2)     
axis(2,mean(ylimw),"% Survival",tick=FALSE,padj=-2)
med_group = rep(NA,n.group)
for(gw in 1:n.group){
    dataw = Fig2e[Fig2e$group==id.group$id[gw],]
    dataw = dataw[order(dataw$n.day,decreasing = FALSE),]
    dataw_day = split(dataw,dataw$n.day)
    dataw_day = dataw_day[order(as.numeric(names(dataw_day)))]
    vect.dayw = as.numeric(names(dataw_day))
    n         = nrow(dataw)
    n.dayw    = length(vect.dayw)
    yw.new    = 100
    xw.new    = 0
    medw      = FALSE
    for(dw in 1:n.dayw){
        yw.old    = yw.new
        xw.old    = xw.new
        xw.new    = vect.dayw[dw]
        segments(xw.old,yw.old,xw.new,yw.old,col=id.group$col[gw],lwd=5)
        yw.new = yw.old - (nrow(dataw_day[[dw]])/n)*100
        segments(xw.new,yw.old,xw.new,yw.new,col=id.group$col[gw],lwd=5)
        if(yw.new<50&medw==FALSE){
            medw = TRUE
            segments(xw.new,0,xw.new,yw.new,col=id.group$col[gw],lwd=0.75,lty=3)
            med_group[gw] = xw.new
            }
        }
    }
abline(h=50,lty=2)


########################################################################################
# p-value
surv_diff = survdiff(Surv(n.day) ~ group, data = Fig2e)
print(surv_diff)



