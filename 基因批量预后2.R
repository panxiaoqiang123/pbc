#futime ?????????????????????
#####??????????????????
rm(list = ls())
library(survival)
library(survminer)
library(dplyr)
pFilter=0.05                                                      #?????????????????????
setwd("D:/dingxiangRyuyan/matrix")
rt=read.table("CDC25A????????????.txt",header=T,sep="\t",check.names=F,row.names=1)     #??????????????????
rt <- as.data.frame(rt)
outTab=data.frame()
sigGenes=c("OS","OS.status")
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
rt <- rt%>% 
  filter(!is.na(OS.status)) %>%
  filter(!is.na(OS))
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(OS, OS.status) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


#???????????????
#??????????????????
rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#????????????
pdf(file="forest.pdf", width = 6,height = 4.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#????????????????????????????????????
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#???????????????
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#????????????KM?????????
for(gene in colnames(uniSigExp)[4:ncol(uniSigExp)]){
  group <- ifelse(uniSigExp[[gene]] > median(uniSigExp[[gene]]),
                  "high", "low")
  diff=survdiff(Surv(OS, OS.status) ~ group, data = uniSigExp)
  pValue=1-pchisq(diff$chisq,df=1)#survdiff?????????P?????????
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS, OS.status) ~ group, data = uniSigExp)
  surPlot=ggsurvplot(fit, 
                     data=uniSigExp,
                     conf.int=TRUE,
                     pval=pValue,
                     pval.size=5,
                     legend.labs=c("High", "Low"),
                     legend.title=gene,
                     xlab="Time (years)",
                     ylab="Overall survival",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("#d7191c", "#2b83ba"),
                     risk.table=T,
                     risk.table.height=.25)
  pdf(file=paste0(gene,".pdf"), onefile=FALSE, width=6.5, height=5.5)
  print(surPlot)
  dev.off()
}