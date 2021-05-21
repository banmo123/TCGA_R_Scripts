rm(list=ls())
setwd("E:\\oral\\DNAMethylation\\")
geneMethylation1<-read.table(file = "geneMethy1.txt",header = T)
geneMethylation2<-read.table(file = "geneMethy2.txt",header = T)
geneMethylation3<-read.table(file = "geneMethy3.txt",header = T)
geneMethylation4<-read.table(file = "geneMethy4.txt",header = T)

alldata<-data.frame()
alldata<-merge(geneMethylation1,geneMethylation2,by="id",all=TRUE)
alldata<-merge(alldata,geneMethylation3,by="id",all=TRUE)
alldata<-merge(alldata,geneMethylation4,by="id",all=TRUE)

alldata<-alldata[,c(1:5,91:98,200:210,279:289,6:90,99:199,211:278,290:381)]

normalNum=34                            #正常样品的数目
tumorNum=346                            #癌症样品的数目
fdrFilter=0.05                          #挑选差异基因的fdr阈值
logFCfilter=1                           #挑选差异基因的logFC阈值

#读取数据
library(limma)
outTab=data.frame()
grade=c(rep(1,normalNum),rep(2,tumorNum))
#Type=c(rep("Normal",normalNum),rep("Tumor",tumorNum))
rt=alldata
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)#去除基因名重复项，β值取均值
data=data[rowMeans(data)>0,]#去除表达量都为0的基因

#矫正数据
data=normalizeBetweenArrays(data)
normalData=cbind(id=row.names(data),data)
write.table(normalData,file="normalizeMethy.txt",sep="\t",row.names=F,quote=F)


#差异分析
for(i in row.names(data)){
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  
  normalGeneMeans=mean(data[i,1:normalNum])
  tumorGeneMeans=mean(data[i,(normalNum+1):ncol(data)])
  logFC=log2(tumorGeneMeans)-log2(normalGeneMeans)  
  pvalue=wilcoxTest$p.value
  normalMed=median(data[i,1:normalNum])
  tumorMed=median(data[i,(normalNum+1):ncol(data)])
  diffMed=tumorMed-normalMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,
                 cbind(gene=i,
                       normalMean=normalGeneMeans,
                       TumorMean=tumorGeneMeans,
                       logFC=logFC,
                       pValue=pvalue))}
  
}

#对p值进行矫正
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

#输出所有基因的甲基化差异情况

#write.table(outTab,file="allGene.xls",sep="\t",row.names=F,quote=F)

#输出差异甲基化的基因
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
#write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)

allsymbol<-read.table(file = "E:\\oral\\finallsymbol.txt",sep = "\t")
dgene<-merge(allsymbol,outDiff,by.x="V1",by.y="gene",all=FALSE)

#survival analysis
survgene<-alldata
rownames(survgene)<-alldata[,"id"]
survgene<-survgene[,-1]
colnames(survgene)<-substr(colnames(survgene),1,12)

barcode <- colnames(survgene)
barcode <- gsub("\\.","-",colnames(survgene))
colnames(survgene)<- barcode
library(survival)
library(survminer)
##融合数据
load(file = "clin.RData")
interbarcode <- intersect(barcode,clin$sample)
gene <- rownames(survgene)
dir.create("surfig")
target_gene<-c()
for(i in 1:nrow(survgene)){
  geneExp <- survgene[gene[i],interbarcode]
  geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame() 
  geneExp$sample <- rownames(geneExp)
  mergdata <- merge(clin,geneExp,by = "sample")
  
  mergdata$Group[mergdata[,gene[i]] <0.2]<-"NONE"
  mergdata$Group[mergdata[,gene[i]] >0.2&mergdata[,gene[i]] <0.6]<-"PART"
  mergdata$Group[mergdata[,gene[i]] >0.6]<-"ALL"
  
  if(length(unique(as.vector(mergdata$Group)))>1){
  diff<-survdiff(Surv(follow, survivalstate) ~ Group, data=mergdata)
  
  pValue<-1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue=signif(pValue,4)
    pValue=format(pValue,scientific=TRUE)
  }else{
    pValue=round(pValue,3)
  }
  if(pValue<0.05){
    target_gene<-c(target_gene,gene[i])
    
    fit <- survfit(Surv(follow, survivalstate) ~ Group, data=mergdata)
    summary(fit)
    filename<-paste(".\\surfig\\",gene[i],"-SurvivalCurve.pdf",sep="")
    pdf(file=filename,width=5.5,height = 5)
    plot(fit,
         lwd=2,
         col=rainbow(3),
         xlab = "Time (Year)",
         mark.time = T,
         ylab = "Survival rate",
         main=paste("Survival curve (p=",pValue,")",sep=""))
    legend("topright",
           paste0(gene[i],c("-NONE","-PART","-ALL")),
           lwd=2,
           col=rainbow(3)
    )
    dev.off()
  }
  }
  
  
}
target_gene<-as.data.frame(target_gene)
sgene<-merge(allsymbol,target_gene,by.x="V1",by.y="target_gene",all=FALSE)
###差异基因与生存基因韦恩图

library(VennDiagram)
survg<-sgene[,1]
diffg<-dgene[,1]
venn.diagram(
  x = list(
    'survival(32)' = survg,
    'diff(5)' = diffg
    
  ),
  filename = 'surv-diff.png',
  col = "black",
  fill = c("dodgerblue", "goldenrod1"),
  alpha = 0.5,
  cex = 0.8,
  cat.col = 'black',
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.05,
  main.cex = 1.2
)

#合并两个分析结果
agene<-merge(sgene,dgene,all=TRUE,by="V1")
agene<-agene[,-c(2:6)]
#保存所有基因
write.table(agene,file="agene.txt",sep="\t",quote=F,row.names=F)

###############富集分析及可视化
#进行一下id转化，将基因名转化为entrezid
library(clusterProfiler)
library(org.Hs.eg.db)



#所有基因富集分析
agene<-bitr(geneID = agene, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb =org.Hs.eg.db)
agene<-na.omit(agene)##删除缺失项，即symbol未能匹配到entrezid，因为一个entrezid可能对应多个symbol

ago_ALL<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="ALL",pvalueCutoff = 0.3,readable = T)
ago_BP<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="BP",pvalueCutoff = 0.3,readable = T)
ago_MF<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="MF",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = T)
ago_CC<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="CC",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = T)
ago_ALLR<-ago_ALL@result
ago_BPR<-ago_BP@result
ago_MFR<-ago_MF@result
ago_CCR<-ago_CC@result

aKEGG<-enrichKEGG(agene$ENTREZID,organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
aKEGGR<-aKEGG@result

barplot(ago_BP,drop=TRUE,showCategory = 10)
barplot(ago_MF,drop=TRUE,showCategory = 10)
barplot(ago_CC,drop=TRUE,showCategory = 10)
barplot(aKEGG,drop=TRUE,showCategory = 10)
b<-aKEGGR[(aKEGGR$p.adjust<0.05)&(aKEGGR$qvalue<0.05),]
write.table(b,file = "M_enrich.txt",sep="\t",quote=F,row.names=F)
