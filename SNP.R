rm(list=ls())
setwd("E:\\oral\\SNP")


#处理临床文件
clin=read.table("E:\\oral\\SNP\\5296cf00-4d8c-4db3-80d7-930a4b44f90d\\clinical.tsv",header = T,sep="\t")
load("E:\\oral\\DNAMethylation\\clin.Rdata")
colnames(clin)<-c("Tumor_Sample_Barcode","vital_status","follow")
#BiocManager::install("maftools")
library(maftools)

laml=read.maf(maf="E:\\oral\\SNP\\TCGA.HNSC.varscan.5296cf00-4d8c-4db3-80d7-930a4b44f90d.DR-10.0.somatic.maf",
              clinicalData =clin,isTCGA = TRUE)

a<-getGeneSummary(laml)
b<-a[a$total>10,]

#画图
plotmafSummary(laml,rmOutlier = TRUE,dashboard = TRUE,
               titvRaw = TRUE,addStat = NULL,showBarcodes = FALSE,fs=1.2,
               textSize = 1.2,color = NULL,titleSize = c(1.5,1.2),
               titvColor = NULL,top=10)
oncoplot(maf = laml,fontSize=1.2,SampleNamefontSize=1.2,titleFontSize=2)
snpdata<-laml@data


sdata<-snpdata[,c("Hugo_Symbol","Variant_Type",
                  "Variant_Classification","Tumor_Sample_Barcode")]
#删除无义突变
sdata<-sdata[which(sdata$Variant_Classification!="Nonsense_Mutation"),
             c("Hugo_Symbol","Variant_Type","Tumor_Sample_Barcode")]
#install.packages("reshape2")
library(reshape2)
rsdata=dcast(sdata,Hugo_Symbol~Tumor_Sample_Barcode)
rownames(rsdata)<-rsdata$Hugo_Symbol


#融合数据
allsymbol<-read.table(file = "E:\\oral\\finallsymbol.txt",sep="\t",header = FALSE)
sgene<-merge(rsdata,allsymbol,by.x = "Hugo_Symbol",by.y="V1",all = FALSE)
rownames(sgene)<-sgene$Hugo_Symbol
sgene<-sgene[,-1]

bgene<-merge(b,allsymbol,by.x = "Hugo_Symbol",by.y="V1",all = FALSE)

#survival analysis
survgene<-sgene
barcode <- colnames(survgene)
library(survival)
library(survminer)
##融合数据
interbarcode <- intersect(barcode,clin$Tumor_Sample_Barcode)
gene <- rownames(survgene)
dir.create("surfig")
target_gene<-c()
for(i in 1:nrow(survgene)){
  geneExp <- survgene[gene[i],interbarcode]
  geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame() 
  geneExp$Tumor_Sample_Barcode <- rownames(geneExp)
  mergdata <- merge(clin,geneExp,by.x= "Tumor_Sample_Barcode")
  if(sum(mergdata[,gene[i]])>0){
  mergdata$Group<-ifelse(mergdata[,gene[i]] >0,"mutation","normal")
  diff<-survdiff(Surv(follow, vital_status) ~ Group, data=mergdata)
    
  pValue<-1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue=signif(pValue,4)
    pValue=format(pValue,scientific=TRUE)
  }else{
    pValue=round(pValue,3)
  }
  if(pValue<0.05){
    target_gene<-c(target_gene,gene[i])
      
    fit <- survfit(Surv(follow, vital_status) ~ Group, data=mergdata)
    summary(fit)
    filename<-paste(".\\surfig\\",gene[i],"-SurvivalCurve.pdf",sep="")
    pdf(file=filename,width=5.5,height = 5)
    plot(fit,
          lwd=2,
          col=rainbow(2),
          xlab = "Time (Year)",
          mark.time = T,
          ylab = "Survival rate",
          main=paste("Survival curve (p=",pValue,")",sep=""))
    legend("topright",
            paste0(gene[i],c("-mutation","-normal")),
            lwd=2,
            col=rainbow(2)
    )
    dev.off()
  }

  }
  
}
target_gene<-as.data.frame(target_gene)
#单个基因的突变数据太少，得到的target_gene没有很大意义


agene<-as.data.frame(gene)
write.table(agene$SYMBOL,file = "agene.txt",sep="\t",quote=F,row.names=F)
write.table(bgene$Hugo_Symbol,file = "bgene.txt",sep="\t",quote=F,row.names=F)
#所有基因富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
bgene<-bitr(geneID = bgene$Hugo_Symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb =org.Hs.eg.db)
bgene<-na.omit(bgene)##删除缺失项，即symbol未能匹配到entrezid，因为一个entrezid可能对应多个symbol


ago_ALL<-enrichGO(bgene$ENTREZID,'org.Hs.eg.db',ont="ALL",pvalueCutoff = 0.3,readable = T)
ago_BP<-enrichGO(bgene$ENTREZID,'org.Hs.eg.db',ont="BP",pvalueCutoff = 0.3,readable = T)
ago_MF<-enrichGO(bgene$ENTREZID,'org.Hs.eg.db',ont="MF",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = T)
ago_CC<-enrichGO(bgene$ENTREZID,'org.Hs.eg.db',ont="CC",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = T)
ago_ALLR<-ago_ALL@result
ago_BPR<-ago_BP@result
ago_MFR<-ago_MF@result
ago_CCR<-ago_CC@result

aKEGG<-enrichKEGG(bgene$ENTREZID,organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
aKEGGR<-aKEGG@result

barplot(ago_BP,drop=TRUE,showCategory = 10)
barplot(ago_MF,drop=TRUE,showCategory = 10)
barplot(ago_CC,drop=TRUE,showCategory = 10)
barplot(aKEGG,drop=TRUE,showCategory = 10)

d<-rbind(ago_BPR[(ago_BPR$p.adjust<0.05)&(ago_BPR$qvalue<0.05),],ago_CCR[(ago_CCR$p.adjust<0.05)&(ago_CCR$qvalue<0.05),])
d<-rbind(d,ago_MFR[(ago_MFR$p.adjust<0.05)&(ago_MFR$qvalue<0.05),])
d<-rbind(d,aKEGGR[(aKEGGR$p.adjust<0.05)&(aKEGGR$qvalue<0.05),])
write.table(d,file = "SNP_enrich.txt",sep="\t",quote=F,row.names=F)
