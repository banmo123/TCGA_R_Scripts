#设置工作路径
setwd("E:\\oral")


##########################1.获取节律基因SYMBOL###########
allsymbol<-read.table(file = "allsymbol.txt",sep = "\t",header=TRUE)
allsymbol<-allsymbol[!duplicated(allsymbol$Gene.name),]
write(allsymbol,file = "finallsymbol.txt",sep = "\t")
allsymbol<-as.data.frame(allsymbol)
colnames(allsymbol)<-"gene_name"
rownames(allsymbol)<-allsymbol[,"gene_name"]

##########################2.处理临床信息################
#读入临床文件名
xmlFileNames<-dir(path="./clin",full.names = T,pattern = "xml$",recursive = T)
library(XML)
AllPaticliniList <- lapply(xmlFileNames,
                           function(x){#x=xmlFileNAMES[1]
                             result <- xmlParse(file=x)
                             rootnode <- xmlRoot(result)
                             xmldataframe <- xmlToDataFrame(rootnode[2])
                             return(t(xmldataframe))})

for(i in 1:365){
  AllPaticliniList[[i]] <- data.frame(t(data.frame(AllPaticliniList[[i]])))
  
}
#设置一个for循环，使Allpaticlinilist中的列表全部转化为数据框形式，因为rbind.fill函数
#只能作用于数据框，365为文件数量
library(plyr)
#整理为数据框的形式的临床信息
AllPaticliniList <- do.call(rbind.fill,AllPaticliniList)

#提取需要的信息
pclinData <-data.frame(Barcode=AllPaticliniList[,"bcr_patient_barcode"],
                       vital_status=AllPaticliniList[,"vital_status"],
                       days_to_death=AllPaticliniList[,"days_to_death"],
                       day_to_last_known_alive=AllPaticliniList[,"days_to_last_known_alive"],
                       last_follow_up_time=AllPaticliniList[,"days_to_last_followup"])

survivaltime <- c()

for(id in row.names(pclinData)) {
  days_to_death=as.numeric(pclinData[id,"days_to_death"])
  last_follow_up_time=as.numeric(pclinData[id,"last_follow_up_time"])
  patients_survival_time <-c()
  if(!is.na(days_to_death)){
    patients_survival_time <- c(days_to_death)
  }else{patients_survival_time <- c(last_follow_up_time)}
  survivaltime <- c(survivaltime,patients_survival_time)
}

getclinicaldata <-function(pclinData){
  tempdata<-pclinData
  finallypclindata<-data.frame(Barcode=tempdata[,"Barcode"],
                               vital_status=tempdata[,"vital_status"],
                               survivaltime=survivaltime
  )
  return(finallypclindata)
}
tidypclindata<-getclinicaldata(pclinData)


#去除数据框中生存状态缺失和为na的项dplyr中的filter函数
library(dplyr)
tidypclindata=tidypclindata%>%filter(is.na(vital_status)==FALSE)
tidypclindata=tidypclindata%>%filter(vital_status!='')                                     

#将筛选过后的生存状态项转变为01，便于生存分析
tidypclindata$vital_status<-ifelse(tidypclindata$vital_status=="Alive",0,1)


#整理一下生存时间和列名
clin<-tidypclindata
clin<-filter(clin,survivaltime!="--")
names(clin)<-c("sample","survivalstate","follow")
clin<-clin[,c("sample","survivalstate","follow")]
clin$follow<-as.numeric(clin$follow)/365

save(clin,file = "clin.Rdata")
##########################3.处理转录组数据################

###处理json文件
library(rjson)
jsonFile<-fromJSON(file="./RNAseqC/metadata.cart.2021-05-13.json")
filesNameToBarcode<-data.frame(filesName=c(),TCGA_Barcode=c())
for(i in 1:length(jsonFile)){
  TCGA_Barcode<-jsonFile[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name<-jsonFile[[i]][["file_name"]]
  filesNameToBarcode<-rbind(filesNameToBarcode,data.frame(fileName=file_name,TCGA_Barcode=TCGA_Barcode))
}
row.names(filesNameToBarcode)<-filesNameToBarcode[,1]

###读入转录组数据
filepath <- dir(path ="./RNAseqC",pattern=".gz$",full.names = TRUE,recursive = T)
exp <- data.frame()
for(wd in filepath){
  read.table(wd,header = F,sep = "\t")
  #每一个循环读取一个文件
  oneSampExp <- read.table(wd,header =FALSE)
  tempPath <- unlist(strsplit(wd,"/"))
  filename <- tempPath[length(tempPath)]
  print(paste0("input your data now:",filename,""))
  #根据filenametobarcode文件中的文件名称与barcode对应关系，命名列名
  colnames(oneSampExp) <- c("Ensembl",filesNameToBarcode[filename,"TCGA_Barcode"])
  if (dim(exp)[1]== 0){
    exp <- oneSampExp
  }
  else 
  {exp <- merge(exp,oneSampExp,by = "Ensembl")}
}
exp<-exp[-c(1:5),]
exp$Ensembl <- substr(exp[,"Ensembl"],1,15)

# 一个函数，通过gtf文件获取Ensemble_ID与基因名称的对应关系
get_map = function(input) {
  if (is.character(input)) {
    if(!file.exists(input)) stop("Bad input file.")
    message("Treat input as file")
    input = data.table::fread(input, header = FALSE)
  } else {
    data.table::setDT(input)
  }
  input = input[input[[3]] == "gene", ]
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  pattern_type = ".*gene_type \"([^;]+)\";.*"
  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  gene_type = sub(pattern_type, "\\1", input[[9]])
  EnsemblTOGenename <- data.frame(Ensembl = gene_id,
                                  gene_name = gene_name,
                                  gene_type = gene_type,
                                  stringsAsFactors = FALSE)
  return(EnsemblTOGenename)
}
EnsemblTOGenename <- get_map(".\\gencode.v37.annotation.gtf") 
#去除版本号ENSG00000223972(.1)
EnsemblTOGenename$Ensembl <- substr(EnsemblTOGenename[,"Ensembl"],1,15)

#融合数据
transdata<-merge(EnsemblTOGenename,exp,by="Ensembl",all = FALSE)
relatedgene<-merge(transdata,allsymbol,by="gene_name",all = FALSE)
relatedgene<-relatedgene[,-c(2:3)]
relatedgene<-relatedgene[!duplicated(relatedgene$gene_name),]
rownames(relatedgene)<-relatedgene[,"gene_name"]


#barcode <- gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(tumorexp))
#colnames(tumorexp) <- barcode 



#肿瘤表达矩阵和正常组织表达矩阵
tumor <- colnames(relatedgene)[as.integer(substr(colnames(relatedgene),14,15)) < 10]
normal <- colnames(relatedgene)[as.integer(substr(colnames(relatedgene),14,15)) >= 10] 
tumor_sample <- relatedgene[,tumor[2:331]]
normal_sample <- relatedgene[,normal[2:33]] 
expSet_by_group <- cbind(normal_sample,tumor_sample) 
##################DESeq2包差异分析
BiocManager::install("DESeq2")
library(DESeq2)
group_list <- c(rep("normal",ncol(normal_sample)),rep("tumor",ncol(tumor_sample)))
condition=factor(group_list)
coldata<-data.frame(row.names = colnames(expSet_by_group),condition)
dds<-DESeqDataSetFromMatrix(countData = expSet_by_group,colData=coldata,design =~condition)
dds$condition<-relevel(dds$condition,ref="normal")#指定normal为对照组
#差异表达矩阵
dds<-DESeq(dds)
allDEG2<-as.data.frame(results(dds))
#添加change列
logFC_cutoff=1
allDEG2$change=as.factor(
  ifelse(allDEG2$pvalue<0.05&abs(allDEG2$log2FoldChange)>logFC_cutoff,
         ifelse(allDEG2$log2FoldChange>logFC_cutoff,"UP","DOWN"),'NOT')#差异倍数大于1，p值小于0.05，肿瘤对基因存在影响
)
#提取基因显著差异的差异矩阵
padj=0.05
foldChange=2
nrDEG_DESeq2_signif=allDEG2[(allDEG2$padj<padj&
                              abs(allDEG2$log2FoldChange)>foldChange),]
nrDEG_DESeq2_signif=nrDEG_DESeq2_signif[order(nrDEG_DESeq2_signif$log2FoldChange),]
nrDEG_DESeq2_signif<-na.omit(nrDEG_DESeq2_signif)




##################edgeR包做差异分析
library(edgeR)
library(limma)
group_list <- c(rep("normal",ncol(normal_sample)),rep("tumor",ncol(tumor_sample)))
group_list=factor(group_list)
design<-model.matrix(~0+group_list)
rownames(design)=colnames(expSet_by_group)
colnames(design)<-levels(group_list)
#差异表达矩阵
DEGlist<-DGEList(counts=expSet_by_group,group=group_list)
keep_gene<-rowSums(cpm(DEGlist)>1)>=2
table(keep_gene)
DEGlist<-DEGlist[keep_gene,,keep.lib.sizes=FALSE]

DEGlist<-calcNormFactors(DEGlist)
DEGlist<-estimateGLMCommonDisp(DEGlist,design )
DEGlist<-estimateGLMTrendedDisp(DEGlist,design)
DEGlist<-estimateGLMTagwiseDisp(DEGlist,design)

fit <- glmFit(DEGlist, design)
results <- glmLRT(fit, contrast = c(-1, 1)) 
nrDEG_edgeR <- topTags(results, n = nrow(DEGlist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
#添加change列
#logFC_cutoff<-with(nrDEG_edgeR,mean(abs(logFC))+2*sd(abs(logFC)))
logFC_cutoff=1
nrDEG_edgeR$change=as.factor(
  ifelse(nrDEG_edgeR$PValue<0.05&abs(nrDEG_edgeR$logFC)>logFC_cutoff,
         ifelse(nrDEG_edgeR$logFC>logFC_cutoff,"UP","DOWN"),'NOT')#差异倍数大于1，p值小于0.05，肿瘤对基因存在影响
)
head(nrDEG_edgeR)

padj = 0.05 # 自定义
foldChange= 2 # 自定义
nrDEG_edgeR_signif  = nrDEG_edgeR[(nrDEG_edgeR$FDR < padj & 
                                     abs(nrDEG_edgeR$logFC)>foldChange),]
                                     
nrDEG_edgeR_signif = nrDEG_edgeR_signif[order(nrDEG_edgeR_signif$logFC),]


 
##################limma包做差异分析
group_list <- c(rep("normal",ncol(normal_sample)),rep("tumor",ncol(tumor_sample)))
group_list=factor(group_list)
design<-model.matrix(~0+group_list)
rownames(design)=colnames(expSet_by_group)
colnames(design)<-levels(group_list)

DGElist <- DGEList( counts = expSet_by_group, group = group_list )
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 # 自定义
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
##limma-trand
logCPM <- cpm(DGElist, log=TRUE, prior.count=3)
design<-model.matrix(~group_list)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

nrDEG_limma_trand = topTable(fit, coef=ncol(design))
nrDEG_limma_trand = na.omit(nrDEG_limma_trand)
head(nrDEG_limma_trand)
#添加change列
logFC_cutoff=1
nrDEG_limma_trand$change=as.factor(
  ifelse(nrDEG_limma_trand$P.Value<0.05&abs(nrDEG_limma_trand$logFC)>logFC_cutoff,
         ifelse(nrDEG_limma_trand$logFC>logFC_cutoff,"UP","DOWN"),'NOT')#差异倍数大于1，p值小于0.05，肿瘤对基因存在影响
)
##limma-voom
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
fit1 <- lmFit(v, design)
fit2 <- eBayes(fit1)
nrDEG_limma_voom = topTable(fit2, coef = ncol(design))
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
#添加change列
logFC_cutoff=1
nrDEG_limma_voom$change=as.factor(
  ifelse(nrDEG_limma_voom$P.Value<0.05&abs(nrDEG_limma_voom$logFC)>logFC_cutoff,
         ifelse(nrDEG_limma_voom$logFC>logFC_cutoff,"UP","DOWN"),'NOT')#差异倍数大于1，p值小于0.05，肿瘤对基因存在影响
)
padj = 0.05 # 自定义
foldChange= 2 # 自定义
nrDEG_limma_trand_signif = nrDEG_limma_trand[(nrDEG_limma_trand$adj.P.Val < padj & 
                                              abs(nrDEG_limma_trand$logFC)>foldChange ),]
nrDEG_limma_trand_signif = nrDEG_limma_trand_signif[order(nrDEG_limma_trand_signif$logFC),]

nrDEG_limma_voom_signif = nrDEG_limma_voom[(nrDEG_limma_voom$adj.P.Val < padj & 
                                                abs(nrDEG_limma_voom$logFC)>foldChange),]
nrDEG_limma_voom_signif = nrDEG_limma_voom_signif[order(nrDEG_limma_voom_signif$logFC),]

#所有方法得到的差异基因
diffgene<-c()
diffgene<-c(rownames(nrDEG_limma_trand_signif),rownames(nrDEG_limma_voom_signif),
            rownames(nrDEG_edgeR_signif),rownames(nrDEG_DESeq2_signif))
diffgene<-diffgene[!duplicated(diffgene),]
diffgene<-as.data.frame(diffgene)
#############画图


####比较三个包差异分析的结果——韦恩图
edgeR = rownames(nrDEG_edgeR_signif)
dim(nrDEG_edgeR_signif)
limma_trand = rownames(nrDEG_limma_trand_signif)
dim(nrDEG_limma_trand_signif)
limma_voom = rownames(nrDEG_limma_voom_signif)
dim(nrDEG_limma_voom_signif)
DESeq2 = rownames(nrDEG_DESeq2_signif)
dim(nrDEG_DESeq2_signif)

library(VennDiagram)
venn.diagram(
  x = list(
    'edgeR(85)' = edgeR,
    'limma_trand(9)' = limma_trand,
    'DESeq2(88)' = DESeq2,
    'limma_voom(8)'=limma_voom
  ),
  filename = 'VN.png',
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1","black"),
  alpha = 0.5,
  cex = 0.8,
  cat.col = 'black',
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.05,
  main.cex = 1.2
)

pca_plot=function(exp,group){
   dat=as.data.frame(t(exp))
   #install.packages("FactoMineR")
   #install.packages("factoextra")
   library(FactoMineR)
   library(factoextra)
   dat.pca<-PCA(dat,graph=FALSE)
   fviz_pca_ind(dat.pca,
                geom.ind="point",
                col.ind=group_list,
                addEllipses=TRUE,
                legend.title="Groups")
 }

#pca分析
dat=log(expSet_by_group+1)
dat<-na.omit(dat)
pca=pca_plot(dat,group_list)


#差异基因与转录组基因融合
#surgene<-merge(etSig,relatedgene,by="gene_name")
#rownames(surgene)<-surgene[,"gene_name"]
#surgene<-surgene[,-c(1:5)]
#surgene<-surgene[,tumor[2:331]]

survgene<-relatedgene
###c处理列名
barcode <- gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(survgene))
colnames(survgene) <- barcode 
######################### 生存分析
library(survival)
library(survminer)
##融合数据
interbarcode <- intersect(barcode,clin$sample)
gene <- rownames(survgene)
dir.create("surfig")
target_gene<-c()
for(i in 1:nrow(survgene)){
  geneExp <- survgene[gene[i],interbarcode]
  geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame() 
  geneExp$sample <- rownames(geneExp)
  mergdata <- merge(clin,geneExp,by = "sample")
  
  mergdata$Group <- ifelse(mergdata[,gene[i]] > median(mergdata[,gene[i]]),"High","Low")
  
  
  
  
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
       col=rainbow(2),
       xlab = "Time (Year)",
       mark.time = T,
       ylab = "Survival rate",
       main=paste("Survival curve (p=",pValue,")",sep=""))
  legend("topright",
         paste0(gene[i],c("-High","-Low")),
         lwd=2,
         col=rainbow(2)
  )
  dev.off()
}
}
target_gene<-as.data.frame(target_gene)
save(target_gene,file = "target_gene.Rdata")
write.table(target_gene,file="target_gene.txt",sep="\t",quote=F,row.names=F)

#看一下差异分析和生存分析基因差异
venn.diagram(
  x = list(
    'diff(93)' = diffgene[,1],
    'survival(39)' = target_gene[,1]
  ),
  filename = 'diff-surv.png',
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
#所有基因
a<-merge(target_gene,diffgene,by.x="target_gene",by.y="diffgene",all=TRUE)
write.table(a,file="a_gene.txt",sep="\t",quote=F,row.names=F)

###############富集分析及可视化
#进行一下id转化，将基因名转化为entrezid
library(org.Hs.eg.db)
library(dplyr)


#所有基因富集分析
agene<-bitr(geneID = a[,1], 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb =org.Hs.eg.db)
agene<-na.omit(agene)##删除缺失项，即symbol未能匹配到entrezid，因为一个entrezid可能对应多个symbol

ago_ALL<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="ALL",pvalueCutoff = 0.3,readable = T)
ago_BP<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="BP",pvalueCutoff = 0.3,readable = T)
ago_MF<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="MF",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = T)
ago_CC<-enrichGO(agene$ENTREZID,'org.Hs.eg.db',ont="CC",pvalueCutoff = 0.3,qvalueCutoff = 0.5,readable = T)
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

