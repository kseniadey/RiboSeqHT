library(statmod)
library(magrittr)
library(dplyr)
library(edgeR)
library(ggplot2)

# Chapter 4 https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

df = read.csv('/Users/kseniadey/Desktop/NGS/HSE_RiboSeq_HT/01. RiboSeq_RNASeq_HCC_counts.tsv', sep = '\t')
df = replace(df, is.na(df), 0)
# Извлечем все RNASeq результаты
RNA = df%>% select(ends_with("RNA"))
RNA = cbind(df$geneSymbol, RNA)
# Извлечем все RiboSeq результаты
RPF = df%>% select(ends_with("RPF"))
RPF = cbind(df$geneSymbol, RPF)
# Сразу подготовим факторы (номер и состояние ткани
patient = factor(c('LC001', 'LC001', 'LC033', 'LC033', 'LC034', 'LC034', 'LC501', 'LC501', 'LC502', 'LC502', 'LC505', 'LC505', 'LC506', 'LC506', 'LC507', 'LC507', 'LC508','LC508', 'LC509', 'LC509'))
tissue= factor(c('Norm', 'Tumor', 'Norm', 'Tumor', 'Norm', 'Tumor', 'Norm', 'Tumor', 'Norm', 'Tumor','Norm', 'Tumor', 'Norm', 'Tumor', 'Norm', 'Tumor', 'Norm', 'Tumor', 'Norm', 'Tumor'))

# Дифф.экспрессия для данных RNASeq
RNA.DE = DGEList(counts = RNA[,2:21], genes = RNA[,1])
data.frame(Sample=colnames(RNA.DE), patient, tissue)
design.RNA = model.matrix(~ patient + tissue)
rownames(design.RNA) = colnames(RNA.DE)
design.RNA
 
RNA.DE = estimateDisp(RNA.DE, design.RNA, robust=TRUE) 
RNA.DE$common.dispersion # квадратный корень common.dispersion = коэффициент биологической вариации. То коэффициент биологической вариации приблизитоельно 0,614
plotBCV(RNA.DE)
fit.RNA = glmFit(RNA.DE, design.RNA)
lrt.RPF = glmLRT(fit.RNA) # likelihood ratio для разницы в генах между нормальной и раковой тканями
topTags(lrt.RPF) 
colnames(design.RNA)
# Посчитаем counts per million
o = order(lrt.RPF$table$PValue)
cpm(RNA.DE)[o[1:10],] # разница ме-ду здоровыми и раковыми образцами явно заметна
summary(decideTests(lrt.RPF))
plotMD(lrt.RPF)
abline(h=c(-1, 1), col="blue")

# Дифф.экспрессия для данных RiboSeq 
# Аналогичная процедура
RPF.DE = DGEList(counts = RPF[,2:21], genes = RPF[,1])
data.frame(Sample=colnames(RPF.DE), patient, tissue)
design.RPF = model.matrix(~ patient + tissue)
rownames(design.RPF ) = colnames(RPF.DE)
design.RPF 
RPF.DE = estimateDisp(RPF.DE, design.RPF , robust=TRUE)
RPF.DE$common.dispersion # 0.4871409 -> коэффициент биологической вариации составляет 0,69795
plotBCV(RPF.DE)
fit.RPF = glmFit(RPF.DE, design.RPF)
lrt.RPF = glmLRT(fit.RPF)
topTags(lrt.RPF)
colnames(design.RPF )
o = order(lrt.RPF$table$PValue)
cpm(RPF.DE)[o[1:10],]
summary(decideTests(lrt.RPF))
plotMD(lrt.RPF)
abline(h=c(-1, 1), col="blue")

# Volcano plots
ggplot(data = df, aes(x = lrt$table$logFC, y = -log10(lrt$table$PValue))) + geom_point()
ggplot(data = df, aes(x = lrt1$table$logFC, y = -log10(lrt1$table$PValue))) + geom_point()
# Профили дифференциальной экспрессии генов отличаются



