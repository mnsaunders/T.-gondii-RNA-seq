library('DESeq2')
library('tidyverse')
library("ggrepel")

data<-neutrocompare %>% remove_rownames %>% column_to_rownames(var="Geneid")
Counts<-filter(data, "RH_1" >0, "RH_2" >0,"RH_3" >0, "RH_4" >0)

condition <- factor(c("RHLPS","RHLPS","RHLPS","RHLPS","RH","RH","RH","RH"))
coldata<-data.frame(row.names=colnames(Counts), condition)

dds<-DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

deseq<-DESeq(dds)

vsdata<-vst(deseq, blind=FALSE)
pcaO<-plotPCA(vsdata, intgroup="condition")
pcaO + geom_text_repel(aes(label = name), size=2.5)+
  geom_point(size=3)+
  theme_bw()

plotDispEsts(deseq)
res<- results(deseq, contrast=c("condition", "RHLPS", "RH"))
res
nona<-na.omit(res)
nona
nona$diffexp<-'insignificant'
nona$diffexp[nona$log2FoldChange>0 & nona$padj<0.05] <-'up'
nona$diffexp[nona$log2FoldChange<0 & nona$padj<0.05] <-'down'
write.csv(nona, "neutrocom.csv")
neutrocom1<-filter(neutrocom, padj<0.05)
names(neutrocom1)[names(neutrocom1) == '...1'] <- 'gene_id'

genedat<-neutrocom1 %>%
  left_join(Human_annotation_table, by = "gene_id")

ggplot(genedat, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexp))+
  geom_point()+ geom_text_repel(aes(label = gene), check_overlap = TRUE, size=3, fontface = "bold")+
  geom_vline(xintercept = c(0,0), col="grey")+ theme_bw()+
  scale_color_manual(values=c("red3","grey", "green3"))+
  labs(x="Gene Regulation", title="RHLPS vs RH neutrophils")

up<-filter(genedat, log2FoldChange> 7)
down<-filter(genedat, log2FoldChange< -3)
top<-rbind(up,down)
top<-filter(top, gene!="NA")


ggplot(top, aes(y=log2FoldChange, x=gene, fill=diffexp))+
  geom_bar(stat="identity")+ geom_hline(yintercept = c(0,0), col="grey")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),axis.text.x=element_text(angle=80,hjust=1,vjust=1, size = 5.5, face="bold.italic"))+
  labs(y="Gene Regulation", title="RHLPS vs RH neutrophils")

