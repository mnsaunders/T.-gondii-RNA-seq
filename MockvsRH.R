library('DESeq2')
library('tidyverse')
library("ggrepel")
library('extrafont')
font_import()
loadfonts(device = "win")

data<-MockvsRH %>% remove_rownames %>% column_to_rownames(var="Geneid")
Counts<-filter(data, "Mock_1" >0, "Mock_2" >0,"Mock_3" >0, "Mock_4" >0)

condition <- factor(c("RH","RH","RH","RH","Mock","Mock","Mock","Mock"))
coldata<-data.frame(row.names=colnames(Counts), condition)

dds<-DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

deseq<-DESeq(dds)

vsdata<-vst(deseq, blind=FALSE)
pcaO<-plotPCA(vsdata, intgroup="condition")
pcaO + geom_text_repel(aes(label = name), size=2.5)+
  geom_point(size=3)+
  theme_bw()

plotDispEsts(deseq)
res<- results(deseq, contrast=c("condition", "Mock", "RH"))
res
nona<-na.omit(res)
nona
nona$diffexp<-'insignificant'
nona$diffexp[nona$log2FoldChange>0 & nona$padj<0.05] <-'up'
nona$diffexp[nona$log2FoldChange<0 & nona$padj<0.05] <-'down'
write.csv(nona, "RHinfections.csv")
RHinfections1<-filter(RHinfections, padj<0.05)
names(RHinfections1)[names(RHinfections1) == '...1'] <- 'gene_id'

comdataup<-subset(res, padj<1e-40)
comdataup<-subset(comdataup, log2FoldChange>1.5)
comdataup<-as.data.frame(comdataup)
comdataup<-rownames_to_column(comdataup)

comdatadown<-subset(res, padj<1e-10)
comdatadown<-subset(comdatadown, log2FoldChange<1.5)
comdatadown<-as.data.frame(comdatadown)
comdatadown<-rownames_to_column(comdatadown)

normalized_counts<-counts(deseq, normalized=T)
Transfromed_counts<-log10(normalized_counts+1)
df<-as.data.frame(Transfromed_counts)
df<-rownames_to_column(df)
up<-df%>%left_join(comdataup, by='rowname')
up<-na.omit(up)
up<-up[,-10:-15]
up<-remove_rownames(up)
up<-column_to_rownames(up)

pheatmap(up, fontsize = 8)

down<-df%>%left_join(comdatadown, by='rowname')
down<-na.omit(down)
down<-down[,-10:-15]
down<-remove_rownames(down)
down<-column_to_rownames(down)

pheatmap(down, fontsize = 8)





genedat<-RHinfections1 %>%
  left_join(Human_annotation_table, by = "gene_id")

ggplot(genedat, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexp))+
  geom_point()+ geom_text_repel(aes(label = gene), check_overlap = TRUE, size=3, fontface = "bold", show.legend = FALSE)+
  geom_vline(xintercept = c(0,0), col="grey")+ theme_bw()+
  scale_color_manual(values=c("red3","grey","navy"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size=13, family="roboto"), legend.title = element_blank(), legend.position="bottom")+
  labs(x="log2FoldChange", y="-log(pvalue)", title= "T. gondii Affects Gene Regulation in Neutrophils")

up<-filter(genedat, log2FoldChange> 6.5)
down<-filter(genedat, log2FoldChange< -3)
top<-rbind(up,down)
top<-filter(top, gene!="NA")


ggplot(top, aes(y=log2FoldChange, x=gene, fill=diffexp))+
  geom_bar(stat="identity")+ geom_hline(yintercept = c(0,0), col="grey")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),axis.text.x=element_text(angle=80,hjust=1,vjust=1, size = 10, face="bold.italic"),legend.title = element_blank(), text = element_text(size=13, family="roboto"), legend.position = "none")+
  scale_fill_manual(values=c("red3","navy"))+
  labs(y="Gene Regulation", x=" ", title= "Toxoplasma gondii Affects Lipopolysaccharide Infected Neutrophils")

