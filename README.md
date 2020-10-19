# main analysis scripts
## :open_file_folder:Input

- Gene annotation:  dmel-all-r6.25.gtf.gz, download it form flybase <ftp://ftp.flybase.net/releases/FB2018_06/dmel_r6.25/gtf/dmel-all-r6.25.gtf.gz>
- Selected genes: `mitoch_geneID.txt`
- The raw microarray data was deposited under the EBI ArrayExpress  [E-MEXP-1817](https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-1817/?query=Mortin).
    - `C1_8-1-08_s1.CEL`
    - `C2_8-1-08_s1.CEL`
    - `C3_8-1-08_s1.CEL`
    - `E1_8-1-08_s1.CEL`
    - `E2_8-1-08_s1.CEL`
    - `E3_8-1-08_s1.CEL`

## :file_folder: output

- Table
    - `DEG_signal.csv`
    - `gsea_C5_output.txt`
    - `gsym.GO.term.down.csv`
    - `gsym.GO.term.up.csv`
    - `limma_output.csv`
    - `signal_gsym.csv`
    - `signal_probeset.csv`
- Figure
    - `GSEA_bubbles.pdf`
    - `heatmap_mitoch_scale.pdf`
    - `valcono.pdf`

## :key:Install packages required

```{r eval=FALSE}
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

devtools::install_git('https://github.com/GuangchuangYu/clusterProfiler')
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("affy")
BiocManager::install("limma")
#Affymetrix Drosophila Genome 2.0 Array annotation data (chip drosophila2)
#https://bioconductor.org/packages/3.8/data/annotation/
BiocManager::install("drosophila2.db")
install.packages("msigdbr")
install.packages("DT")
```

```{r}
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
```

## :one: QC

```{r, fig.width=6, fig.height=6}
library(affy)
(AffyData <- ReadAffy())

#Image of the log intensities.
image(AffyData)

# Histogram of PM intensities
hist(AffyData)

# RNA degradation plots
deg <- AffyRNAdeg(AffyData)
plotAffyRNAdeg(deg) 

#Pairwise MA plots
MAplot(AffyData, pairs=TRUE, plot.method="smoothScatter")
```

## :two: Convert raw probe intensities to expression values

```{r}
normmethod <- "rma"
if (normmethod == "lw"){
    eset <- expresso(AffyData, normalize.method="qspline",
               bgcorrect.method="rma",pmcorrect.method="pmonly",
               summary.method="liwong")
    
    exp.eset <- exprs(eset) 
    exp.eset <- log2(exp.eset)
} else if (normmethod == "rma") {
  eset <- rma(AffyData)
  exp.eset <- exprs(eset) 
}

par(mfrow = c(1,2), las = 2, 
    xpd = T, mar = par()$mar + c(5,0,0,0))
boxplot(AffyData)
title("before")
boxplot(exp.eset)
title("after")
par(xpd=FALSE)

write.csv(exp.eset, file="Output/Table/signal_probeset.csv", quote = F)
```

Add gene symbol to expression values

```{r}
eset <- read.csv("Output/Table/signal_probeset.csv", row.names = 1)
dim(eset)
library(AnnotationDbi)
library(drosophila2.db)
probe.gsym <- select(drosophila2.db, keys=keys(drosophila2.db), 
       columns = "SYMBOL")
head(probe.gsym)
idvec <- probe.gsym$SYMBOL
names(idvec) <- probe.gsym$PROBEID
eset$SYMBOL <- idvec[rownames(eset)]

# remove all probes that do not have an Entrez Gene ID and Symbol
HasSymbol <- !is.na(eset$SYMBOL)
eset <- eset[HasSymbol,]
dim(eset)

# combine probe sets with same gene symbol
eset_uniq <- aggregate(.~SYMBOL,eset,median)
dim(eset_uniq)

# set gene symbol as rowname
rownames(eset_uniq) <- eset_uniq$SYMBOL
eset_uniq <- subset(eset_uniq, select = -SYMBOL)
write.csv(eset_uniq, file="Output/Table/signal_gsym.csv",row.names = T,quote = F)
```

## :three: Statistics for differential expression

```{r}
library(limma)
eset_uniq <- read.csv("Output/Table/signal_gsym.csv", row.names = 1)

mulg.design <- model.matrix(~ 0 + factor(c(1,1,1,2,2,2))) 
colnames(mulg.design) <- c("C", "E")
mulg.fit <- lmFit(eset_uniq, mulg.design)
mulg.matrix <- makeContrasts(E-C, levels = mulg.design)
mulg.fit <- contrasts.fit(mulg.fit, mulg.matrix)
mulg.fit <- eBayes(mulg.fit)

results <- decideTests(mulg.fit)
summary(results)

write.csv(topTable(mulg.fit, coef=1, number = nrow(mulg.fit), adjust="BH"), file="Output/Table/limma_output.csv")
```

Differential expression genes

```{r}
# check the number of differentially expressed genes
allgeneslimma <- read.csv("Output/Table/limma_output.csv", row.names = 1)

allgeneslimma$threshold <- as.factor(
  ifelse(allgeneslimma$adj.P.Val <= 0.05 & abs(allgeneslimma$logFC) >= 1, "DEG", "Ambiguous"))

DEG <- allgeneslimma[allgeneslimma$threshold == "DEG",]
dim(DEG)
head(DEG)

DEG_sortbyFC <- DEG[order(abs(DEG$logFC), decreasing = T),]

eset_DEG <- eset_uniq[rownames(DEG_sortbyFC),] #sorted by abs(FC) and p <= 0.05
head(eset_DEG)
write.csv(eset_DEG, "Output/Table/DEG_signal.csv", quote = F, row.names = T)
```

How many genes in Drosophila

```{r}
dmgtf <- read.table("Input/dmel-all-r6.25.gtf.gz", sep = "\t")
dmgene <- dmgtf[dmgtf$V3 == "gene",]
dim(dmgene)
12874/nrow(dmgene)
```

Add limma output to GO and GO term

```{r}
library(msigdbr)

C5 <- msigdbr(species = "Drosophila melanogaster", 
             category = "C5") 
head(C5)
library(dplyr)
gsym.GO.term <- select(C5, gene_symbol, gs_name, gs_subcat, human_gene_symbol)

allgeneslimma <- read.csv("Output/Table/limma_output.csv", row.names = 1)
allgeneslimma$gene_symbol <- row.names(allgeneslimma)

# add logFC
idvec <- allgeneslimma$logFC
names(idvec) <- row.names(allgeneslimma)
gsym.GO.term$logFC <- idvec[gsym.GO.term$gene_symbol]

# add AveExpr
idvec <- allgeneslimma$AveExpr
names(idvec) <- row.names(allgeneslimma)
gsym.GO.term$AveExpr <- idvec[gsym.GO.term$gene_symbol]

# add Pvalue
idvec <- allgeneslimma$adj.P.Val
names(idvec) <- row.names(allgeneslimma)
gsym.GO.term$adj.P.Val <- idvec[gsym.GO.term$gene_symbol]

# add threshold
gsym.GO.term$threshold <- as.factor(
  ifelse(gsym.GO.term$adj.P.Val <= 0.05 & abs(gsym.GO.term$logFC) >= 1, 
         ifelse(gsym.GO.term$logFC >= 1, "Increase", "Decrease"), "Ambiguous"))

summary(gsym.GO.term$threshold)

# filter increase and decrease
gsym.GO.term.up <- gsym.GO.term[grep(pattern="Increase", gsym.GO.term$threshold),]
gsym.GO.term.down <- gsym.GO.term[grep(pattern="Decrease", gsym.GO.term$threshold),]

# Table S1
write.csv(gsym.GO.term.up, "Output/Table/gsym.GO.term.up.csv", quote = F, row.names = F)
write.csv(gsym.GO.term.down, "Output/Table/gsym.GO.term.down.csv", quote = F, row.names = F)
```

## :four: GSEA

```{r}
allgeneslimma <- read.csv("Output/Table/limma_output.csv", row.names = 1)
dim(allgeneslimma)
gsym.fc <- data.frame(gsym = row.names(allgeneslimma), logFC = allgeneslimma$logFC)

gsym.fc <- gsym.fc[order(gsym.fc$logFC, decreasing = T),]
head(gsym.fc)
gene.expr <- gsym.fc$logFC
names(gene.expr) <- gsym.fc$gsym
```

GSEA with C5

```{r}
library(clusterProfiler)
library(msigdbr)
category <- "C5"
gmt <- msigdbr(species = "Drosophila melanogaster", 
             category = category)
gmt2 <- gmt%>%
  dplyr::select(gs_name, gene_symbol, gs_cat)

y <- GSEA(gene.expr, pvalueCutoff = 1,
          TERM2GENE = gmt2[, 1:2],
          TERM2NAME = gmt2[, c(1, 3)])

dim(y@result[y@result$pvalue <= 0.01, ])

# Table S2
write.table(y@result[y@result$pvalue <= 0.01, ], "Output/Table/gsea_C5_output.txt" , sep = "\t", quote = F, row.names = F)
```

Draw prety GSEA bubbles

```{r, fig.width=12, fig.height=20}
library(stringr)
library(ggplot2)

x <- read.table("Output/Table/gsea_C5_output.txt", sep = "\t", header = T, row.names = 1)
x <- x[order(x$enrichmentScore, decreasing = T),]
x <- x[c(1:15,(nrow(x)-14):nrow(x)),]

x$group <- cut(x$enrichmentScore, breaks = c(-Inf, 0, Inf),labels = c("Decrease","Increase"))
x$pathway <- str_replace_all(x$Description, "GO_", "")
x$pathway <- tolower(str_replace_all(x$pathway, "_", " "))

sortx <- x[order(x$enrichmentScore, decreasing = F),]
sortx$pathway <- factor(sortx$pathway, levels = sortx$pathway)

ggplot(sortx, aes(enrichmentScore, pathway, colour = enrichmentScore ))+
  geom_point(aes(size=(-log10(pvalue))))+ 
  scale_color_gradientn(colours=c("blue","red"))+ 
  scale_size_continuous(range = c(2,10))+ 
  #scale_x_continuous(limits = c(0.2,1))+ 
  theme_bw(base_size = 18)+
  xlab("Enrichment score") + ylab("") +
  theme(legend.position=c(1,0),legend.justification = c(1,0))+

  theme(legend.background = element_blank())+
  labs(colour = "Enrichment Score", size = "-Log10 (P-value)") +
  theme(legend.key = element_blank())

ggsave("Output/Figure/GSEA_bubbles.pdf", width = 10, height = 12)
ggsave("Output/Figure/GSEA_bubbles.png", width = 12, height = 10)
```

![image](https://github.com/ying-ge/clbn/Output/Figure/GSEA_bubbles.png)

## :five: Heatmap for mitochondial genes from Dong Li

```{r}
library(pheatmap)
selectgenesID <- read.table("Input/mitoch_geneID.txt")

eset_DEG <- read.csv("Output/Table/DEG_signal.csv", row.names = 1)
eset_selectgenes <- eset_DEG[row.names(eset_DEG) %in% selectgenesID$V1, ]

setdiff(selectgenesID$V1, rownames(eset_selectgenes))
selectgenesID <- union(selectgenesID$V1, rownames(eset_selectgenes))

eset_uniq <- read.csv("Output/Table/signal_gsym.csv", row.names = 1)
eset_signaling <- eset_uniq[row.names(eset_uniq) %in% selectgenesID, ]
dim(eset_signaling)

annotation_col <- data.frame(strain = c(rep("WT",3),rep("MU",3)))
rownames(annotation_col) = colnames(eset_signaling)

ann_colors = list(
  strain = c(WT = "red", MU = "blue")
)

pheatmap <- pheatmap(eset_signaling,cellwidth = 20, cellheight = 12, fontsize = 10,
         method="spearman", 
         scale="row", 
         cluster_rows=T,
         cluster_cols=T,
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         #treeheight_row = "0",
         #treeheight_col = "0",
         border_color = "NA",
         silent=T, 
         filename = "Output/Figure/heatmap_mitoch_scale.png")
```

![image](https://github.com/ying-ge/clbn/Output/Figure/heatmap_mitoch_scale.png)

## :six: Draw volcano plot with mitochondial genes

```{r, fig.width=6, fig.height=4}
library(ggrepel)

allgeneslimma <- read.csv("Output/Table/limma_output.csv", row.names = 1)

# differentially expressed genes
allgeneslimma$threshold <- as.factor(
  ifelse(allgeneslimma$adj.P.Val <= 0.05 & abs(allgeneslimma$logFC) >= 1, 
         ifelse(allgeneslimma$logFC >= 1, "Increase", "Decrease"), "Ambiguous"))
table(allgeneslimma$threshold)

allgeneslimma$label <- row.names(allgeneslimma)

# mitochondial genes
selectgenes <- allgeneslimma[row.names(allgeneslimma) %in% selectgenesID,]

p <- ggplot(allgeneslimma, aes(logFC, -log10(P.Value), color = threshold, label = label)) + 
  geom_point(alpha = 0.3, size = 2, shape = 16) +
  scale_color_manual(values = c("lightgrey", "blue", "red"), guide = FALSE) +
  
  geom_vline(xintercept = c(-1, 1),lty = 2,col = "darkgrey", lwd = 0.6) +
  geom_hline(yintercept = -log10(0.05), lty = 2,col = "darkgrey", lwd = 0.6) +

  # DEG
  geom_point(data = subset(allgeneslimma, 
                           threshold == "Increase" | threshold == "Decrease"), 
             alpha = 1, size = 4, shape = 16) +
  
  # mitochondial genes
  geom_point(data = selectgenes, alpha = 1, size = 4.2, shape = 1, 
             stroke = 1, 
             color = "black") +

  ylab(bquote(~-Log[10]~adjusted~italic("P-value"))) +
  xlab(bquote(~Log[2]~"(fold change)")) +
  
  scale_x_continuous(
    breaks = c(-6, -3, -1, 0, 1, 3, 6), 
    labels = c(-6, -3, -1, 0, 1, 3, 6),
    limits = c(-5.9, 5.9) 
  )  + 
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank())


# label mitochondial genes
set.seed(123)
pvol <- p + geom_text_repel(data = subset(selectgenes, logFC> 0),
                          nudge_x = 3.5 - subset(selectgenes, logFC> 0)$logFC,
                          segment.size = 0.2, 
                          hjust = 0,
                          direction = "y",
                          colour="darkred",
                          size = 4,
                          show.legend = FALSE) +
  geom_text_repel(data = subset(selectgenes, logFC < 0),
                          nudge_x = -5 - subset(selectgenes, logFC < 0)$logFC,
                          segment.size = 0.2, 
                          hjust = 0,
                          direction = "y",
                          colour="darkred",
                          size = 4,
                          show.legend = FALSE)
pvol
ggsave("Output/Figure/valcono.pdf", width = 6, height = 4)
ggsave("Output/Figure/valcono.png", width = 6, height = 4)
```

![image](https://github.com/ying-ge/clbn/Output/Figure/valcono.png)

## :+1: Citation

If the code is helpful to you, please kindly cite our paper as:  

    @article{10.1371/journal.pgen.1009140,
    author = {Dai, Zhaoxia AND Li, Dong AND Du, Xiao AND Ge, Ying AND Hursh, Deborah A. AND Bi, Xiaolin},
    journal = {PLOS Genetics},
    publisher = {Public Library of Science},
    title = {Drosophila Caliban preserves intestinal homeostasis and lifespan through regulating mitochondrial dynamics and redox state in enterocytes},
    year = {2020},
    month = {10},
    volume = {16},
    url = {https://doi.org/10.1371/journal.pgen.1009140},
    pages = {1-26},
    abstract = {Author summary Self-renewal and differentiation of somatic stem cells are critical for tissue homeostasis. In Drosophila, intestinal homeostasis is maintained by tightly controlled proliferation and differentiation of intestinal stem cells (ISCs). In this study, we demonstrate that maintenance of mitochondrial dynamics and redox balance in enterocytes (ECs) by Caliban (Clbn) is critical for ISCs proliferation and intestinal homeostasis. We show that clbn mutant flies have shortened lifespan, involving disruption of intestinal homeostasis. We find that Clbn is highly expressed and localized to the outer membrane of the mitochondria in enterocytes. Clbn is important for mitochondrial dynamics, as loss of clbn in enterocytes results in mitochondrial fragmentation. Mitochondrial defects cause accumulation of reactive oxygen species (ROS) production and cellular damage, which in turn leads to activation of oxidative stress and promotion of ISCs over-proliferation. We further show that depletion of clbn promotes tumor growth in gut generated by activated Ras in intestinal progenitor cells. Our results establish Clbn as a modulator of mitochondrial dynamics in enterocytes, and highlight the importance of mitochondrial dynamics in the regulation of somatic stem cells activity and tissue homeostasis.},
    number = {10},
    doi = {10.1371/journal.pgen.1009140}
    }
