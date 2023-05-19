
#Question 1
#########################################################################
#########################################################################

library(Biobase)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE14333", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#Ensuring the column names are valid
colnames(gset) <- make.names(colnames(gset))

# group membership for all samples
gsms <- paste0("0X111011000010101XX010XX11100X11X0000111X01X11X0X1",
               "10010110101100010100101001011101001100101100110X11",
               "X00X0000110X000100001X1XXX1X00X10X1100000000011110",
               "1X1001101110011101001001111010111001110X11X010010X",
               "10011110XX00X0X00X0X100110010X11110111101X0001X00X",
               "10000000X10100110011111011X1X101XX101X0X")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X"), samples such as rectal etc.
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]


# extract gene expression data
ex <- exprs(gset)

# check if log2 transformation is necessary based on data distribution
qx <- quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

# if necessary, perform log2 transformation
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  ex <- log2(ex)
  exprs(gset) <- ex
}

#Q2
#########################################################################
#########################################################################

# check for missing values
sum(is.na(ex))

# check for outliers
boxplot(ex)

library(limma)
ex_norm <- normalizeBetweenArrays(ex, method="quantile")
# Create an ExpressionSet object
eset <- ExpressionSet(assayData = ex_norm)

# Extract data attributes
pdata <- pData(eset)
fdata <- fData(eset)

# List data attributes
attributes(eset)
attributes(pdata)
attributes(fdata)

#Q4
#########################################################################
#########################################################################


# assign samples to groups 
gs <- factor(sml)
groups <- make.names(c("Right","Left"))
levels(gs) <- groups
gset$group <- gs

# Load required libraries
library(stats)
library(ggplot2)

# perform t-test and store results
ttest <- apply(ex_norm, 1, function(x) t.test(x[groups=="Right"], x[groups=="Left"]))
ttest.df <- data.frame(
  gene = rownames(ex_norm),
  pvalue = sapply(ttest, function(x) x$p.value),
  estimate = sapply(ttest, function(x) x$estimate[1]-x$estimate[2])
)

# correct p-values with Holm method
ttest.df$padj <- p.adjust(ttest.df$pvalue, method = "holm")

# create volcano plot with colored points
ggplot(ttest.df, aes(x=estimate, y=-log10(pvalue), color=padj < 0.05)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("blue", "red")) +
  theme_classic() +
  labs(x="log2(Fold Change)", y="-log10(p-value)") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", alpha=0.5)


#Q5
#########################################################################
#########################################################################

# set up design matrix
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="holm", sort.by="B", number=250)



tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# remove outliers
upper_fc <- quantile(tT$logFC, 0.99)
lower_fc <- quantile(tT$logFC, 0.01)
upper_padj <- quantile(tT$adj.P.Val, 0.99)
tT_filtered <- tT[tT$logFC <= upper_fc & tT$logFC >= lower_fc & tT$adj.P.Val <= upper_padj, ]

# create volcano plot with filtered data
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)),
            xlim=c(-3,3), ylim=c(0,5))

################################################################

#Q7
#########################################################################
#########################################################################

library(clusterProfiler)

# get the vector of gene symbols from the topTable output
gene_symbol <- tT$Gene.symbol

# set the background gene set to be all genes in the ExpressionSet
library(org.Hs.eg.db)
# perform gene set enrichment analysis using the "org.Hs.eg.db" database
ego <- enrichGO(gene = gene_symbol, 
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod= "holm",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.1,
                      readable = TRUE)

# print the results
ego@result

#Q8
#########################################################################
#########################################################################

###Dot plot 
library(ggplot2)
library(dplyr)

# Subset the top 20 pathways
top_pathways <- head(ego@result, 60)

# Plot the top 20 enriched pathways using ggplot2
g <- ggplot(top_pathways, aes(x = -log10(pvalue), y = Description))
g + geom_point(size = 3) + xlab("-log10(p-value)") + ylab("Pathway")


###Bar plot
library(ggplot2)
library(dplyr)
# Get the top 20 enriched pathways
top_pathways <- head(ego@result, 20)

# Create a bar plot of the enrichment scores
ggplot(top_pathways, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "#69b3a2") +
  xlab("Pathway") + ylab("-log10(p-value)") +
  ggtitle("Enrichment Scores") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


#Q9
#########################################################################
#########################################################################

library(enrichR)

# get the vector of gene symbols from the topTable output
gene_symbols <- tT$Gene.symbol


# Perform KEGG pathway analysis using enrichKEGG function
kegg_enrich <- enrichr(gene = gene_symbols, database= "KEGG_2019_Human")

# print the results
kegg_enrich

summary(kegg_enrich)

#top pathways from the analysis
main_path <- head(ego@result, n=15)

paths <- main_path$Description

print(paths)
