# microbiome

---
title: "Toadsih-DADA2"
author: "Anthony Bonacolta"
date: "7/15/22"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{bash}
ls *_L001_R1_001.fastq.gz | cut -f 1-2 -d "_" > samples
for sample in $(cat samples)
do
echo "On sample: $sample"
cutadapt -a ^GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC \
-A ^GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC \
-m 200:200 -M 300:300 --discard-untrimmed \
-o ${sample}_L001_R1_001_trimmed.fastq.gz -p ${sample}_L001_R2_001_trimmed.fastq.gz \
${sample}_L001_R1_001.fastq.gz ${sample}_L001_R2_001.fastq.gz >> cutadapt_primer_trimming_stats.txt 2>&1
done
```

```{bash}
paste samples <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")
```


```{r}
library(dada2); packageVersion("dada2")
```

```{r}
setwd("~/Desktop/Toadfish_workflow/data/BonacoltaEMPV4")
path <- "~/Desktop/Toadfish_workflow/data/BonacoltaEMPV4" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```
```{r}
fnFs <- sort(list.files(path, pattern="_L001_R1_001_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001_trimmed.fastq.gz", full.names = TRUE))
```

```{r}
sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
```

```{r}
plotQualityProfile(fnFs[2:4])
plotQualityProfile(fnRs[2:4])
```


```{r}
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=10, truncLen=c(250,250),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```


```{r}
errF <- learnErrors(filtFs, multithread=2)
errR <- learnErrors(filtRs, multithread=2)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```

```{r}
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names
ddFs <- dada(derepFs, err=NULL, selfConsist=TRUE)
ddRs <- dada(derepRs, err=NULL, selfConsist=TRUE)
plotErrors(ddFs)
plotErrors(ddRs)
dadaFs <- dada(derepFs, err=ddFs[[1]]$err_out, pool=TRUE, multithread=2)
dadaRs <- dada(derepRs, err=ddRs[[1]]$err_out, pool=TRUE, multithread=2)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
```


```{r}
seqtab.all <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab.all)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths
seqtab.chim <- getSequences(seqtab.all)[!getSequences(seqtab.all) %in% getSequences(seqtab)]
```

```{r}
ref_fasta <- "~/Desktop/Toadfish_workflow/silva_nr99_v138_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab, refFasta=ref_fasta, multithread=2)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxa <- addSpecies(taxa, "~/Desktop/Toadfish_workflow/silva_species_assignment_v138.fa.gz")
taxa.print <- taxa
```

```{r}
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(ALDEx2);packageVersion("ALDEx2")
library(CoDaSeq);packageVersion("CoDaSeq")
library(vegan); packageVersion("vegan")
library(viridis); packageVersion("viridis")
library(tidyr); packageVersion("tidyr")
library("dplyr"); packageVersion("dplyr")
library("scales"); packageVersion("scales")
library("grid"); packageVersion("grid")
library("reshape2"); packageVersion("reshape2")
library("edgeR"); packageVersion("edgeR")
library("plyr"); packageVersion("plyr")
library("tidyr"); packageVersion("tidyr")
library("viridis"); packageVersion("viridis")
library("DESeq2"); packageVersion("DESeq2")
library("gridExtra"); packageVersion("gridExtra")
library("ampvis2"); packageVersion("ampvis2")

```

```{r}
setwd("~/Desktop/Toadfish_workflow/analysis")
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
 asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)

write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

asv_tabtax <- cbind(asv_tab,asv_tax)
row.names(asv_tabtax) <- sub(">", "", asv_headers)
write.table(asv_tabtax, "ASVs_counts_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
```

```{r}
setwd("~/Desktop/Toadfish_workflow/analysis")
SV <- read.table("ASVs_counts.tsv", row.names = 1, check.names= FALSE, sep = "\t", header = TRUE)
tax <-as.matrix(read.table("ASVs_taxonomy.tsv", row.names = 1, header = TRUE, sep = "\t"))
map <- read.table("toadfish_metadata.txt", row.names = 1, sep ="\t", header = TRUE)
ps = phyloseq(otu_table(SV, taxa_are_rows=TRUE), 
               sample_data(map), 
               tax_table(tax))
ps
```

```{r}
psf <- prune_samples(sample_sums(ps)>=250, ps)
f1<- filterfun_sample(topf(0.99))
wh1 <- genefilter_sample(psf, f1, A=2)
psf <- prune_taxa(wh1, psf)
psf_toadfish<- subset_samples(psf, psf@sam_data$Source=='Toadfish')
psf_toadfish<- subset_samples(psf, psf@sam_data$TANK!='QC')
psf_toadfish
```

```{r misc, include=FALSE}
setwd("~/Desktop/Toadfish_workflow/analysis")
phyla_csv <- read.csv("phyla.csv",header=F)
phyla <- phyla_csv %>% pull(V1)
classes_csv <- read.csv("class.csv",header=F)
classes <- classes_csv %>% pull(V1)
families_csv <- read.csv("family.csv",header=F)
families <- families_csv %>% pull(V1)
orders_csv <- read.csv("order.csv",header=F)
orders <- orders_csv %>% pull(V1)
```


```{r}
library(randomcoloR)
n <- 60
palette <- distinctColorPalette(n)
```

# Figures
## Bubble by Region & Heat-Resistance (Euk & Bac)
```{r}
type_taxa <- psf_toadfish %>% 
  tax_glom(taxrank = "Class") %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>% 
  arrange(Class) 

p <- ggplot(type_taxa, aes(x = Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=palette) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance \n") +
  ggtitle("Class Composition of Toadfish Samples") + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Sample")
p

pdf("~/Desktop/Toadfish_workflow/analysis/figures/ra_Class.pdf", height = 7, width=13)
plot(p)
dev.off()

```

```{r}
type_taxa <- psf_toadfish %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>% 
  arrange(Phylum) 

p <- ggplot(type_taxa, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=palette) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance \n") +
  ggtitle("Phylum Composition of Toadfish Samples") + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Sample")
p

pdf("~/Desktop/Toadfish_workflow/analysis/figures/ra_Phylum.pdf", height = 7, width=13)
plot(p)
dev.off()

```
```{r}
type_taxa <- psf_toadfish %>% 
  tax_glom(taxrank = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>% 
  arrange(Order) 

p <- ggplot(type_taxa, aes(x = Sample, y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=palette) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance \n") +
  ggtitle("Order Composition of Toadfish Samples") + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Sample")
p

pdf("~/Desktop/Toadfish_workflow/analysis/figures/ra_Order.pdf", height = 7, width=15)
plot(p)
dev.off()
```

```{r}
type_taxa <- psf_toadfish %>% 
  tax_glom(taxrank = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)*100} )  
type_taxa <- type_taxa %>% 
  psmelt()  %>%
	arrange(Location)

sample_names(psf_toadfish)

xx <- ggplot(type_taxa, aes(x = factor(Sample, levels=c("1-9A-Ant", "26-9C-Ant", "15-9B-Ant", "13-35A-Ant", "25-35C-Ant", "14-35B-Ant",  "12-60C-Ant", "24-60B-Ant", "9-60A-Ant", "11-QC-Post", "23-9A-Post", "16-9B-Post", "18-35B-Post", "20-35C-Post", "5-36A-Post", "22-60B-Post","3-60A-Post", "7-60C-Post", "17-35B-Prec", "19-35A-Prec", "2-60C-Prec", "21-60A-Prec", "4-60B-Prec", "43-9B-Fluid", "6-9A-Fluid", "27-9C-Fluid", "40-35C-Fluid", "42-35A-Fluid","10-36B-Fluid", "8-60B-Fluid")), y = factor(Order, levels=orders))) + geom_point(aes(size = Abundance, fill=Phylum), alpha = 0.9, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "Sample", y = "Order", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Orders within toadfish") 
xx


pdf("~/Desktop/Toadfish_workflow/analysis/figures/bub_orders.pdf", height =11, width=10)
plot(xx)
dev.off()
```

```{r}
type_taxa <- psf_toadfish %>% 
  tax_glom(taxrank = "Class") %>% 
  transform_sample_counts(function(x) {x/sum(x)*100} )  
type_taxa <- type_taxa %>% 
  psmelt()  %>%
	arrange(Location)

sample_names(psf_toadfish)

xx <- ggplot(type_taxa, aes(x = factor(Sample, levels=c("1-9A-Ant", "26-9C-Ant", "15-9B-Ant", "13-35A-Ant", "25-35C-Ant", "14-35B-Ant",  "12-60C-Ant", "24-60B-Ant", "9-60A-Ant", "11-QC-Post", "23-9A-Post", "16-9B-Post", "18-35B-Post", "20-35C-Post", "5-36A-Post", "22-60B-Post","3-60A-Post", "7-60C-Post", "17-35B-Prec", "19-35A-Prec", "2-60C-Prec", "21-60A-Prec", "4-60B-Prec", "43-9B-Fluid", "6-9A-Fluid", "27-9C-Fluid", "40-35C-Fluid", "42-35A-Fluid","10-36B-Fluid", "8-60B-Fluid")), y = factor(Class, levels=classes))) + geom_point(aes(size = Abundance, fill=Phylum), alpha = 0.9, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "Sample", y = "Class", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Classes within toadfish") 
xx


pdf("~/Desktop/Toadfish_workflow/analysis/figures/bub_class.pdf", height =7, width=10)
plot(xx)
dev.off()
```


```{r}
bub <- merge_samples(psf_toadfish, "Treatment") 
type_taxa <- bub %>% 
  tax_glom(taxrank = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)*100} )  
type_taxa <- type_taxa %>% 
  psmelt()  %>%
	arrange(Location)

xx <- ggplot(type_taxa, aes(x = factor(Sample, levels=c("QC", "9", "35", "36", "60")), y = factor(Order, levels=orders))) + geom_point(aes(size = Abundance, fill=Phylum), alpha = 0.9, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "Treatment", y = "Order", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Orders between toadfish treatments") 
xx


pdf("~/Desktop/Toadfish_workflow/analysis/figures/bub_order_treatment.pdf", height =7, width=10)
plot(xx)
dev.off()
```

```{r}
bub <- merge_samples(psf_toadfish, "Location") 
type_taxa <- bub %>% 
  tax_glom(taxrank = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)*100} )  
type_taxa <- type_taxa %>% 
  psmelt()  %>%
	arrange(Location)

xx <- ggplot(type_taxa, aes(x = factor(Sample, levels=c("Anterior", "Posterior", "Precipitates", "Fluid", "Water")), y = factor(Order, levels=orders))) + geom_point(aes(size = Abundance, fill=Phylum), alpha = 0.9, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "Sample Location", y = "Order", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Orders between toadfish location") 
xx


pdf("~/Desktop/Toadfish_workflow/analysis/figures/bub_order_location.pdf", height =7, width=10)
plot(xx)
dev.off()
```

```{r}
variable1 = as.character(get_variable(psf_toadfish, "Treatment"))
variable2 = as.character(get_variable(psf_toadfish, "Location"))
sample_data(psf_toadfish)$NewPastedVar <- mapply(paste0, variable1, variable2, 
    collapse = "_")

bub <- merge_samples(psf_toadfish, "NewPastedVar")
type_taxa <- bub %>% 
  tax_glom(taxrank = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)*100} )

type_taxa <- type_taxa %>% 
  psmelt() %>% 
  arrange(Abundance) 

xx <- ggplot(type_taxa["Abundance" > 0], aes(x = factor(Sample, levels=c("9Anterior", "9Posterior", "9Fluid", "9Precipitates", "35Anterior", "35Posterior", "35Fluid", "35Precipitates", "60Anterior", "60Posterior", "60Fluid", "60Precipitates")), y = factor(Order, levels=orders))) + geom_point(aes(size = Abundance, fill=Phylum), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "Sample", y = "Order", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Orders between toadfish locations & treatments")
xx

pdf("~/Desktop/Toadfish_workflow/analysis/figures/bub_order_location_treatment.pdf", height =10, width=10)
plot(xx)
dev.off()
```
# Beta-diversity
```{r}
ps_clr <- microbiome::transform(psf_toadfish, 'clr', shift = 1)
psr_clr.ord <- ordinate(ps_clr, "RDA", "euclidean")
PCA = plot_ordination(ps_clr, psr_clr.ord,
                                shape="Treatment",
                                color="Location",
                                title="Aitchison Distance PCA of Toadfish Prokaryome") + theme_classic() + geom_point(aes(color = Location), alpha = 1, size = 3) + scale_colour_manual(values=palette) + labs(color="Location") + stat_ellipse()

PCA

pdf("~/Desktop/Toadfish_workflow/analysis/figures/PCA_Location_Treatment.pdf",  width=10, height=6.5)
plot(PCA)
dev.off()
```



# Alpha-Diversity
```{r}
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)


variable1 = as.character(get_variable(psf_toadfish, "Treatment"))
variable2 = as.character(get_variable(psf_toadfish, "Location"))
sample_data(psf_toadfish)$NewPastedVar <- mapply(paste0, variable1, variable2, 
    collapse = "_")

ad <- prune_taxa(taxa_sums(psf_toadfish) > 0, psf_toadfish) #Use whichever phyloseq object you want here.  
tab <- microbiome::alpha(ad, index = "all")
kable(head(tab))
ps1.meta <- meta(ad)
kable(head(ps1.meta))
ps1.meta$Shannon <- tab$diversity_shannon # We want to look at the shannon-weiner diversity index.

p1 <- ggboxplot(ps1.meta, x = "NewPastedVar", y = "Shannon",
 add = "boxplot", fill = "Location", order=c("9Anterior", "9Posterior", "9Fluid", "9Precipitates", "35Anterior", "35Posterior", "35Fluid", "35Precipitates", "60Anterior", "60Posterior", "60Fluid", "60Precipitates")) + theme_classic() + scale_fill_viridis_d()

pdf("~/Desktop/Toadfish_workflow/analysis/figures/alpha.pdf",  width=12, height=5)
plot(p1)
dev.off()
```
