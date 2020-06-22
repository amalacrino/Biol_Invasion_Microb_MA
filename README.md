# Biological invasions alter environmental microbiomes: a meta-analysis  

### *Antonino Malacrinò, Victoria A. Sadowski, Tvisha K. Martin, Nathalia Cavichiolli de Oliveira, Ian J. Brackett, James D. Feller, Kristian J. Harris, Orlando Combita Heredia, Rosa Vescio, Alison E. Bennett*

#### Journal, 20xx. Preprint DOI: [10.1101/2020.06.17.157883](https://www.biorxiv.org/content/10.1101/2020.06.17.157883v1) Published article DOI: XXX

## Acknowledgements

This study is the result of the joint effort of the SP2020 class EEOB 8896.12 "Microbiome meta-analysis" held at The Ohio State University. Along with the pure ecological aspects, we focused our attention on the concept of Open Science and its positive impact on science and society. From this our decision to:

1. use publicly available data to test our hypothesis;
2. use open source tools for data analysis;
3. share the scripts we coded for data analysis;
4. share a pre-print version of our manuscript;
5. submit our manuscript to journal managed by a gold open-access publisher with open peer-review.

## Abstract 

Biological invasions impact both agricultural and natural systems. The damage can be quantified in terms of both economic loss and reduction of biodiversity. Although the literature is quite rich about the impact of invasive species on plant and animal communities, their impact on environmental microbiomes is underexplored. Here, we re-analyze publicly available data using a common framework to create a global synthesis of the effects of biological invasions on environmental microbial communities. Our findings suggest that non-native species are responsible for the loss of microbial diversity and shifts in the structure of microbial populations. Therefore, the impact of biological invasions on native ecosystems might be more pervasive than previously thought, influencing both macro- and micro-biomes. We also identified gaps in the literature which encourage research on a wider variety of environments and invaders, and the influence of invaders across seasons and geographical ranges.  


## 0. Data

| **Study  ID** | Link to raw data                                             | **Reference**          | DOI                          |
| ------------- | ------------------------------------------------------------ | ---------------------- | ---------------------------- |
| MPG13011      | http://metagenomics.anl.gov/linkin.cgi?project=13011         | Gibbons  et al. 2017   | 10.1128/mSystems.00178-16    |
| MPG87547      | https://www.mg-rast.org/mgmain.html?mgpage=project&project=mgp87547 | Wehr et  al. 2019      | 10.1038/s41598-019-48922-    |
| PRJNA296487   | https://www.ncbi.nlm.nih.gov/bioproject/PRJNA296487          | Rodrigues  et al. 2015 | 10.1371/journal.pone.0141424 |
| PRJNA320310   | https://www.ncbi.nlm.nih.gov/bioproject/PRJNA320310/         | Collins  et al. 2016   | 10.1111/1365-2745.12616      |
| PRJNA385848   | https://www.ncbi.nlm.nih.gov/bioproject/PRJNA385848          | Denef et  al. 2017     | 10.1128/mSphere.00189-17     |


## 1. Download data

Data from NCBI SRA was downloaded using SRA Toolkit 2.10.4. 

```bash
prefetch --option-file SRA_Acc_List.txt
for i in *; do fasterq-dump --split-3 $i; done
find . ! -name '*.fastq' -delete
```

data from MG-RAST was directly downloaded from the website.

## 2. Merge PE reads

PE reads were merged using FLASH 1.2.11 (Magoč and Salzberg, 2011)

```bash
for i in *_1.fastq; do FILENAME=`basename $i _1.fastq`; echo $FILENAME; flash ${FILENAME}_1.fastq ${FILENAME}_2.fastq --output-prefix=${FILENAME}.merged; done

find . ! -name '*.merged.extendedFrags.fastq' -delete
for file in *; do mv ${file} ${file/.merged.extendedFrags/}; done
```

## 3. Process data with QIIME

Data handling was carried out using QIIME 1.9.1 (Caporaso et al. 2010, Caporaso et al. 2012). Taxonomy was assigned to each OTU through the BLAST method using SILVA v132 database (Quast et al. 2012). 

Get SILVA v132 reference database:

```bash
mkdir $PRJ_DIR/Taxonomy
cd $PRJ_DIR/Taxonomy
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
unzip Silva_132_release.zip
rm -r Silva_132_release.zip
cd $PRJ_DIR
```

Run the pipeline:

```bash
mkdir $PRJ_DIR/1_Split_Libraries
multiple_split_libraries_fastq.py -i $PRJ_DIR/0_Data -o $PRJ_DIR/1_Split_Libraries --sampleid_indicator .fastq

mkdir $PRJ_DIR/2_Pick_OTUs
pick_open_reference_otus.py -i $PRJ_DIR/1_Split_Libraries/seqs.fna -o $PRJ_DIR/2_Pick_OTUs -f -a -O 32 --suppress_step4 --suppress_taxonomy_assignment --suppress_align_and_tree

mkdir $PRJ_DIR/3_Pick_Repset
pick_rep_set.py -i $PRJ_DIR/2_Pick_OTUs/final_otu_map.txt -f $PRJ_DIR/1_Split_Libraries/seqs.fna -m most_abundant -o $PRJ_DIR/3_Pick_Repset/repset.fasta

mkdir $PRJ_DIR/4_Assign_Taxonomy
parallel_assign_taxonomy_blast.py -i $PRJ_DIR/3_Pick_Repset/repset.fasta -o $PRJ_DIR/4_Assign_Taxonomy -O 32 -r $PRJ_DIR/Taxonomy/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -t $PRJ_DIR/Taxonomy/SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt

mkdir $PRJ_DIR/QIIME1/5_OTU_Table
make_otu_table.py -i $PRJ_DIR/2_Pick_OTUs/final_otu_map.txt -o $PRJ_DIR/5_OTU_Table/otutable.biom -t $PRJ_DIR/4_Assign_Taxonomy/repset_tax_assignments.txt

mkdir $PRJ_DIR/QIIME1/6_Phylogeny
mafft --thread 32 $PRJ_DIR/3_Pick_Repset/repset.fasta > $PRJ_DIR/6_Phylogeny/rep_set_aligned.fasta
make_phylogeny.py -i $PRJ_DIR/6_Phylogeny/rep_set_aligned.fasta -o $PRJ_DIR/6_Phylogeny/tree.tre
```

For the data analysis we used:

```bash
$PRJ_DIR/5_OTU_Table/otutable.biom 
$PRJ_DIR/6_Phylogeny/tree.tre
```

## 4. Data analysis

Load packages

```R
library("phyloseq") #McMurdie & Holmes (2013) 
library("ggplot2") #Wickham (2011)
library("DESeq2") #Love et al. (2014)
library("data.table") #https://cran.r-project.org/web/packages/data.table/index.html
library("limma") #Smyth (2005)
library("ggsignif") #Ahlmann-Eltze (2017)
library("ggpubr") #Kassambara (2017)
library("car") #Fox et al. (2012)
library("lme4") #Bates et al. (2007)
library("scales") #https://cran.r-project.org/web/packages/scales/index.html
```

Load data

```R
biom <- import_biom(BIOMfilename ="otutable.biom")
map <- import_qiime_sample_data("metadata.tsv")
tree <- read.tree("tree.tre")
tree <- root(tree, 1, resolve.root = T)
MicroDF <- merge_phyloseq(biom, map, tree)
```

Format taxonomy

```R
tax_table(MicroDF) <- tax_table(MicroDF)[,c(1:7)]
tax_table(MicroDF) <- gsub(".*_","",tax_table(MicroDF))
colnames(tax_table(MicroDF)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
```

Clean database

```R
MicroDF <- subset_taxa(MicroDF, Class !="Chloroplast") 
MicroDF <- subset_taxa(MicroDF, Order !="Chloroplast")
MicroDF <- filter_taxa(MicroDF, function (x) {sum(x > 0) > 1}, prune=TRUE)
MicroDF <- prune_samples(sample_sums(MicroDF) >= 5000, MicroDF)
```

Normalize data

```R
diagdds = phyloseq_to_deseq2(MicroDF, ~ Sample_type)
gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
diagdds.c <- removeBatchEffect(diagvst)
diagdds.c[diagdds.c<0] <- 0
MicroDF2 <- MicroDF
otu_table(MicroDF2) <- otu_table(diagdds.c, taxa_are_rows = TRUE)
```

PERMANOVA

```R
sampledf <- data.frame(sample_data(MicroDF2))
dist.mat <- phyloseq::distance(MicroDF2, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~ Sample_type * Organism, strata = c("Study_ID", "Environment", "Invasive_species"), data = sampledf, permutations = perm)
pmv
```

RDA ordination

```R
RDA_ord <- ordinate(physeq = MicroDF2, method = "RDA", distance = dist.mat, formula = ~ Sample_type)

RDA_plot <- plot_ordination(physeq = MicroDF2, ordination = RDA_ord, axes = c(1,2)) +
  theme_bw(base_size = 14) +
  stat_ellipse(aes(fill = Sample_type), alpha = 0.4, geom = "polygon", colour = "black", size=0.3, show.legend=T) +
  geom_point(aes(colour=Sample_type, fill=Sample_type), shape = 21, size = 5, colour = "black") +
  theme(legend.title=element_blank(), 
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.01),
        legend.justification = c(0.99, 0.01)) +
  scale_fill_manual(name = "Legend", values=c("#0570b0", "#d7301f"), breaks = c("Control", "Invaded")) +
  scale_colour_manual(name = "Legend", values=c("#0570b0", "#d7301f"), breaks = c("Control", "Invaded"))
RDA_plot
```

Phylogenetic diversity

```R
diversity <- estimate_richness(MicroDF, split = TRUE, measures = c("Shannon"))
diversity <- cbind(sample_data(MicroDF), diversity)

model <- lmer(Shannon ~ Sample_type * Organism * (1|Study_ID) * (1|Environment), data = diversity)
Anova(model)

div_plot <- ggplot(diversity, aes(x = Sample_type, y = Shannon, fill = Sample_type)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(y = "Shannon diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="none") +
  scale_fill_manual(name = "Legend", values=c("#0570b0", "#d7301f"), breaks = c("Control", "Invaded")) +
  scale_colour_manual(name = "Legend", values=c("#0570b0", "#d7301f"), breaks = c("Control", "Invaded")) +
  geom_signif(comparisons = list(c("Control", "Invaded")), annotations="*", y_position = 510, tip_length = 0.03)
div_plot
```

Relative abundance

```R
phy <- transform_sample_counts(MicroDF2, function(x) x/sum(x))
glom <- tax_glom(phy, taxrank = 'Class')
dat <- data.table(psmelt(glom))
dat$Class <- as.character(dat$Class)
dat[, mean := mean(Abundance, na.rm = TRUE), by = "Class"]
dat[(mean <= 0.01), Class := "Others"]
dat <- dat[which(dat$Class !="Others")]
                               
boxplot.genus <- ggplot(dat[Abundance > 0], aes(x=Class, y=Abundance, fill = Sample_type)) + 
  geom_boxplot(outlier.shape=NA) + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(0.0001, 1)) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y=element_blank(),
        axis.text.y =element_text(face="italic"),
        legend.title=element_blank(), 
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="transparent"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.01),
        legend.justification = c(0.99, 0.01)) +
  ylab("Abundance (log scale)") +
  scale_x_discrete(limits = rev(unique(sort(dat$Class)))) +
  scale_fill_manual(name = "Legend", values=c("#0570b0", "#d7301f"), breaks = c("Control", "Invaded")) +
  scale_colour_manual(name = "Legend", values=c("#0570b0", "#d7301f"), breaks = c("Control", "Invaded"))
boxplot.genus

list.bact <-unique(c(as.character(dat$Class)))
model_calculator <- sapply(list.bact,  
                           function(x){
                             data.s <- dat[which(dat$Class==x),]
                             model <- lmer(Abundance ~ Sample_type * (1|Study_ID) * (1|Environment) * (1|Invasive_species), data = data.s)
                             res <-  Anova(model)
                             return(res)},
                           simplify = FALSE,USE.NAMES = TRUE)

res <- do.call(rbind, model_calculator)
res <- setDT(res, keep.rownames = TRUE)[]
```

Plot Figure 1

```R
p1 <- ggarrange(div_plot, RDA_plot, ncol = 2, nrow = 1,  align = "hv", widths = c(0.5, 1), heights = 1, labels = c("A", "B"))
p2 <- ggarrange(p1, boxplot.genus, ncol = 1, nrow = 2, widths = c(1, 1), heights = 1, labels = c("", "C"))
p2

ggsave(p2, filename = "Plots.pdf", dpi = 600, device = cairo_pdf, width = 7, height = 6.5, units = "in", family="Arial Unicode MS")

```

## References

Ahlmann-Eltze, C. (2017). ggsignif: Significance Brackets for ‘ggplot2’. *R package version 0.4. 0*.

Bates, D., Sarkar, D., Bates, M. D., & Matrix, L. (2007). The lme4 package. *R package version*, *2*(1), 74.

Caporaso, J. G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F. D., Costello, E. K., ... & Huttley, G. A. (2010). QIIME allows analysis of high-throughput community sequencing data. *Nature methods*, *7*(5), 335.

Caporaso, J. G., Lauber, C. L., Walters, W. A., Berg-Lyons, D., Huntley, J., Fierer, N., ... & Gormley, N. (2012). Ultra-high-throughput microbial community analysis on the Illumina HiSeq and MiSeq platforms. *The ISME journal*, *6*(8), 1621-1624.

Fox, J., Weisberg, S., Adler, D., Bates, D., Baud-Bovy, G., Ellison, S., ... & Heiberger, R. (2012). Package ‘car’. *Vienna: R Foundation for Statistical Computing*.

Kassambara, A. (2017). ggpubr:“ggplot2” based publication ready plots. *R package version 0.1*, *6*.

Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H., Ackerly, D. D., ... & Webb, C. O. (2010). Picante: R tools for integrating phylogenies and ecology. *Bioinformatics*, *26*(11), 1463-1464.

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome biology*, *15*(12), 550.

Magoč, T., & Salzberg, S. L. (2011). FLASH: fast length adjustment of short reads to improve genome assemblies. *Bioinformatics*, *27*(21), 2957-2963.

McMurdie, P. J., & Holmes, S.  (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PloS one*, *8*(4).

Quast, C., Pruesse, E., Yilmaz, P.,  Gerken, J., Schweer, T., Yarza, P., ... & Glöckner, F. O. (2012).  The SILVA ribosomal RNA gene database project: improved data processing  and web-based tools. *Nucleic acids research*, *41*(D1), D590-D596.

Smyth, G. K. (2005). Limma: linear models for microarray data. In *Bioinformatics and computational biology solutions using R and Bioconductor* (pp. 397-420). Springer, New York, NY.

Wickham, H. (2011). ggplot2. *Wiley Interdisciplinary Reviews: Computational Statistics*, *3*(2), 180-185.
