

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")

##if (!require("BiocManager", quietly = TRUE))
  ##install.packages("BiocManager")
##BiocManager::install()

## if (!requireNamespace("BiocManager", quietly = TRUE))
  ##install.packages("BiocManager")
## BiocManager::install("dada2", version = "3.10")

## Block 1 from Lee tutorial

library(dada2)
library(here)
library(DECIPHER)
packageVersion("dada2") # 1.11.5 when this was initially put together, though might be different in the binder or conda installation, that's ok!
library(decontam)
packageVersion("decontam") # 1.1.2 when this was put together
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(phyloseq) ; packageVersion("phyloseq") # 1.22.3
library(vegan) ; packageVersion("vegan") # 2.5.4
library(DESeq2) ; packageVersion("DESeq2") # 1.18.1
library(dendextend) ; packageVersion("dendextend") # 1.10.0
library(viridis) ; packageVersion("viridis") # 0.5.1

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

# Cram says I don't have to do this.
#setwd("~/dada2_amplicon_ex_workflow") Using data from Dr. Sosa's study

list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan(here("samples"), what = "character")

# one holding the file names of all the forward reads

## I don't have forward and reverse reads, so we're just doing this
## once and calling everything "forward" reads.

forward_reads <- paste0(samples, ".fastq.gz")
# # and one with the reverse
# reverse_reads <- paste0(samples, "_sub_R2.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, ".filtered.fq.gz")
#filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")

## Second block of code
plotQualityProfile(here("DATA", forward_reads))
#plotQualityProfile(reverse_reads)
# and just plotting the last 4 samples of the reverse reads
plotQualityProfile(here("DATA", forward_reads)[c(1)])

##Mike's third block of code' This step is in order to reduce the number of base pairs for the forward reads
filtered_out <- filterAndTrim(here ("DATA" , forward_reads) ,
                              here("Filtered" , filtered_forward_reads),
                              maxEE=c(2),
                              rm.phix=TRUE, minLen=175, truncLen=c(285))

## Mike's fourth block of code.This was to show how we were able to trim and filter reads down to smaller base pairs
plotQualityProfile(here("Filtered" , filtered_forward_reads[c(6:2)]))

## Mike's fifth block of code. This helps me find and uncertainity in my model/data set.
err_forward_reads <- learnErrors(here("Filtered" , filtered_forward_reads))

save(err_forward_reads, file = here("SavedRData", "err_forward_reads.RData"))

## Mike's sixth block of code.Allows me to see the error rates for transitioning A-C and A-G
plotErrors(err_forward_reads, nominalQ=TRUE)

## Mike's seventh block. Allows to see what my true biological data once all the filtering has been conducted above.
derep_forward <- derepFastq(here( "Filtered", filtered_forward_reads), verbose=TRUE)
names(derep_forward) <- samples

## Mike's 8 block 1. Was able to see how much time it took to process this command. But allowed for samples to filtered with higher intensity in order to not miss out on rare variants.
pt0 <- proc.time()
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo")
pt1 <- proc.time()
pt1- pt0

## Mike's 9th block of code. Gave us the number of observations per sample
seqtab <- makeSequenceTable(dada_forward)
class(seqtab)
dim(seqtab)

## Mike's 10 block of code. Removed any chimera data (data which is a combination of reads that will not be needed for this study )
## .92 alittle over 4000 bimers were found in this step
pt0 <- proc.time()
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
pt1 <- proc.time()
pt1- pt0
sum(seqtab.nochim)/sum(seqtab)

## Mike's 10 block of code. 
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples,
                          dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2],
                          dada_f=sapply(dada_forward, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
summary_tab
save(summary_tab, file = here( "SavedRData" , "summary_tab.RData"))

save(dada_forward, file = here( "SavedRData" , "dada_forward.RData"))
save(seqtab.nochim, file = here( "SavedRData", "seqtab.nochim.RData"))

## Mike's 11th block of code
write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)

## Mike's 12th block of code DECIPHER. Step was needed in order to run/assign taxonmy to my data set. 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("S4Vectors")

## Mike's 13th block of code. Ussed his training sets as a reference point for my own data. 

## loading reference taxonomy object
##load("dada_forward.RData")
## skipping this codeblock for time, and it will not run in the binder environment
## downloading DECIPHER-formatted SILVA v138 reference
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")

## loading DECIPHER
library(DECIPHER)
# packageVersion("DECIPHER") # v2.6.0 when this was initially put together, though might be different in the binder or conda installation, that's ok!

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

## Mike's 14th block of code. Loaded taxonomic info in order to find ASV for the next step. 
load("tax-info.RData")

# Mike's 15th block of code

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## Mike's 16th block of code

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

## Mike's 17th block of code

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)


## Mike's 18th block of code. Showed me the names of bacteria all the way down to the species. Help with identification
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)


save(tax_info, file = here("SavedRData", "tax_info.RData"))
save(asv_tab, file = here("SavedRData", "asv_tab.RData"))
save(dna, file = here("SavedRData", "dna.RData"))

# load("tax_info.RData")
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
# asv_tax <- t(sapply(tax_info, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa

colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

save(asv_tax, file = here("SavedRData", "asv_tax.RData"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

## mikes 20th block of code
colnames(asv_tab) # our blanks are the first 4 of 20 samples in this case
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 16))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 6 as contaminants

## don't worry if the numbers vary a little, this might happen due to different versions being used 
## from when this was initially put together

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

## mikes 21st block of code
asv_tax[row.names(asv_tax) %in% contam_asvs, ]

## mikes 22nd block of code


# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

save(contam_df, file = here("SavedRData", "contam_df.RData"))

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]
save(asv_fasta_no_contam, file = here("SavedRData", "asv_tab_no_contam.RData"))

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]
save(asv_tax_no_contam, file = here("SavedRData", " asv_tax_no_contam.RData"))


## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)


# don't worry if versions are different from what's listed here, shown are are just what was used when this was initially put together
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(phyloseq) ; packageVersion("phyloseq") # 1.22.3
library(vegan) ; packageVersion("vegan") # 2.5.4
library(DESeq2) ; packageVersion("DESeq2") # 1.18.1
library(dendextend) ; packageVersion("dendextend") # 1.10.0
library(viridis) ; packageVersion("viridis") # 0.5.1
library(ggplot2)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggplot2")

## NOTE ##
# if you loaded the saved R data above with `load("amplicon_dada2_ex.RData")` and
# didn't run all the steps, you may want to write out these files here before
# clearing the environment. It won't hurt if you did it already, so better to just
# run these here :)

write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)


# ok, now moving on

rm(list=ls())

df_count_table <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")[ , -c(1:4)]

save(count_tab, file = here("SavedRData", "count_tab.RData"))

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

save(tax_tab, file = here("SavedRData", "tax_tab.RData"))

plot(count_tab)

# sample_info_tab$color <- as.character(sample_info_tab$)
ample_info_tab$color <- as.character(samp$color)

count_tab$color <- as.character(count_tab$color)


###sample_info_tab <- read.table("sample_info_tab.tsv", header=T,row.names=1,
                              ##check.names=F, sep="\t")

sample_info_0 <- tibble(filename = colnames(count_tab))
sample_info_1 <- sample_info_0 %>%
  mutate(environment = str_extract(filename, "(?<=(on|in)_).*(?=\\.filtered)")) %>%
  mutate(color = case_when(
    environment == "glass" ~ "black",
    environment == "water" ~ "blue",
    environment == "Polystyrene" ~ "red",
    environment == "PLA" ~ "orange",
    environment == "Polypropyelene" ~ "pink",
    TRUE ~ "yellow" 
  ))

rarecurve(t(count_tab), col=sample_info_1$color, step=100, lwd=2, ylab= "ASVs", label=F)
legend(x=1, legend = colnames(count_tab))

set.seed(3)
grp<- factor(sample(seq_along(col),nrow(count_tab), replace=TRUE))
cols <- col[grp]

# first we need to make a DESeq2 object

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_1, design = ~environment) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. You could run:
#deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
#now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

euc_clust <- hclust(euc_dist, method="ward.D2")

plot(euc_clust)
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_1$color[order.dendrogram(euc_dend)])
dend_labs <- as.character(sample_info_1$environment[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols
labels(euc_dend) <- dend_labs

plot(euc_dend, ylab="VST Euc. dist.")
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
##<- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_1_phy = ~environment)
library("phyloseq")

phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="phylum")) 

phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="phylum"))[,"phylum"]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)

phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

# now we'll remove the Proteobacteria, so we can next add them back in
# broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

# making count table broken down by class (contains classes beyond the
# Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="class")) 

# making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="class")) 

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)

# making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$phylum == "Proteobacteria", "class"])

# changing the row names like above so that they correspond to the taxonomy,
# rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$class) 

# making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

# there are also possibly some some sequences that were resolved to the level
# of Proteobacteria, but not any further, and therefore would be missing from
# our class table
# we can find the sum of them by subtracting the proteo class count table
# from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

# now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)

# and to check we didn't miss any other sequences, we can compare the column
# sums to see if they are the same
# if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab)) 

# now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
# we see there are currently 42 rows, which might be a little busy for a
# summary figure
# many of these taxa make up a very small percentage, so we're going to
# filter some out
# this is a completely arbitrary decision solely to ease visualization and
# intepretation, entirely up to your data and you
# here, we'll only keep rows (taxa) that make up greater than 5% in any
# individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
# now we have 12, much more manageable for an overview figure

# though each of the filtered taxa made up less than 5% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

## don't worry if the numbers or taxonomy vary a little, this might happen due to different versions being used 
## from when this was initially put together


# first let's make a copy of our table that's safe for manipulating
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(filt_major_taxa_proportions_tab_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()

# take a look at the new table and compare it with the old one
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)
# manipulating tables like this is something you may need to do frequently in R

# now we want a table with "color" and "characteristics" of each sample to
# merge into our plotting table so we can use that more easily in our plotting
# function
# here we're making a new table by pulling what we want from the sample
# information table
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_1), "environment"=sample_info_1$environment, "color"=sample_info_1$environment, stringsAsFactors=F)
mutate(color = case_when(
  environment == "glass" ~ "black",
  environment == "water" ~ "blue",
  environment == "Polystyrene" ~ "red",
  environment == "PLA" ~ "orange",
  environment == "Polypropyelene" ~ "pink",
  TRUE ~ "yellow" 
))


# and here we are merging this table with the plotting table we just made
# (this is an awesome function!)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

# and now we're ready to make some summary figures with our wonderfully
# constructed table

## a good color scheme can be hard to find, i included the viridis package
## here because it's color-blind friendly and sometimes it's been really
## helpful for me, though this is not demonstrated in all of the following :/ 

# one common way to look at this is with stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")


ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(char), shape=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$color[order(filt_major_taxa_proportions_tab_for_plot.g2$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")


# let's set some helpful variables first:
Water_sample_IDs <- row.names(sample_info_1)[sample_info_1$type == "environment"]
PLA_sample_IDs <- row.names(sample_info_tab)[sample_info_1$type == "PLA"]

anova(betadisper(euc_dist, sample_info_1$environment))
