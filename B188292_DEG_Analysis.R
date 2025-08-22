#loading the necessary r packages
library(ggrepel)
library(gprofiler2)
library(edgeR)
library(UpSetR)
library(gridExtra)
library(fgsea)
library(plotly)
library(scatterplot3d)  
library(Rsubread)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(biomaRt)
library(msigdbr)
library(tidyverse)
library(plotly)
library(enrichR)
library(STRINGdb)
library(RUVSeq)
library(tseries)
library(patchwork)


################################################################################
#                                                                              #
#                            PART 1.A: Data Handling                             #
#                                                                              #
################################################################################

#1.a loading the prevously generated count data
load("GSE113343_features")
load("GSE1139957_features_2")

GSE113957_rna_seq<- (GSE1139957_features_2$counts)
GSE113343_rna_seq<- (GSE113343_features$counts)

GSE113957_meta <- read.table("GSE113957_meta_v2.csv", sep=',', row.names=1,fill=T,header=T)
GSE113343_meta <- read.table("GSE113343_meta.csv", sep=',', row.names=1,fill=T,header=T)


#1.b meta data handling
#Changing the values of "Normal" to "Control"
GSE113957_meta$disease[GSE113957_meta$disease == "Normal"] <- "Control"

#changing the "disease column in GSE113957 meta to Condition
colnames(GSE113957_meta)[colnames(GSE113957_meta) == "disease"] <- "Condition"
colnames(GSE113957_meta)[colnames(GSE113957_meta) == "sex"] <- "Sex"
colnames(GSE113343_meta)[colnames(GSE113343_meta) == "Disease_Status"] <- "Condition"

#adding the group column and assigning colours
group_colours <- c(Control_male = "#4A90E2",
                  HGPS_male = "#1F497D",
                  Control_female = "#FF8C42",
                  HGPS_female = "#B35E00")

GSE113957_meta$Group <- with(GSE113957_meta, paste(Condition, Sex, sep = "_"))
GSE113957_meta$Group_Colour <- group_colours[GSE113957_meta$Group]

GSE113343_meta$Group <- with(GSE113343_meta, paste(Condition, Sex, sep = "_"))
GSE113343_meta$Group_Colour <- group_colours[GSE113343_meta$Group]

#fixing the column names
colnames(GSE113957_rna_seq) <- c("SRR7093875", "SRR7093881", "SRR7093882", "SRR7093883",
                                 "SRR7093884", "SRR7093887", "SRR7093888", "SRR7093942",
                                 "SRR7093943", "SRR7093944", "SRR7093945", "SRR7093946",
                                 "SRR7093947", "SRR7093948", "SRR7093949", "SRR7093950",
                                 "SRR7093951")
colnames(GSE113343_rna_seq) <- c("SRR7030736", "SRR7030737", "SRR7030738", "SRR7030739")

#accessing column names
colnames_GSE113957 <- colnames(GSE113957_rna_seq)
colnames_GSE113343 <- colnames(GSE113343_rna_seq)

#Ordering the elements in colnames_rnaseq to the rownames of the metadata
GSE113957_idx<-match(colnames_GSE113957, rownames(GSE113957_meta))
GSE113343_idx<-match(colnames_GSE113343, rownames(GSE113343_meta))

#Renaming the column names of table_rnaseq to the Sample name from adf for conciseness
GSE113957_meta <- GSE113957_meta[GSE113957_idx,]
GSE113343_meta <- GSE113343_meta[GSE113343_idx,]

rownames(GSE113957_meta) <-colnames(GSE113957_rna_seq)
rownames(GSE113343_meta) <-colnames(GSE113343_rna_seq)

#Adding a Sample column 
GSE113957_meta$Sample <- rownames(GSE113957_meta)
GSE113343_meta$Sample <- rownames(GSE113343_meta)


#Because we are handling two data sets here, I am going to combine the two, and remove any genes that are not shared between the data sets. 
group_genes <- intersect(rownames(GSE113957_rna_seq), rownames(GSE113343_rna_seq))
GSE113957_data <- GSE113957_rna_seq[group_genes, , drop=FALSE]
GSE113343_data <- GSE113343_rna_seq[group_genes, , drop=FALSE]


################################################################################
#                                                                              #
#                            PART 2: DEG ANALYSIS                              #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
#                         Part 2.A: GSE113957 Analysis                         #
#                                                                              #
################################################################################

################################################################################
#                            2.A.i: DEseq2 Analysis                              #
################################################################################

#creating the matrix
GSE113957_dseq <- DESeqDataSetFromMatrix(countData = GSE113957_rna_seq,
                                        colData = GSE113957_meta, 
                                        design = ~ Group)
GSE113957_dseq <- GSE113957_dseq[rowSums(counts(GSE113957_dseq) >= 10) >= 3, ]
nrow(GSE113957_dseq)
GSE113957_dseq <- estimateSizeFactors(GSE113957_dseq)
print(sizeFactors(GSE113957_dseq))

GSE113957_dseq_data <- DESeq(GSE113957_dseq)
head(GSE113957_dseq_data)
resultsNames(GSE113957_dseq_data)

#normalising the data to cpm/fpm
GSE113957_dseq_cpm <- fpm(GSE113957_dseq_data, robust = TRUE)
GSE113957_dseq_cpm_cols <- colnames(GSE113957_dseq_cpm)

#plotting the data
GSE113957_dseq_cpm_df <- as.data.frame(GSE113957_dseq_cpm)
GSE113957_dseq_cpm_df <- rownames_to_column(GSE113957_dseq_cpm_df, "Gene")
GSE113957_dseq_cpm_df <- pivot_longer(GSE113957_dseq_cpm_df, cols = -Gene, names_to = "Sample", values_to = "Count")
GSE113957_dseq_cpm_df <- left_join(GSE113957_dseq_cpm_df, GSE113957_meta, by = "Sample")

#Boxplot displaying lcount per sample with cpm transformation
cpm_plot <- GSE113957_dseq_cpm_df%>%
  ggplot( aes(x = Count, y = Sample, fill = Group)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.4)  +
  scale_fill_manual(values = group_colours) +
  labs(title = "Boxplot of CPM for each sample", x = "Sample", y = "CPM") + 
  theme(legend.position="none",
        plot.title = element_text(size=11)) +
  theme_minimal() 
  
#plotting the log2 normalisation
GSE113957_dseq_log <- as.data.frame(log2(1+counts(GSE113957_dseq, normalized=TRUE)))
GSE113957_dseq_log <- rownames_to_column(GSE113957_dseq_log, "Gene")
GSE113957_dseq_log<- pivot_longer(GSE113957_dseq_log, cols = -Gene, names_to = "Sample", values_to = "Count")
GSE113957_dseq_log <- left_join(GSE113957_dseq_log, GSE113957_meta, by = "Sample")

##Boxplot displaying lcount per sample with log2 transformation
log_plot <- plot1 <- GSE113957_dseq_log %>%
  ggplot( aes(x = Count, y = Sample, fill = Group)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.4)  +
  scale_fill_manual(values = group_colours) +
  labs(title = "Boxplot of log2 for each sample", x = "Sample", y = "log2(Count)") + 
  theme(legend.position="none",
        plot.title = element_text(size=11)) +
  theme_minimal() 

(cpm_plot | log_plot)

GSE113957_effects <- DESeq(GSE113957_dseq)
resultsNames(GSE113957_effects)

# Creating contrasts for the analysis
GSE113957_hgps_cntrl_male_vs_fem <- results(GSE113957_effects, name = "Group_Control_male_vs_Control_female")
GSE113957_hgps_hgps_male_vs_fem <- results(GSE113957_effects, contrast = c("Group", "HGPS_male", "HGPS_female"))
GSE113957_male_hgps_vs_cntrl <- results(GSE113957_effects, contrast = c("Group", "HGPS_male", "Control_male"))
GSE113957_female_hgps_vs_cntrl <- results(GSE113957_effects, name = "Group_HGPS_female_vs_Control_female")

# Analysis between control males and control females #
#Filtering and data cleaning
GSE113957_hgps_cntrl_male_vs_fem <- as.data.frame(subset(GSE113957_hgps_cntrl_male_vs_fem, !is.na(padj)))
GSE113957_hgps_cntrl_male_vs_fem$Gene <- rownames(GSE113957_hgps_cntrl_male_vs_fem)

#checking for normality
ks.test(rank(GSE113957_hgps_cntrl_male_vs_fem$log2FoldChange, na.last = "keep"), "pnorm")
#using ranked z_scores because the log2fold change is nor normally distributed
GSE113957_hgps_cntrl_male_vs_fem$rank_zscore <- qnorm(rank(GSE113957_hgps_cntrl_male_vs_fem$log2FoldChange, na.last = "keep") / 
                                             (sum(!is.na(GSE113957_hgps_cntrl_male_vs_fem$log2FoldChange)) + 1))

GSE113957_hgps_cntrl_male_vs_fem$Gene <- rownames(GSE113957_hgps_cntrl_male_vs_fem)
GSE113957_hgps_cntrl_male_vs_fem_idx <- colData(GSE113957_dseq)$Group[c("Control_male", "Control_female")]


GSE113957_hgps_cntrl_male_vs_fem_sig <- subset(GSE113957_hgps_cntrl_male_vs_fem, padj < 0.05)
GSE113957_hgps_cntrl_male_vs_fem_sig <- GSE113957_hgps_cntrl_male_vs_fem_sig[order(GSE113957_hgps_cntrl_male_vs_fem_sig$padj, decreasing = FALSE), ]
nrow(GSE113957_hgps_cntrl_male_vs_fem_sig)
#write.csv(hgps_cntrl_male_vs_fem_sig, file = "hgps_cntrl_male_vs_fem_sig_dseq.csv", row.names = TRUE)

#isolating the genes upregulated in males
GSE113957_cntrl_male_up1_dseq <- subset(GSE113957_hgps_cntrl_male_vs_fem, padj < 0.05 & log2FoldChange >= 0.58)
nrow(GSE113957_cntrl_male_up1_dseq)

#isolating the genes upregulated in females
GSE113957_cntrl_fem_up1_dseq <- subset(GSE113957_hgps_cntrl_male_vs_fem, padj < 0.05 & log2FoldChange < -0.58)
nrow(GSE113957_cntrl_fem_up1_dseq)


#Following code is adapted from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
# GSE113957_hgps_cntrl_male_vs_fem$diffexpressed <- "Not Significant"
# GSE113957_hgps_cntrl_male_vs_fem$diffexpressed[GSE113957_hgps_cntrl_male_vs_fem$log2FoldChange >= 0.58 & GSE113957_hgps_cntrl_male_vs_fem$padj < 0.05] <- "Upregulated in Control Males"
# GSE113957_hgps_cntrl_male_vs_fem$diffexpressed[GSE113957_hgps_cntrl_male_vs_fem$log2FoldChange < -0.58 & GSE113957_hgps_cntrl_male_vs_fem$padj < 0.05] <- "Upregulated in Control Females"
# 
# #identifying the top 10 upregulated genes in each category
# top_cntrl_fem <- GSE113957_hgps_cntrl_male_vs_fem[GSE113957_hgps_cntrl_male_vs_fem$diffexpressed == "Upregulated in Control Males", ]
# top_cntrl_fem <- top_cntrl_fem[order(top_cntrl_fem$log2FoldChange, decreasing = TRUE), ]
# top_cntrl_fem <- top_cntrl_fem[1:20, ]
# 
# top_cntrl_male <- GSE113957_hgps_cntrl_male_vs_fem[GSE113957_hgps_cntrl_male_vs_fem$diffexpressed == "Upregulated in Control Males", ]
# top_cntrl_male <- top_cntrl_male[order(top_cntrl_male$log2FoldChange), ]
# top_cntrl_male <- top_cntrl_male[1:20, ]
# 
# top_cntrl <- bind_rows(top_cntrl_fem, top_cntrl_male)
# #label only the top genes in females
# GSE113957_hgps_cntrl_male_vs_fem$delabel <- NA
# GSE113957_hgps_cntrl_male_vs_fem$delabel[rownames(GSE113957_hgps_cntrl_male_vs_fem) %in% rownames(GSE113957_hgps_cntrl_male_vs_fem)] <- rownames(GSE113957_hgps_cntrl_male_vs_fem)[rownames(GSE113957_hgps_cntrl_male_vs_fem) %in% rownames(top_cntrl)]


#making a boxplot for degs in the control comparison 

#getting numbers of degs
degs_cntrl_males1 <- nrow(GSE113957_cntrl_male_up1_dseq)
degs_cntrl_fem1 <- nrow(GSE113957_cntrl_fem_up1_dseq)
degs_cntrl_tot <- degs_cntrl_males1 - degs_cntrl_fem1
#making a dataframe of the degs
GSE113957_cntrl_deg_counts <- data.frame(Category = c("Total", "Control Male-Biased", "Control Female-Biased"),
                                   Count = c(degs_cntrl_tot, degs_cntrl_males1, degs_cntrl_fem1))
GSE113957_cntrl_deg_counts$Category <- factor(GSE113957_cntrl_deg_counts$Category, levels = c("Total", "Control Male-Biased", "Control Female-Biased", "Shared"))


ggplot(GSE113957_cntrl_deg_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1800)) +
  scale_y_continuous(breaks = seq(0, 1800, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between Control Males and Control Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "Control Male-Biased" = "#4A90E2", "Control Female-Biased" = "#FF8C42", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

# Analysis between hgps males and hgps females #

#Filtering and data cleaning
GSE113957_hgps_male_vs_fem <- as.data.frame(subset(GSE113957_hgps_hgps_male_vs_fem, !is.na(padj)))
GSE113957_hgps_male_vs_fem$Gene <- rownames(GSE113957_hgps_male_vs_fem)

#checking for normality
jarque.bera.test(GSE113957_hgps_male_vs_fem$log2FoldChange)
#using ranked z_scores because the log2fold change is nor normally distributed
GSE113957_hgps_male_vs_fem$rank_zscore <- qnorm(rank(GSE113957_hgps_male_vs_fem$log2FoldChange, na.last = "keep") / 
                                          (sum(!is.na(GSE113957_hgps_male_vs_fem$log2FoldChange)) + 1))

GSE113957_hgps_male_vs_fem <- subset(GSE113957_hgps_male_vs_fem, !is.na(padj))
GSE113957_hgps_male_vs_fem_idx <- colData(GSE113957_dseq)$Group[c("HGPS_male", "HGPS_female")]
GSE113957_hgps_male_vs_fem_idx <- colData(GSE113957_dseq)$Group %in% c("HGPS_male", "HGPS_female")
GSE113957_hgps_male_vs_fem_sig <- subset(GSE113957_hgps_male_vs_fem, padj < 0.05)
GSE113957_hgps_male_vs_fem_sig <- GSE113957_hgps_male_vs_fem_sig[order(GSE113957_hgps_male_vs_fem_sig$padj, decreasing = FALSE), ]
nrow(GSE113957_hgps_male_vs_fem_sig)
#write.csv(hgps_hgps_male_vs_fem_sig, file = "hgps_hgps_male_vs_fem_sig_dseq.csv", row.names = TRUE)

#isolating the genes upregulated in males
GSE113957_hgps_male_up1_dseq <- subset(GSE113957_hgps_male_vs_fem, padj < 0.05 & log2FoldChange >= 0.58)
nrow(GSE113957_hgps_male_up1_dseq)

#isolating the genes upregulated in hgps females
GSE113957_hgps_fem_up1_dseq <- subset(GSE113957_hgps_male_vs_fem, padj < 0.05 & log2FoldChange < -0.58)
nrow(GSE113957_hgps_fem_up1_dseq)


#getting numbers of degs
degs_hgps_male1 <- nrow(GSE113957_hgps_male_up1_dseq)
degs_hgps_fem1 <- nrow(GSE113957_hgps_fem_up1_dseq )
degs_hgps_tot <- degs_hgps_male1 + degs_hgps_fem1

#making a dataframe of the degs
GSE113957_hgps_deg_counts <- data.frame(Category = c("Total", "HGPS Male-Biased", "HGPS Female-Biased"),
                                       Count = c(degs_hgps_tot, degs_hgps_male1, degs_hgps_fem1))
GSE113957_hgps_deg_counts$Category <- factor(GSE113957_hgps_deg_counts$Category, levels = c("Total", "HGPS Male-Biased", "HGPS Female-Biased"))



ggplot(GSE113957_hgps_deg_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1800)) +
  scale_y_continuous(breaks = seq(0, 1800, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between HGPS Males and HGPS Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "HGPS Female-Biased" = "#B35E00", "HGPS Male-Biased" = "#1F497D", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))


# Analysis between hgps and control females #
#Filtering and data cleaning
GSE113957_female_hgps_vs_cntrl <- as.data.frame(subset(GSE113957_female_hgps_vs_cntrl, !is.na(padj)))
GSE113957_female_hgps_vs_cntrl$Gene <- rownames(GSE113957_female_hgps_vs_cntrl)

#checking for normality
ks.test(rank(GSE113957_female_hgps_vs_cntrl$log2FoldChange, na.last = "keep"), "pnorm")
#using ranked z_scores because the log2fold change is nor normally distributed
GSE113957_female_hgps_vs_cntrl$rank_zscore <- qnorm(rank(GSE113957_female_hgps_vs_cntrl$log2FoldChange, na.last = "keep") / 
                                         (sum(!is.na(GSE113957_female_hgps_vs_cntrl$log2FoldChange)) + 1))

GSE113957_female_hgps_vs_cntrl <- subset(GSE113957_female_hgps_vs_cntrl, !is.na(padj))
GSE113957_female_hgps_vs_cntrl_idx <- colData(GSE113957_dseq)$Group[c("HGPS_female", "Control_female")]
GSE113957_female_hgps_vs_cntrl_idx <- colData(GSE113957_dseq)$Group %in% c("HGPS_female", "Control_female")


#Filtering for significance
GSE113957_female_hgps_vs_cntrl_sig <- subset(GSE113957_female_hgps_vs_cntrl, padj < 0.05)
nrow(GSE113957_female_hgps_vs_cntrl_sig)
GSE113957_female_hgps_vs_cntrl_sig <- GSE113957_female_hgps_vs_cntrl_sig[order(GSE113957_female_hgps_vs_cntrl_sig$padj, decreasing = FALSE), ]
#write.csv(female_hgps_vs_cntrl_sig, file = "female_hgps_vs_cntrl_sig_dseq.csv", row.names = TRUE)
head(GSE113957_female_hgps_vs_cntrl_sig)

#isolating the genes upregulated in hgps females
GSE113957_hgps_fem_up2_dseq <- subset(GSE113957_female_hgps_vs_cntrl_sig, padj < 0.05 & log2FoldChange >= 0.58)
nrow(GSE113957_hgps_fem_up2_dseq)

#isolating the genes upregulated in females
GSE113957_cntrl_fem_up2_dseq <- subset(GSE113957_female_hgps_vs_cntrl_sig, padj < 0.05 &  log2FoldChange < -0.58)
nrow(GSE113957_cntrl_fem_up2_dseq)

#getting numbers of degs
degs_hgps_fem2 <- nrow(GSE113957_hgps_fem_up2_dseq)
degs_cntrl_fem2 <- nrow(GSE113957_cntrl_fem_up2_dseq)
degs_fem_tot <- degs_hgps_fem2 + degs_cntrl_fem2

#chi-squared test to see if the numbers of sig degs is significantly different
# fill in your numbers
n1 <- nrow(GSE113957_female_hgps_vs_cntrl_sig);  x1 <- degs_hgps_fem2     # males: 214 / 16 833 DEGs
n2 <- nrow(GSE113957_female_hgps_vs_cntrl_sig);  x2 <- degs_cntrl_fem2   # females: 147 / 16 833 DEGs

tab <- matrix(c(x1, n1 - x1,
                x2, n2 - x2),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("HGPS", "Control"),
                              c("Significant", "Not_significant")))

# or, without the correction (gives the classic z statistic)
prop.test(tab, correct = FALSE)

#making a dataframe of the degs
GSE113957_fem_deg_counts <- data.frame(Category = c("Total", "HGPS Female-Biased", "Control Female-Biased"),
                                         Count = c(degs_fem_tot, degs_hgps_fem2, degs_cntrl_fem2))
GSE113957_fem_deg_counts$Category <- factor(GSE113957_fem_deg_counts$Category, levels = c("Total", "HGPS Female-Biased", "Control Female-Biased"))



ggplot(GSE113957_fem_deg_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1800)) +
  scale_y_continuous(breaks = seq(0, 1800, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between HGPS Females and Control Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "HGPS Female-Biased" = "#B35E00", "Control Female-Biased" = "#FF8C42", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))


# Analysis between hgps and control males #
#Filtering and data cleaning
GSE113957_male_hgps_vs_cntrl <- as.data.frame(subset(GSE113957_male_hgps_vs_cntrl, !is.na(padj)))
GSE113957_male_hgps_vs_cntrl$Gene <- rownames(GSE113957_male_hgps_vs_cntrl)

#checking for normality
jarque.bera.test(GSE113957_male_hgps_vs_cntrl$log2FoldChange)
#using ranked z_scores because the log2fold change is nor normally distributed
GSE113957_male_hgps_vs_cntrl$rank_zscore <- qnorm(rank(GSE113957_male_hgps_vs_cntrl$log2FoldChange, na.last = "keep") / 
                                          (sum(!is.na(GSE113957_male_hgps_vs_cntrl$log2FoldChange)) + 1))

GSE113957_male_hgps_vs_cntrl_idx <- colData(GSE113957_dseq)$Group[c("HGPS_male", "Control_male")]
GSE113957_male_hgps_vs_cntrl_idx <- colData(GSE113957_dseq)$Group %in% c("HGPS_male", "Control_male")

#Filtering for significance
GSE113957_male_hgps_vs_cntrl_sig <- subset(GSE113957_male_hgps_vs_cntrl, padj < 0.05)
GSE113957_male_hgps_vs_cntrl_sig <- GSE113957_male_hgps_vs_cntrl_sig[order(GSE113957_male_hgps_vs_cntrl_sig$padj, decreasing = FALSE), ]
#write.csv(male_hgps_vs_cntrl_sig, file = "male_hgps_vs_cntrl_sig_dseq.csv", row.names = TRUE)
head(GSE113957_male_hgps_vs_cntrl_sig)     
                               
#isolating the genes upregulated in hgps males
GSE113957_hgps_male_up2_dseq <- subset(GSE113957_male_hgps_vs_cntrl_sig, padj < 0.05 & log2FoldChange >= 0.58)
nrow(GSE113957_hgps_male_up2_dseq)

#isolating the genes upregulated in control males
GSE113957_cntrl_male_up2_dseq <- subset(GSE113957_male_hgps_vs_cntrl_sig, padj < 0.05 & rank_zscore < -0.58)
nrow(GSE113957_cntrl_male_up2_dseq)

#getting numbers of degs
degs_hgps_male2 <- nrow(GSE113957_hgps_male_up2_dseq)
degs_cntrl_male2 <- nrow(GSE113957_cntrl_male_up2_dseq)
degs_male_tot <- degs_hgps_male2 + degs_cntrl_male2 

#making a dataframe of the degs
GSE113957_male_deg_counts <- data.frame(Category = c("Total", "HGPS Male-Biased", "Control Male-Biased"),
                                       Count = c(degs_male_tot, degs_hgps_male2, degs_cntrl_male2))
GSE113957_male_deg_counts$Category <- factor(GSE113957_male_deg_counts$Category, levels = c("Total", "HGPS Male-Biased", "Control Male-Biased"))


ggplot(GSE113957_male_deg_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1800)) +
  scale_y_continuous(breaks = seq(0, 1800, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between HGPS Males and Control Males",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "HGPS Male-Biased" = "#1F497D", "Control Male-Biased" = "#4A90E2", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))


GSE113957_hgps_fem_uniq1 <- GSE113957_hgps_fem_up1_dseq[!rownames(GSE113957_hgps_fem_up1_dseq) %in% rownames(GSE113957_cntrl_fem_up1_dseq), ]
GSE113957_hgps_fem_uniq1 <- GSE113957_hgps_fem_uniq1[order(GSE113957_hgps_fem_uniq1$padj, decreasing = FALSE), ]
nrow(GSE113957_hgps_fem_uniq1)
head(GSE113957_hgps_fem_uniq1, 10)
GSE113957_hgps_fem_uniq1_names <- rownames(GSE113957_hgps_fem_uniq1)[1:20]

#plotting the data
hgps_male_fem_dseq_samples <- colnames(GSE113957_dseq_cpm)[GSE113957_hgps_male_vs_fem_idx]
hgps_male_fem_dseq_annot <- data.frame(Group = GSE113957_meta$Group[match(hgps_male_fem_dseq_samples, rownames(GSE113957_meta))],
                                       row.names = hgps_male_fem_dseq_samples)
hgps_male_fem_dseq_annot_cols <- list(Group = c("HGPS_female" = "#B35E00", "HGPS_male" = "#1F497D"))

pheatmap(GSE113957_dseq_cpm[GSE113957_hgps_fem_uniq1_names, GSE113957_hgps_male_vs_fem_idx],
         scale="row",
         show_rownames=T,
         main="Top 20 Genes Upregulated in HGPS Females Relative to HGPS Males",
         fontfamily = "Arial", 
         fontsize = 20,
         annotation_col= hgps_male_fem_dseq_annot,
         annotation_colors = hgps_male_fem_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)

GSE113957_hgps_fem_uniq2 <- GSE113957_hgps_fem_up2_dseq[!rownames(GSE113957_hgps_fem_up2_dseq) %in% rownames(GSE113957_cntrl_fem_up1_dseq), ]
nrow(GSE113957_hgps_fem_uniq2)
GSE113957_hgps_fem_uniq2 <- GSE113957_hgps_fem_uniq2[!rownames(GSE113957_hgps_fem_uniq2) %in% rownames(GSE113957_hgps_male_up2_dseq), ]
nrow(GSE113957_hgps_fem_uniq2)
GSE113957_hgps_fem_uniq2 <- GSE113957_hgps_fem_uniq2[order(GSE113957_hgps_fem_uniq2$padj, decreasing = FALSE), ]
head(GSE113957_hgps_fem_uniq2, 10)

#plotting the data
fem_hgps_cntrl_dseq_samples <- colnames(GSE113957_dseq_cpm)[GSE113957_female_hgps_vs_cntrl_idx]
fem_hgps_cntrl_dseq_annot <- data.frame(Group = GSE113957_meta$Group[match(fem_hgps_cntrl_dseq_samples, rownames(GSE113957_meta))],
                                       row.names = fem_hgps_cntrl_dseq_samples)
fem_hgps_cntrl_dseq_annot_cols <- list(Group = c("HGPS_female" = "#B35E00", "Control_female" = "#FF8C42"))

GSE113957_hgps_fem_uniq2_names <- rownames(GSE113957_hgps_fem_uniq2)[1:20]
pheatmap(GSE113957_dseq_cpm[GSE113957_hgps_fem_uniq2_names, GSE113957_female_hgps_vs_cntrl_idx],
         scale="row",
         show_rownames=T,
         main="Top 20 Genes Upregulated in HGPS Females Relative to Control Females",
         fontfamily = "Arial", 
         fontsize = 20,
         #fontsize_col = 20,
         annotation_col= fem_hgps_cntrl_dseq_annot,
         annotation_colors = fem_hgps_cntrl_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)

#unique male genes
GSE113957_hgps_male_uniq1 <- GSE113957_hgps_male_up1_dseq[!rownames(GSE113957_hgps_male_up1_dseq) %in% rownames(GSE113957_cntrl_male_up1_dseq), ]
GSE113957_hgps_male_uniq1 <- GSE113957_hgps_male_uniq1[order(GSE113957_hgps_male_uniq1$padj, decreasing = FALSE), ]
nrow(GSE113957_hgps_male_uniq1)
head(GSE113957_hgps_male_uniq1, 10)

#plotting the data
GSE113957_hgps_male_uniq1_names <- rownames(GSE113957_hgps_male_uniq1)[1:20]
#plotting the data
hgps_male_fem_dseq_samples <- colnames(GSE113957_dseq_cpm)[GSE113957_hgps_male_vs_fem_idx]
hgps_male_fem_dseq_annot <- data.frame(Group = GSE113957_meta$Group[match(hgps_male_fem_dseq_samples, rownames(GSE113957_meta))],
                                       row.names = hgps_male_fem_dseq_samples)
hgps_male_fem_dseq_annot_cols <- list(Group = c("HGPS_female" = "#B35E00", "HGPS_male" = "#1F497D"))

pheatmap(GSE113957_dseq_cpm[GSE113957_hgps_male_uniq1_names, GSE113957_hgps_male_vs_fem_idx],
         scale="row",
         show_rownames=T,
         main="Top 20 Genes Upregulated in HGPS Females Relative to HGPS Males",
         fontfamily = "Arial", 
         fontsize = 20,
         annotation_col= hgps_male_fem_dseq_annot,
         annotation_colors = hgps_male_fem_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)



GSE113957_hgps_male_uniq2 <-  GSE113957_hgps_male_up2_dseq[!rownames(GSE113957_hgps_male_up2_dseq) %in% rownames(GSE113957_cntrl_male_up1_dseq), ]
nrow(GSE113957_hgps_male_uniq2)
GSE113957_hgps_male_uniq2 <- GSE113957_hgps_male_uniq2[!rownames(GSE113957_hgps_male_uniq2) %in% rownames(GSE113957_hgps_fem_up2_dseq), ]
nrow(GSE113957_hgps_male_uniq2)
GSE113957_hgps_male_uniq2 <- GSE113957_hgps_male_uniq2[order(GSE113957_hgps_male_uniq2$padj, decreasing = FALSE), ]
head(GSE113957_hgps_male_uniq2, 10)

#plotting the data
male_hgps_cntrl_dseq_samples <- colnames(GSE113957_dseq_cpm)[GSE113957_male_hgps_vs_cntrl_idx]
male_hgps_cntrl_dseq_annot <- data.frame(Group = GSE113957_meta$Group[match(male_hgps_cntrl_dseq_samples, rownames(GSE113957_meta))],
                                        row.names = male_hgps_cntrl_dseq_samples)
male_hgps_cntrl_dseq_annot_cols <- list(Group = c("HGPS_male" = "#1F497D", "Control_male" = "#4A90E2"))

GSE113957_hgps_male_uniq2_names <- rownames(GSE113957_hgps_male_uniq2)[1:20]
pheatmap(GSE113957_dseq_cpm[GSE113957_hgps_male_uniq2_names, GSE113957_male_hgps_vs_cntrl_idx],
         scale="row",
         show_rownames=T,
         main="Top 20 Genes Upregulated in HGPS Males Relative to Control Males",
         fontfamily = "Arial", 
         fontsize = 20,
         #fontsize_col = 20,
         annotation_col= male_hgps_cntrl_dseq_annot,
         annotation_colors = male_hgps_cntrl_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)

#checking that there's no overlap between the two data sets
length(intersect(GSE113957_hgps_male_uniq1$Gene, GSE113957_hgps_fem_uniq1$Gene))
length(intersect(GSE113957_hgps_male_uniq2$Gene, GSE113957_hgps_fem_uniq2$Gene))

################################################################################
#                             2.A.ii: edgeR Analysis                              #
################################################################################
# After running RUVSeq to get W matrix
GSE113957_counts_matrix <- GSE113957_rna_seq
GSE113957_group <- factor(GSE113957_meta$Group)

GSE113957_edgr <- DGEList(counts = GSE113957_counts_matrix, group = GSE113957_group)
nrow(GSE113957_edgr)

#removing lowly expressed genes
GSE113957_edgr_keep <- filterByExpr(GSE113957_edgr, group = GSE113957_group)
length(GSE113957_edgr_keep)
GSE113957_edgr_counts <- GSE113957_edgr[GSE113957_edgr_keep, , keep.lib.sizes = FALSE]
nrow(GSE113957_edgr_counts)

#normalization
GSE113957_edgr_counts<- calcNormFactors(GSE113957_edgr_counts)

#Visualising the data
GSE113957_edgr_cpm <- cpm(GSE113957_edgr_counts, log = FALSE, prior.count=1)

GSE113957_edgr_cpm_df <- as.data.frame(GSE113957_edgr_cpm)
GSE113957_edgr_cpm_df <- rownames_to_column(GSE113957_edgr_cpm_df, "Gene")
GSE113957_edgr_cpm_df <- pivot_longer(GSE113957_edgr_cpm_df, cols = -Gene, names_to = "Sample", values_to = "Count")
GSE113957_edgr_cpm_df <- left_join(GSE113957_edgr_cpm_df, GSE113957_meta, by = "Sample")

#Boxplot displaying lcount per sample with cpm transformation
cpm_plot <- GSE113957_edgr_cpm_df%>%
  ggplot( aes(x = Count, y = Sample, fill = Group)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.4)  +
  scale_fill_manual(values = group_colours) +
  labs(title = "Boxplot of CPM for each sample", x = "Sample", y = "CPM") + 
  theme(legend.position="none",
        plot.title = element_text(size=11)) +
  theme_minimal() 

#creating the desing matrix
GSE113957_edgr_investigation <- model.matrix(~0 + Group, data = GSE113957_meta)
colnames(GSE113957_edgr_investigation) <- levels(GSE113957_group)

#estimating dispersion
GSE113957_edgr_disp <- estimateDisp(GSE113957_edgr_counts, GSE113957_edgr_investigation)

#fitting the model
GSE113957_edgr_fit <- glmFit(GSE113957_edgr_disp, GSE113957_edgr_investigation)

# Contrasts for the analysis
GSE113957_edgr_contrasts <- makeContrasts(fem_hgps_vs_control_e = HGPS_female - Control_female,
                                     male_hgps_vs_control_e = HGPS_male - Control_male,
                                     hgps_fem_vs_male_e = HGPS_male - HGPS_female,
                                     control_fem_vs_male_e = Control_male - Control_female,
                                     levels = GSE113957_edgr_investigation)

#control male vs control female
cntrl_fem_vs_male_edgr <- glmLRT(GSE113957_edgr_fit, contrast = GSE113957_edgr_contrasts[, "control_fem_vs_male_e"])
cntrl_fem_vs_male_edgr_res <- topTags(cntrl_fem_vs_male_edgr, n = Inf)$table
cntrl_fem_vs_male_edgr_res <- as.data.frame(subset(cntrl_fem_vs_male_edgr_res, !is.na(FDR)))
cntrl_fem_vs_male_edgr_res$Gene <- rownames(cntrl_fem_vs_male_edgr_res)

#checking for normality
ks.test(rank(cntrl_fem_vs_male_edgr_res$logFC, na.last = "keep"), "pnorm")
jarque.bera.test(cntrl_fem_vs_male_edgr_res$logFC)
#using ranked z_scores because the log2fold change is nor normally distributed
cntrl_fem_vs_male_edgr_res$rank_zscore <- qnorm(rank(cntrl_fem_vs_male_edgr_res$logFC, na.last = "keep") / 
                                       (sum(!is.na(cntrl_fem_vs_male_edgr_res$logFC)) + 1))

#filtering for globally significant genes in this analysis
cntrl_fem_vs_male_edgr_res_sig <- subset(cntrl_fem_vs_male_edgr_res, FDR < 0.05)
nrow(cntrl_fem_vs_male_edgr_res_sig)

#filtering for genes upregulated in control females
cntrl_fem_up1_edgr <- subset(cntrl_fem_vs_male_edgr_res, FDR < 0.05 & logFC < -0.58)
nrow(cntrl_fem_up1_edgr)
#filtering for genes upregulated in control males 2
cntrl_male_up1_edgr <- subset(cntrl_fem_vs_male_edgr_res, FDR < 0.05 & logFC >= 0.58)
nrow(cntrl_male_up1_edgr)

degs_cntrl_male1_e <- nrow(cntrl_male_up1_edgr)
degs_cntrl_fem1_e <- nrow(cntrl_fem_up1_edgr)
degs_cntrl_tot_e <- degs_cntrl_male1_e + degs_cntrl_fem1_e

#making a dataframe of the degs
GSE113957_cntrl_deg_counts_e <- data.frame(Category = c("Total", "Control Male-Biased", "Control Female-Biased"),
                                        Count = c(degs_cntrl_tot_e, degs_cntrl_male1_e, degs_cntrl_fem1_e))
GSE113957_cntrl_deg_counts_e$Category <- factor(GSE113957_cntrl_deg_counts_e$Category, levels = c("Total", "Control Male-Biased", "Control Female-Biased"))


ggplot(GSE113957_cntrl_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1600)) +
  scale_y_continuous(breaks = seq(0, 1600, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between Control Males and Control Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "Control Female-Biased" = "#FF8C42", "Control Male-Biased" = "#4A90E2", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

# hgps male vs hgps female
hgps_fem_vs_male_edgr <- glmLRT(GSE113957_edgr_fit, contrast = GSE113957_edgr_contrasts[, "hgps_fem_vs_male_e"])
hgps_fem_vs_male_edgr_res <- topTags(hgps_fem_vs_male_edgr, n = Inf)$table
hgps_fem_vs_male_edgr_res <- as.data.frame(subset(hgps_fem_vs_male_edgr_res, !is.na(FDR)))
hgps_fem_vs_male_edgr_res$Gene <- rownames(hgps_fem_vs_male_edgr_res)

#checking for normality
ks.test(rank(hgps_fem_vs_male_edgr_res$logFC, na.last = "keep"), "pnorm")
jarque.bera.test(hgps_fem_vs_male_edgr_res$logFC)
#using ranked z_scores because the log2fold change is nor normally distributed
hgps_fem_vs_male_edgr_res$rank_zscore <- qnorm(rank(hgps_fem_vs_male_edgr_res$logFC, na.last = "keep") / 
                                                  (sum(!is.na(hgps_fem_vs_male_edgr_res$logFC)) + 1))


hgps_fem_vs_male_edgr_res_sig <- subset(hgps_fem_vs_male_edgr_res, FDR < 0.05)
nrow(hgps_fem_vs_male_edgr_res_sig)

#filtering for genes upregulated in hgps females
hgps_fem_up1_edgr <- subset(hgps_fem_vs_male_edgr_res, FDR < 0.05 & logFC < -0.58)
nrow(hgps_fem_up1_edgr)
#filtering for genes upregulated in control males 2
hgps_male_up1_edgr <- subset(hgps_fem_vs_male_edgr_res, FDR < 0.05 & logFC >= 0.58)
nrow(hgps_male_up1_edgr)

#counting the degs
degs_hgps_male1_e <- nrow(hgps_male_up1_edgr)
degs_hgps_fem1_e <- nrow(hgps_fem_up1_edgr)
degs_hgps_tot_e <- degs_hgps_male1_e + degs_hgps_fem1_e

#making a dataframe of the degs
GSE113957_hgps_deg_counts_e <- data.frame(Category = c("Total", "HGPS Male-Biased", "HGPS Female-Biased"),
                                           Count = c(degs_hgps_tot_e, degs_hgps_male1_e, degs_hgps_fem1_e))
GSE113957_hgps_deg_counts_e$Category <- factor(GSE113957_hgps_deg_counts_e$Category, levels = c("Total", "HGPS Male-Biased", "HGPS Female-Biased"))


ggplot(GSE113957_hgps_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1600)) +
  scale_y_continuous(breaks = seq(0, 1600, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between HGPS Males and HGPS Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "HGPS Female-Biased" = "#B35E00", "HGPS Male-Biased" = "#1F497D", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

# hgps female vs control female
fem_hgps_vs_cntrl_edgr <- glmLRT(GSE113957_edgr_fit, contrast = GSE113957_edgr_contrasts[, "fem_hgps_vs_control_e"])
fem_hgps_vs_cntrl_edgr_res <- topTags(fem_hgps_vs_cntrl_edgr, n = Inf)$table
fem_hgps_vs_cntrl_edgr_res <- as.data.frame(subset(fem_hgps_vs_cntrl_edgr_res, !is.na(FDR)))
fem_hgps_vs_cntrl_edgr_res$Gene <- rownames(fem_hgps_vs_cntrl_edgr_res)

#checking for normality
ks.test(rank(fem_hgps_vs_cntrl_edgr_res$logFC, na.last = "keep"), "pnorm")
jarque.bera.test(fem_hgps_vs_cntrl_edgr_res$logFC)
#using ranked z_scores because the log2fold change is nor normally distributed
fem_hgps_vs_cntrl_edgr_res$rank_zscore <- qnorm(rank(fem_hgps_vs_cntrl_edgr_res$logFC, na.last = "keep") / 
                                                 (sum(!is.na(fem_hgps_vs_cntrl_edgr_res$logFC)) + 1))


fem_hgps_vs_cntrl_edgr_res_sig <- subset(fem_hgps_vs_cntrl_edgr_res, FDR < 0.05)

#filtering for genes upregulated in hgps females
hgps_fem_up2_edgr <- subset(fem_hgps_vs_cntrl_edgr_res, FDR < 0.05 & logFC >= 0.58)
#filtering for genes upregulated in control males 2
cntrl_fem_up2_edgr <- subset(fem_hgps_vs_cntrl_edgr_res, FDR < 0.05 & logFC < -0.58)

#counting the degs
degs_hgps_fem2_e <- nrow(hgps_fem_up2_edgr)
degs_cntrl_fem2_e <- nrow(cntrl_fem_up2_edgr)
degs_fem_tot_e <- nrow(hgps_fem_up2_edgr) + nrow(cntrl_fem_up2_edgr)

#making a dataframe of the degs
GSE113957_fem_deg_counts_e <- data.frame(Category = c("Total", "HGPS Female-Biased", "Control Female-Biased"),
                                          Count = c(degs_fem_tot_e, degs_hgps_fem2_e, degs_cntrl_fem2_e))
GSE113957_fem_deg_counts_e$Category <- factor(GSE113957_fem_deg_counts_e$Category, levels = c("Total", "HGPS Female-Biased", "Control Female-Biased"))


ggplot(GSE113957_fem_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1600)) +
  scale_y_continuous(breaks = seq(0, 1600, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between HGPS Females and Control Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "HGPS Female-Biased" = "#B35E00", "Control Female-Biased" = "#FF8C42", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

# hgps male vs control male #
male_hgps_vs_cntrl_edgr <- glmLRT(GSE113957_edgr_fit, contrast = GSE113957_edgr_contrasts[, "male_hgps_vs_control_e"])
male_hgps_vs_cntrl_edgr_res <- topTags(male_hgps_vs_cntrl_edgr, n = Inf)$table
male_hgps_vs_cntrl_edgr_res <- as.data.frame(subset(male_hgps_vs_cntrl_edgr_res, !is.na(FDR)))
male_hgps_vs_cntrl_edgr_res$Gene <- rownames(male_hgps_vs_cntrl_edgr_res)

#checking for normality
ks.test(rank(male_hgps_vs_cntrl_edgr_res$logFC, na.last = "keep"), "pnorm")
jarque.bera.test(male_hgps_vs_cntrl_edgr_res$logFC)
#using ranked z_scores because the log2fold change is nor normally distributed
male_hgps_vs_cntrl_edgr_res$rank_zscore <- qnorm(rank(male_hgps_vs_cntrl_edgr_res$logFC, na.last = "keep") / 
                                                  (sum(!is.na(male_hgps_vs_cntrl_edgr_res$logFC)) + 1))

male_hgps_vs_cntrl_edgr_res_sig <- subset(male_hgps_vs_cntrl_edgr_res, FDR < 0.05)

#filtering for genes upregulated in hgps males
hgps_male_up2_edgr <- subset(male_hgps_vs_cntrl_edgr_res, FDR < 0.05 & logFC >= 0.58)
#filtering for genes upregulated in control males 2
cntrl_male_up2_edgr <- subset(male_hgps_vs_cntrl_edgr_res, FDR < 0.05 & logFC < 0.58)

#counting the degs
degs_hgps_male2_e <- nrow(hgps_male_up2_edgr)
degs_cntrl_male2_e <- nrow(cntrl_male_up2_edgr)
degs_male_tot_e <- degs_hgps_male2_e + degs_cntrl_male2_e

#making a dataframe of the degs
GSE113957_male_deg_counts_e <- data.frame(Category = c("Total", "HGPS Male-Biased", "Control Male-Biased"),
                                         Count = c(degs_male_tot_e, degs_hgps_male2_e, degs_cntrl_male2_e))
GSE113957_male_deg_counts_e$Category <- factor(GSE113957_male_deg_counts_e$Category, levels = c("Total", "HGPS Male-Biased", "Control Male-Biased"))


ggplot(GSE113957_male_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 1600)) +
  scale_y_continuous(breaks = seq(0, 1600, by = 100)) +
  labs(title = "Number of Differentially Expressed Genes Between HGPS Males and Control Males",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "HGPS Male-Biased" = "#1F497D", "Control Male-Biased" = "#4A90E2", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))


#combining the data frames
#version 1: only removing the control fem genes from the hgps fem1 analysis
hgps_fem_specific_edgr1 <- hgps_fem_up1_edgr[!rownames(hgps_fem_up1_edgr) %in% rownames(cntrl_fem_up1_edgr), ]
hgps_fem_specific_edgr1 <- hgps_fem_specific_edgr1[order(hgps_fem_specific_edgr1$FDR, decreasing = FALSE), ]
nrow(hgps_fem_specific_edgr1)
head(hgps_fem_specific_edgr1)

#version 2: combining the two hgps fem analyses and removing duplicates
hgps_fem_specific_edgr2 <-  hgps_fem_up2_edgr[!rownames(hgps_fem_up2_edgr) %in% rownames(cntrl_fem_up1_edgr), ]
nrow(hgps_fem_specific_edgr2)
hgps_fem_specific_edgr2 <- hgps_fem_specific_edgr2[!rownames(hgps_fem_specific_edgr2) %in% rownames(hgps_male_up2_edgr), ]
nrow(hgps_fem_specific_edgr2)
hgps_fem_specific_edgr2 <- hgps_fem_specific_edgr2[order(hgps_fem_specific_edgr2$FDR, decreasing = FALSE), ]
head(hgps_fem_specific_edgr2, 10)


#version 1: only removing the control fem genes from the hgps fem1 analysis
hgps_male_specific_edgr1 <- hgps_male_up1_edgr[!rownames(hgps_male_up1_edgr) %in% rownames(cntrl_male_up1_edgr), ]
hgps_male_specific_edgr1 <- hgps_male_specific_edgr1[order(hgps_male_specific_edgr1$FDR, decreasing = FALSE), ]
nrow(hgps_male_specific_edgr1)
head(hgps_male_specific_edgr1, 10)

#version 2: combining the two hgps fem analyses and removing duplicates
hgps_male_specific_edgr2 <- hgps_male_up2_edgr[!rownames(hgps_male_up2_edgr) %in% rownames(cntrl_male_up1_edgr), ]
nrow(hgps_male_specific_edgr2)
hgps_male_specific_edgr2 <- hgps_male_specific_edgr2[!rownames(hgps_male_specific_edgr2) %in% rownames(hgps_fem_up2_edgr), ]
nrow(hgps_male_specific_edgr2)
hgps_male_specific_edgr2 <- hgps_male_specific_edgr2[order(hgps_male_specific_edgr2$FDR, decreasing = FALSE), ]
head(hgps_male_specific_edgr2, 10)

################################################################################
#                 2.C: Combining the deseq and edgr analysis                   #
################################################################################
#Getting the intersect of the uniquely expressed genes in HGPS females in the deseq2 and edgeR analyses
hgps_fem_specific_common1 <- inner_join(GSE113957_hgps_fem_uniq1, hgps_fem_specific_edgr1, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(hgps_fem_specific_common1)

hgps_fem_specific_common2 <- inner_join(GSE113957_hgps_fem_uniq2, hgps_fem_specific_edgr2, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(hgps_fem_specific_common2)

#reorganising the dataframe
hgps_fem_specific_common1 <- hgps_fem_specific_common1[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                         "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]
hgps_fem_specific_common2 <- hgps_fem_specific_common2[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                           "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]


hgps_fem_specific_common1 <- hgps_fem_specific_common1[!duplicated(hgps_fem_specific_common1$Gene), ]
rownames(hgps_fem_specific_common1) <- hgps_fem_specific_common1$Gene
nrow(hgps_fem_specific_common1)

hgps_fem_specific_common2 <- hgps_fem_specific_common2[!duplicated(hgps_fem_specific_common2$Gene), ]
rownames(hgps_fem_specific_common2) <- hgps_fem_specific_common2$Gene
nrow(hgps_fem_specific_common2)

hgps_fem_specific_common3 <- hgps_fem_specific_common1[ hgps_fem_specific_common1$Gene %in% hgps_fem_specific_common2$Gene, ]
nrow(hgps_fem_specific_common3)

#organise by log2foldchange
nrow(hgps_fem_specific_common1)
hgps_fem_specific_common1 <- hgps_fem_specific_common1[order(hgps_fem_specific_common1$log2FoldChange, decreasing = TRUE), ]
hgps_fem_specific_common1_names <- rownames(hgps_fem_specific_common1)
length(hgps_fem_specific_common1_names)
head(hgps_fem_specific_common1, 10)

hgps_fem_specific_common1_annot <- data.frame(Group = GSE113957_meta$Group[match(hgps_male_fem_dseq_samples, rownames(GSE113957_meta))],
                                              row.names = hgps_male_fem_dseq_samples)
hgps_male_fem_dseq_annot_cols <- list(Group = c("HGPS_female" = "#B35E00", "HGPS_male" = "#1F497D"))
pheatmap(GSE113957_dseq_cpm[hgps_fem_specific_common1_names, GSE113957_hgps_male_vs_fem_idx],
         scale="row",
         show_rownames=T,
         main="Top 20 Genes Upregulated in HGPS Females Relative to HGPS Males",
         fontfamily = "Arial", 
         fontsize = 20,
         annotation_col= hgps_male_fem_dseq_annot,
         annotation_colors = hgps_male_fem_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)

nrow(hgps_fem_specific_common2)
hgps_fem_specific_common2 <- hgps_fem_specific_common2[order(hgps_fem_specific_common2$log2FoldChange, decreasing = TRUE), ]
head(hgps_fem_specific_common2, 10)
hgps_fem_specific_common2_names <- rownames(hgps_fem_specific_common2)

hgps_fem_specific_common2_annot <- data.frame(Group = GSE113957_meta$Group[match(fem_hgps_cntrl_dseq_samples, rownames(GSE113957_meta))],
                                              row.names = fem_hgps_cntrl_dseq_samples)
pheatmap(GSE113957_dseq_cpm[hgps_fem_specific_common2_names, GSE113957_female_hgps_vs_cntrl_idx],
         scale="row",
         show_rownames=F,
         main="Top 20 Genes Upregulated in HGPS Females Relative to Control Females",
         fontfamily = "Arial", 
         fontsize = 20,
         annotation_col= fem_hgps_cntrl_dseq_annot,
         annotation_colors = fem_hgps_cntrl_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)


hgps_fem_specific_common3 <- hgps_fem_specific_common1[ hgps_fem_specific_common1$Gene %in% hgps_fem_specific_common2$Gene, ]
nrow(hgps_fem_specific_common3)
head(hgps_fem_specific_common3)


#Getting the intersect of the uniquely expressed genes in HGPS females in the deseq2 and edgeR analyses
hgps_male_specific_common1 <- inner_join(GSE113957_hgps_male_uniq1, hgps_male_specific_edgr1, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(hgps_male_specific_common1)

hgps_male_specific_common2 <- inner_join(GSE113957_hgps_male_uniq2, hgps_male_specific_edgr2, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(hgps_male_specific_common2)

#reorganising the dataframe
hgps_male_specific_common1 <- hgps_male_specific_common1[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                         "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]
hgps_male_specific_common2 <- hgps_male_specific_common2[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                             "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]

hgps_male_specific_common1 <- hgps_male_specific_common1[!duplicated(hgps_male_specific_common1$Gene), ]
rownames(hgps_male_specific_common1) <- hgps_male_specific_common1$Gene

hgps_male_specific_common2 <- hgps_male_specific_common2[!duplicated(hgps_male_specific_common2$Gene), ]
rownames(hgps_male_specific_common2) <- hgps_male_specific_common2$Gene

#organise by padj
hgps_male_specific_common1 <- hgps_male_specific_common1[order(hgps_male_specific_common1$log2FoldChange, decreasing = TRUE), ]
head(hgps_male_specific_common1, 10)

hgps_male_specific_common1_names <- rownames(hgps_male_specific_common1)
length(hgps_male_specific_common1_names)
head(hgps_male_specific_common1, 10)

hgps_male_specific_common1_annot <- data.frame(Group = GSE113957_meta$Group[match(hgps_male_fem_dseq_samples, rownames(GSE113957_meta))],
                                              row.names = hgps_male_fem_dseq_samples)
hgps_male_fem_dseq_annot_cols <- list(Group = c("HGPS_female" = "#B35E00", "HGPS_male" = "#1F497D"))
pheatmap(GSE113957_dseq_cpm[hgps_male_specific_common1_names, GSE113957_hgps_male_vs_fem_idx],
         scale="row",
         show_rownames=T,
         main="Top 20 Genes Upregulated in HGPS Males Relative to HGPS Females",
         fontfamily = "Arial", 
         fontsize = 20,
         annotation_col= hgps_male_fem_dseq_annot,
         annotation_colors = hgps_male_fem_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)


hgps_male_specific_common2 <- hgps_male_specific_common2[order(hgps_male_specific_common2$log2FoldChange, decreasing = TRUE), ]
head(hgps_male_specific_common2, 10)
hgps_male_specific_common2_names <- rownames(hgps_male_specific_common2)

hgps_male_specific_common2_annot <- data.frame(Group = GSE113957_meta$Group[match(male_hgps_cntrl_dseq_samples, rownames(GSE113957_meta))],
                                              row.names = male_hgps_cntrl_dseq_samples)
pheatmap(GSE113957_dseq_cpm[hgps_male_specific_common2_names, GSE113957_male_hgps_vs_cntrl_idx],
         scale="row",
         show_rownames=F,
         main="Top 20 Genes Upregulated in HGPS Males Relative to Control Males",
         fontfamily = "Arial", 
         fontsize = 20,
         annotation_col= male_hgps_cntrl_dseq_annot,
         annotation_colors = male_hgps_cntrl_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)


hgps_male_specific_common3 <- hgps_male_specific_common1[ hgps_male_specific_common1$Gene %in% hgps_male_specific_common2$Gene, ]
nrow(hgps_male_specific_common3)
head(hgps_male_specific_common3)


################################################################################
# Communal dseq2 and edgr analysis to make the volacno plots
#dseq data set
GSE113957_hgps_male_vs_fem
GSE113957_hgps_cntrl_male_vs_fem

#edger data set
hgps_fem_vs_male_edgr_res
cntrl_fem_vs_male_edgr_res

cntrl_common <- inner_join(GSE113957_hgps_cntrl_male_vs_fem, cntrl_fem_vs_male_edgr_res, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
rownames(cntrl_common) <- cntrl_common$Gene

hgps_common <- inner_join(GSE113957_hgps_male_vs_fem, hgps_fem_vs_male_edgr_res, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(hgps_common)
rownames(hgps_common) <- hgps_common$Gene
#hgps_common <- hgps_common[!rownames(hgps_common) %in% rownames(cntrl_common), ]
hgps_common <- hgps_common[!rownames(hgps_common) %in% rownames(cntrl_fem_up1_edgr), ]
nrow(hgps_common)
hgps_common <- hgps_common[!rownames(hgps_common) %in% rownames(GSE113957_cntrl_fem_up1_dseq), ]
nrow(hgps_common)
hgps_common <- hgps_common[!rownames(hgps_common) %in% rownames(cntrl_male_up1_edgr), ]
nrow(hgps_common)
hgps_common <- hgps_common[!rownames(hgps_common) %in% rownames(GSE113957_cntrl_male_up1_dseq), ]
nrow(hgps_common)
hgps_common <- hgps_common[order(hgps_common$log2FoldChange, decreasing = TRUE), ]
nrow(hgps_common)
hgps_common <- hgps_common[order(hgps_common$log2FoldChange, decreasing = TRUE), ]

#chi-squared test
hgps_male_deg1 <- nrow(subset(hgps_common, padj < 0.05 & FDR < 0.05 & log2FoldChange >= 0.58 & logFC >= 0.58))
hgps_fem_deg1 <- nrow(subset(hgps_common, padj < 0.05 & FDR < 0.05 & log2FoldChange < -0.58 & logFC < -0.58))
hgps_deg <- nrow(subset(hgps_common, padj < 0.05 & FDR < 0.05))

# creating the matrix
n1 <- hgps_deg;  x1 <- hgps_male_deg1
n2 <- hgps_deg;  x2 <- hgps_fem_deg1    

tab <- matrix(c(x1, n1 - x1,
                x2, n2 - x2),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("Male", "Female"),
                              c("Significant", "Not_significant")))

#testing for significance
prop.test(tab, correct = TRUE)


hgps_common_df <- hgps_common %>%
  # calculate -log10 padj from DESeq2
  mutate(
    negLog10Padj = -log10(padj),
    # call up/down by DESeq2 direction
    status = case_when(
      padj < 0.05 & FDR < 0.05 & log2FoldChange >=  0.58 ~ "Up in HGPS Males",
      padj < 0.05 & FDR < 0.05 & log2FoldChange <  -0.58 ~ "Up in HGPS Females",
      TRUE                        ~ "Neutral"
    ),
    status = factor(status, levels = c("Up in HGPS Males", "Up in HGPS Females", "Neutral"))
  )
rownames(hgps_common_df) <- hgps_common_df$Gene
hgps_common_df <- hgps_common_df[order(hgps_common_df$log2FoldChange, decreasing = FALSE), ]

#removing some outliers for aesthetic plot purposes
hgps_common_df <- subset(hgps_common_df, log2FoldChange > -9)


# 2) choose a few top genes to label (by absolute DESeq2 LFC)
# how many top genes you want
n_top <- 20

# with slice_max (recommended)
hgps_common_label <- hgps_common_df %>%
  slice_max(order_by = abs(log2FoldChange), n = n_top)

# pick top-by-direction
lab_fem <- hgps_common_df %>%
  filter(status == "Up in HGPS Females") %>%
  slice_max(order_by = log2FoldChange, n = 9, with_ties = TRUE)

lab_male <- hgps_common_df %>%
  filter(status == "Up in HGPS Males") %>%
  slice_min(order_by = log2FoldChange, n = 20, with_ties = FALSE)  # most negative

set.seed(123)  # reproducible label placement
ggplot(hgps_common_df,
       aes(x = log2FoldChange, y = negLog10Padj, color = status)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # label each group separately
  geom_text_repel(data = lab_fem, aes(label = Gene),
                  size = 2.6, max.overlaps = Inf, box.padding = 0.3, segment.size = 0.2,
                  show.legend = FALSE) +
  geom_text_repel(data = lab_male, aes(label = Gene),
                  size = 2.6, max.overlaps = Inf, box.padding = 0.3, segment.size = 0.2,
                  show.legend = FALSE) +
  
  scale_color_manual(values = c(
    "Up in HGPS Females" = "#B35E00",
    "Up in HGPS Males"   = "#1F497D",
    "Neutral"            = "grey70"
  )) +
  labs(
    title = "Volcano Plot: Common DEGs in HGPS Females vs HGPS Males",
    x = expression(log[2]~Fold~Change~"(DESeq2)"),
    y = expression(-log[10]~Adjusted~P~Value),
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


female_common <- inner_join(GSE113957_female_hgps_vs_cntrl, fem_hgps_vs_cntrl_edgr_res, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
rownames(female_common) <- female_common$Gene
nrow(female_common)
female_common <- female_common[!rownames(female_common) %in% rownames(cntrl_fem_up1_edgr), ]
nrow(female_common)
female_common <- female_common[!rownames(female_common) %in% rownames(hgps_male_up2_edgr), ]
nrow(female_common)
female_common <- female_common[!rownames(female_common) %in% rownames(GSE113957_cntrl_fem_up1_dseq), ]
nrow(female_common)
female_common <- female_common[!rownames(female_common) %in% rownames(GSE113957_hgps_male_up2_dseq), ]
nrow(female_common)
female_common <- female_common[order(female_common$log2FoldChange, decreasing = TRUE), ]


#chi-squared test
hgps_fem_deg2 <- nrow(subset(female_common, padj < 0.05 & FDR < 0.05 & log2FoldChange >= 0.58  & logFC >= 0.58))
cntrl_fem_deg1 <- nrow(subset(female_common, padj < 0.05 & FDR < 0.05 & log2FoldChange < -0.58  & logFC < -0.58))
female_deg <- nrow(subset(female_common, padj < 0.05 & FDR < 0.05))

# creating the matrix
n1 <- female_deg;  x1 <- hgps_fem_deg2
n2 <- female_deg;  x2 <- cntrl_fem_deg1    

tab <- matrix(c(x1, n1 - x1,
                x2, n2 - x2),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("Male", "Female"),
                              c("Significant", "Not_significant")))

#testing for significance
prop.test(tab, correct = TRUE)


female_common_df <- female_common %>%
  # calculate -log10 padj from DESeq2
  mutate(
    negLog10Padj = -log10(padj),
    # call up/down by DESeq2 direction
    status = case_when(
      padj < 0.05 & FDR < 0.05 & log2FoldChange >=  0.58 ~ "Up in HGPS Females",
      padj < 0.05 & FDR < 0.05 & log2FoldChange <  -0.58 ~ "Up in Control Females",
      TRUE                        ~ "Neutral"
    ),
    status = factor(status, levels = c("Up in HGPS Females", "Up in Control Females", "Neutral"))
  )
rownames(female_common_df) <- female_common_df$Gene
female_common_df <- female_common_df[order(female_common_df$log2FoldChange, decreasing = TRUE), ]


# with slice_max (recommended)
female_common_label <- female_common_df %>%
  slice_max(order_by = abs(log2FoldChange), n = n_top)

# pick top-by-direction
lab_fem_h <- female_common_df %>%
  filter(status == "Up in HGPS Females") %>%
  slice_max(order_by = log2FoldChange, n = 20, with_ties = TRUE)

lab_fem_c <- female_common_df %>%
  filter(status == "Up in Control Females") %>%
  slice_min(order_by = log2FoldChange, n = 20, with_ties = FALSE)  # most negative

set.seed(123)  # reproducible label placement
ggplot(female_common_df,
       aes(x = log2FoldChange, y = negLog10Padj, color = status)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # label each group separately
  geom_text_repel(data = lab_fem_h, aes(label = Gene),
                  size = 2.6, max.overlaps = Inf, box.padding = 0.3, segment.size = 0.2,
                  show.legend = FALSE) +
  geom_text_repel(data = lab_fem_c, aes(label = Gene),
                  size = 2.6, max.overlaps = Inf, box.padding = 0.3, segment.size = 0.2,
                  show.legend = FALSE) +
  
  scale_color_manual(values = c(
    "Up in HGPS Females" = "#B35E00",
    "Up in Control Females"   = "#FF8C42",
    "Neutral"            = "grey70"
  )) +
  labs(
    title = "Volcano Plot: Common DEGs in HGPS Females vs Control Females",
    x = expression(log[2]~Fold~Change~"(DESeq2)"),
    y = expression(-log[10]~Adjusted~P~Value),
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


male_common <- inner_join(GSE113957_male_hgps_vs_cntrl, male_hgps_vs_cntrl_edgr_res, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
rownames(male_common) <- male_common$Gene
nrow(male_common)
male_common <- male_common[!rownames(male_common) %in% rownames(cntrl_male_up1_edgr), ]
nrow(male_common)
male_common <- male_common[!rownames(male_common) %in% rownames(hgps_fem_up2_edgr), ]
nrow(male_common)
male_common <- male_common[!rownames(male_common) %in% rownames(GSE113957_cntrl_male_up1_dseq), ]
nrow(male_common)
male_common <- male_common[!rownames(male_common) %in% rownames(GSE113957_hgps_fem_up2_dseq), ]
nrow(male_common)
male_common <- male_common[order(male_common$log2FoldChange, decreasing = TRUE), ]

hgps_male_deg2 <- nrow(subset(male_common, padj < 0.05 & FDR < 0.05 & log2FoldChange >= 0.58 & logFC >= 0.58))
cntrl_male_deg1 <- nrow(subset(male_common, padj < 0.05 & FDR < 0.05 & log2FoldChange < -0.58  & logFC < -0.58))
male_deg <- nrow(subset(male_common, padj < 0.05 & FDR < 0.05))

# creating the matrix
n1 <- male_deg;  x1 <- hgps_male_deg2
n2 <- male_deg;  x2 <- cntrl_male_deg1    

tab <- matrix(c(x1, n1 - x1,
                x2, n2 - x2),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("HGPS Male", "Control Male"),
                              c("Significant", "Not_significant")))

#testing for significance
prop.test(tab, correct = TRUE)


male_common_df <- male_common %>%
  # calculate -log10 padj from DESeq2
  mutate(
    negLog10Padj = -log10(padj),
    # call up/down by DESeq2 direction
    status = case_when(
      padj < 0.05 & FDR < 0.05 & log2FoldChange >=  0.58 ~ "Up in HGPS Males",
      padj < 0.05 & FDR < 0.05 & log2FoldChange <  -0.58 ~ "Up in Control Males",
      TRUE                        ~ "Neutral"
    ),
    status = factor(status, levels = c("Up in HGPS Males", "Up in Control Males", "Neutral"))
  )
rownames(male_common_df) <- male_common_df$Gene
male_common_df <- male_common_df[order(male_common_df$log2FoldChange, decreasing = TRUE), ]

#removing some outliers for aesthetic plot purposes
male_common_df <- subset(male_common_df, log2FoldChange > -9)

# with slice_max (recommended)
male_common_label <- male_common_df %>%
  slice_max(order_by = abs(log2FoldChange), n = n_top)

# pick top-by-direction
lab_male_h <- male_common_df %>%
  filter(status == "Up in HGPS Males") %>%
  slice_max(order_by = log2FoldChange, n = 20, with_ties = FALSE)

lab_male_c <- male_common_df %>%
  filter(status == "Up in Control Males") %>%
  slice_min(order_by = log2FoldChange, n = 20, with_ties = FALSE)  # most negative

set.seed(123)  # reproducible label placement
ggplot(male_common_df,
       aes(x = log2FoldChange, y = negLog10Padj, color = status)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # label each group separately
  geom_text_repel(data = lab_fem_h, aes(label = Gene),
                  size = 2.6, max.overlaps = Inf, box.padding = 0.3, segment.size = 0.2,
                  show.legend = FALSE) +
  geom_text_repel(data = lab_fem_c, aes(label = Gene),
                  size = 2.6, max.overlaps = Inf, box.padding = 0.3, segment.size = 0.2,
                  show.legend = FALSE) +
  
  scale_color_manual(values = c(
    "Up in HGPS Males" = "#1F497D",
    "Up in Control Males"   = "#4A90E2",
    "Neutral"            = "grey70"
  )) +
  labs(
    title = "Volcano Plot: Common DEGs in HGPS Males vs Control Males",
    x = expression(log[2]~Fold~Change~"(DESeq2)"),
    y = expression(-log[10]~Adjusted~P~Value),
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


cor(hgps_common$padj, hgps_common$FDR, method = "pearson")
cor.test(hgps_common$padj, hgps_common$FDR)
ggplot(hgps_common) +
  aes(x = padj, y = FDR) +
  geom_point(colour = "#0c4c8a") +
  theme_minimal()

cor(female_common$padj, female_common$FDR, method = "pearson")
cor.test(female_common$padj, female_common$FDR)
ggplot(female_common) +
  aes(x = padj, y = FDR) +
  geom_point(colour = "#0c4c8a") +
  theme_minimal()

cor(male_common$padj, male_common$FDR, method = "pearson")
cor.test(male_common$padj, male_common$FDR)
ggplot(male_common) +
  aes(x = padj, y = FDR) +
  geom_point(colour = "#0c4c8a") +
  theme_minimal()

################################################################################
################################################################################
#                     2.A: Setting up the Biomart Library                      #
################################################################################
#Code adapted from https://igordot.github.io/msigdbr/articles/msigdbr-intro.html
gene_sets_raw <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
#data_cleaning
#remove duplicates
gene_sets <- gene_sets_raw[!duplicated(gene_sets_raw[c("entrez_gene", "gs_name", "gs_exact_source")]), ]
#remove NAs
gene_sets <- gene_sets %>% mutate(pathway_label = paste0(gs_name, " (", gs_exact_source, ")"))

gene_sets_go <- gene_sets %>% group_by(pathway_label) %>% 
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  { set_names(.$genes, .$pathway_label) }

#gene_sets_go <- split(gene_sets_go$entrez_gene, gene_sets_go$gs_name)
#gene_sets <- lapply(gene_sets, function(ids) ids[!is.na(ids)])
#gene_sets <- gene_sets[sapply(gene_sets, length) > 0]

h_mart <- useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))
#h_mart
ens2entrez <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"), mart= h_mart)
head(ens2entrez, 10)

################################################################################
#                                  2.C: FGSEA                                  #
################################################################################
#HGPS Males 1
hgps_males_genes1 <- inner_join(hgps_male_specific_common1, ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
hgps_males_genes1 <- hgps_males_genes1[order(hgps_males_genes1$log2FoldChange, decreasing = TRUE), ]
nrow(hgps_males_genes1)


#clean the data set
#remove NAs
hgps_male_genes_cleaned1 <- hgps_males_genes1[!is.na(hgps_males_genes1$entrezgene_id), ]
nrow(hgps_male_genes_cleaned1)

#collapse duplicates
hgps_male_genes_cleaned1 <- hgps_male_genes_cleaned1[ !duplicated(hgps_male_genes_cleaned1$entrezgene_id), ]
nrow(hgps_male_genes_cleaned1)

#set the row names to entrez gene id
rownames(hgps_male_genes_cleaned1) <- hgps_male_genes_cleaned1$entrezgene_id

hgps_male_gene_names1 <- rownames(hgps_male_genes_cleaned1)
hgps_male_gene_stats1 <- setNames(hgps_male_genes_cleaned1$stat, hgps_male_genes_cleaned1$Gene)

#creating the fgsea table
hgps_males_fgsea1 <- fgsea(pathways = gene_sets_go, 
                          stats = hgps_male_gene_stats1,
                          minSize = 5, 
                          maxSize = length(hgps_male_gene_stats1) - 1,
                          scoreType = "pos")

nrow(hgps_males_fgsea1)
hgps_males_fgsea1 <- hgps_males_fgsea1[order(hgps_males_fgsea1$padj, decreasing = FALSE), ]
hgps_males_fgsea1 <- hgps_males_fgsea1 %>%
  extract(col   = pathway,
          into  = c("term_name", "GO_code"),
          regex = "^(.*) \\((GO:\\d+)\\)$",
          remove = TRUE)
hgps_males_fgsea1 <- hgps_males_fgsea1 %>%
  mutate(term_name = str_remove(term_name, "^GOBP_"))
head(hgps_males_fgsea1, 20)

hgps_male_com_path1_fgsea_plot <- hgps_males_fgsea1 %>%
  slice_head(n = 10) %>%
  mutate(log10_padj = -log10(padj),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_male_com_path1_fgsea_plot, aes(x = term_name, y = log10_padj, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#030e1e", "#1F497D"))(2),  # 20-step gradient
    name = expression(NES)
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, family = "Serif"),
    plot.title = element_text(hjust = 0.5, face = "bold"))


#HGPS Males 2
hgps_males_genes2 <- inner_join(hgps_male_specific_common2, ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
hgps_males_genes2 <- hgps_males_genes2[order(hgps_males_genes2$log2FoldChange, decreasing = TRUE), ]
nrow(hgps_males_genes2)

#clean the data set
#remove NAs
hgps_male_genes_cleaned2 <- hgps_males_genes2[!is.na(hgps_males_genes2$entrezgene_id), ]
hgps_male_genes_cleaned2 <- hgps_male_genes_cleaned2 %>%
  filter(!is.na(Gene), is.finite(stat)) %>%
  distinct(Gene, .keep_all = TRUE) 

#collapse duplicates
hgps_male_genes_cleaned2 <- hgps_male_genes_cleaned2[ !duplicated(hgps_male_genes_cleaned2[c("entrezgene_id", "Gene")]),]
nrow(hgps_male_genes_cleaned2)

#set the row names to entrez gene id
rownames(hgps_male_genes_cleaned2) <- hgps_male_genes_cleaned2$entrezgene_id

hgps_male_gene_names2 <- rownames(hgps_male_genes_cleaned2)
hgps_male_gene_stats2 <- setNames(hgps_male_genes_cleaned2$stat, hgps_male_genes_cleaned2$Gene)


#creating the fgsea table
hgps_males_fgsea2 <- fgsea(pathways = gene_sets_go, 
                           stats = hgps_male_gene_stats2, 
                           minSize = 5, 
                           maxSize = length(hgps_male_gene_stats2) - 1,
                           scoreType = "pos")

nrow(hgps_males_fgsea2)
hgps_males_fgsea2 <- hgps_males_fgsea2[order(hgps_males_fgsea2$padj, decreasing = FALSE), ]
hgps_males_fgsea2 <- hgps_males_fgsea2 %>%
  extract(col   = pathway,
          into  = c("term_name", "GO_code"),
          regex = "^(.*) \\((GO:\\d+)\\)$",
          remove = TRUE)
hgps_males_fgsea2 <- hgps_males_fgsea2 %>%
  mutate(term_name = str_remove(term_name, "^GOBP_"))
head(hgps_males_fgsea2, 20)

hgps_male_com_path2_fgsea_plot <- hgps_males_fgsea2 %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(padj),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_male_com_path2_fgsea_plot, aes(x = term_name, y = log10_padj, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#030e1e", "#1F497D"))(2),  # 20-step gradient
    name = expression(-log[10]("padj"))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, family = "Serif"),
    plot.title = element_text(hjust = 0.5, face = "bold"))


#HGPS Females 1
hgps_fem_genes1 <- inner_join(hgps_fem_specific_common1, ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
hgps_fem_genes1 <- hgps_fem_genes1[order(hgps_fem_genes1$log2FoldChange, decreasing = TRUE), ]
nrow(hgps_fem_genes1)

#clean the data set
#remove NAs
hgps_fem_genes_cleaned1 <- hgps_fem_genes1[!is.na(hgps_fem_genes1$entrezgene_id), ]
nrow(hgps_fem_genes_cleaned1)

#collapse duplicates
hgps_fem_genes_cleaned1 <- hgps_fem_genes_cleaned1[ !duplicated(hgps_fem_genes_cleaned1$entrezgene_id), ]
nrow(hgps_fem_genes_cleaned1)

#set the row names to entrez gene id
rownames(hgps_fem_genes_cleaned1) <- hgps_fem_genes_cleaned1$entrezgene_id

hgps_fem_gene_names1 <- rownames(hgps_fem_genes_cleaned1)
hgps_fem_gene_stats1 <- setNames(hgps_fem_genes_cleaned1$stat, hgps_fem_genes_cleaned1$Gene)

#creating the fgsea table
hgps_fem_fgsea1 <- fgsea(pathways = gene_sets_go, 
                           stats = hgps_fem_gene_stats1,
                           minSize = 5, 
                           maxSize = length(hgps_fem_gene_stats1) - 1,
                           scoreType = "pos")

nrow(hgps_fem_fgsea1)
hgps_fem_fgsea1 <- hgps_fem_fgsea1[order(hgps_fem_fgsea1$padj, decreasing = FALSE), ]
hgps_fem_fgsea1 <- hgps_fem_fgsea1 %>%
  extract(col   = pathway,
          into  = c("term_name", "GO_code"),
          regex = "^(.*) \\((GO:\\d+)\\)$",
          remove = TRUE)
hgps_fem_fgsea1 <- hgps_fem_fgsea1 %>%
  mutate(term_name = str_remove(term_name, "^GOBP_"))
head(hgps_fem_fgsea1, 20)

hgps_fem_com_path1_fgsea_plot <- hgps_fem_fgsea1 %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(padj),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_fem_com_path1_fgsea_plot, aes(x = term_name, y = log10_padj, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#B35E00"))(2),  # 20-step gradient
    name = expression(NES)
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, family = "Serif"),
    plot.title = element_text(hjust = 0.5, face = "bold"))


#HGPS Females 2
hgps_fem_genes2 <- inner_join(hgps_fem_specific_common2, ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
hgps_fem_genes2 <- hgps_fem_genes2[order(hgps_fem_genes2$log2FoldChange, decreasing = TRUE), ]
nrow(hgps_fem_genes2)

#clean the data set
#remove NAs
hgps_fem_genes_cleaned2 <- hgps_fem_genes2[!is.na(hgps_fem_genes2$entrezgene_id), ]
hgps_fem_genes_cleaned2 <- hgps_fem_genes_cleaned2 %>%
  filter(!is.na(Gene), is.finite(stat)) %>%
  distinct(Gene, .keep_all = TRUE) 
nrow(hgps_fem_genes_cleaned2)

#collapse duplicates
hgps_fem_genes_cleaned2 <- hgps_fem_genes_cleaned2[ !duplicated(hgps_fem_genes_cleaned2[c("entrezgene_id", "Gene")]),]
nrow(hgps_fem_genes_cleaned2)

#set the row names to entrez gene id
rownames(hgps_fem_genes_cleaned2) <- hgps_fem_genes_cleaned2$entrezgene_id

hgps_fem_gene_names2 <- rownames(hgps_fem_genes_cleaned2)
hgps_fem_gene_stats2 <- setNames(hgps_fem_genes_cleaned2$stat, hgps_fem_genes_cleaned2$Gene)


#creating the fgsea table
hgps_fem_fgsea2 <- fgsea(pathways = gene_sets_go, 
                           stats = hgps_fem_gene_stats2, 
                           minSize = 5, 
                           maxSize = length(hgps_fem_gene_stats2) -1,
                           scoreType = "pos")

nrow(hgps_fem_fgsea2)
hgps_fem_fgsea2 <- hgps_fem_fgsea2[order(hgps_fem_fgsea2$padj, decreasing = FALSE), ]
hgps_fem_fgsea2 <- hgps_fem_fgsea2 %>%
  extract(col   = pathway,
          into  = c("term_name", "GO_code"),
          regex = "^(.*) \\((GO:\\d+)\\)$",
          remove = TRUE)
hgps_fem_fgsea2 <- hgps_fem_fgsea2 %>%
  mutate(term_name = str_remove(term_name, "^GOBP_"))
head(hgps_fem_fgsea2, 20)

hgps_fem_com_path2_fgsea_plot <- hgps_fem_fgsea2 %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(padj),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_fem_com_path2_fgsea_plot, aes(x = term_name, y = log10_padj, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#FF8C42","#B35E00"))(2),  # 20-step gradient
    name = expression(NES)
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, family = "Serif"),
    plot.title = element_text(hjust = 0.5, face = "bold"))




# hgps_females_genes2 <- inner_join(hgps_fem_specific_common2, ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
# hgps_females_genes2 <- hgps_females_genes2[order(hgps_females_genes2$padj, decreasing = FALSE), ]
# nrow(hgps_females_genes2)
# 
# #clean the data set
# #remove NAs
# hgps_female_genes_cleaned2 <- hgps_females_genes2[!is.na(hgps_females_genes2$entrezgene_id), ]
# nrow(hgps_female_genes_cleaned2)
# 
# #collapse duplicates
# hgps_female_genes_cleaned2 <- hgps_female_genes_cleaned2[ !duplicated(hgps_female_genes_cleaned2$entrezgene_id), ]
# nrow(hgps_female_genes_cleaned2)
# 
# #set the row names to entrez gene id
# rownames(hgps_female_genes_cleaned2) <- hgps_female_genes_cleaned2$entrezgene_id
# 
# hgps_female_gene_names2 <- rownames(hgps_female_genes_cleaned2)
# hgps_female_gene_stats2 <- setNames(hgps_female_genes_cleaned2$stat, hgps_female_genes_cleaned2$entrezgene_id)
# hgps_female_gene_stats2 <- hgps_female_gene_stats2[is.finite(hgps_female_gene_stats2)]
# 
# #creating the fgsea table
# hgps_females_fgsea2 <- fgsea(pathways = gene_sets, 
#                              stats = hgps_female_gene_stats2, 
#                              minSize = 1, 
#                              maxSize = 500)
# head(hgps_females_fgsea2, 10)
# hgps_females_fgsea2 <- hgps_females_fgsea2[order(hgps_females_fgsea2$padj, decreasing = FALSE), ]
# nrow(hgps_females_fgsea2)


#Visualising the data
#male pathways
hgps_male_path_2<- ggplot(hgps_males_fgsea2, aes(x = reorder(pathway, -padj), 
                                                        y = padj,
                                                        fill = padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(name = "p.adj", 
                       option = "magma") +
  labs(x = NULL, 
       y = "Normalized Enrichment Score",
       title = "Top 20 Enriched Pathways in HGPS Males in Order of Decreasing Significance", 
       subtitle = "System Analysed: FGSEA") +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                                        axis.title.x    = element_text(size = 14),
                                        axis.title.y    = element_text(size = 14),
                                        axis.text.x     = element_text(size = 12),
                                        axis.text.y     = element_text(size = 12),
                                        legend.title    = element_text(size = 13),
                                        legend.text     = element_text(size = 11))

#female pathways
hgps_female_path_2 <- ggplot(hgps_females_fgsea2, aes(x = reorder(pathway, -padj), 
                                                             y = padj,
                                                             fill = padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(name = "p.adj", 
                       option = "magma") +
  labs(x = NULL, 
       y = "Normalized Enrichment Score",
       title = "Top 20Enriched Pathways in HGPS Females in Order of Decreasing Significance",
       subtitle = "System Analysed: FGSEA") +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                                        axis.title.x    = element_text(size = 14),
                                        axis.title.y    = element_text(size = 14),
                                        axis.text.x     = element_text(size = 12),
                                        axis.text.y     = element_text(size = 12),
                                        legend.title    = element_text(size = 13),
                                        legend.text     = element_text(size = 11))
(hgps_male_path_2 | hgps_female_path_2)
plot_annotation(
  title = "Comparison of Enriched Pathways in HGPS Males and HGPS Females)",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)

#HGPS tottal 1
hgps_genes <- inner_join(subset(hgps_common, padj < 0.05 & FDR < 0.05), ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
hgps_genes <- hgps_genes[order(hgps_genes$log2FoldChange, decreasing = TRUE), ]
nrow(hgps_genes)

#clean the data set
#remove NAs
hgps_genes_cleaned <- hgps_genes[!is.na(hgps_genes$entrezgene_id), ]
nrow(hgps_genes_cleaned)

#collapse duplicates
hgps_genes_cleaned <- hgps_genes_cleaned[ !duplicated(hgps_genes_cleaned$entrezgene_id), ]
nrow(hgps_genes_cleaned)

#set the row names to entrez gene id
rownames(hgps_genes_cleaned) <- hgps_genes_cleaned$entrezgene_id

hgps_gene_names <- rownames(hgps_genes_cleaned)
hgps_gene_stats <- setNames(hgps_genes_cleaned$stat, hgps_genes_cleaned$entrezgene_id)
hgps_gene_stats <- hgps_gene_stats[is.finite(hgps_male_stats)]

#creating the fgsea table
hgps_fgsea1 <- fgsea(pathways = gene_sets, 
                           stats = hgps_gene_stats, 
                           minSize = 1, 
                           maxSize = 500)

head(hgps_fgsea1, 20)
hgps_fgsea <- hgps_fgsea1[order(hgps_fgsea1$padj, decreasing = FALSE), ]


################################################################################
#                                  2.D: enrichR                                #
################################################################################
enrichrdbs <- listEnrichrDbs()
colnames(enrichrdbs)
enrichrdbs <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025", "GO_Biological_Process_2025")

#investigating the unique hgps male genes
#first we need to only save the list of genes - we don't want stats
hgps_males_specific_common_list1 <- hgps_male_specific_common1$Gene

enriched_hgps_males1 <- enrichr(hgps_males_specific_common_list1, enrichrdbs)
head(enriched_hgps_males1)

e_hgps_males_processes1 <- enriched_hgps_males1[["GO_Biological_Process_2025"]]
e_hgps_males_processes1$Count <- as.numeric(sub("/.*", "", e_hgps_males_processes1$Overlap))
e_hgps_males_processes1 <- e_hgps_males_processes1[order(e_hgps_males_processes1$Adjusted.P.value, decreasing = FALSE), ]
View(e_hgps_males_processes1)

e_hgps_males_processes1 <- e_hgps_males_processes1 %>%
  separate(Term, into = c("Pathway", "GO_code"), sep = " \\(", remove = TRUE) %>%
  mutate(GO_code = sub("\\)", "", GO_code))
top_20_hgps_male_enrichr1 <- e_hgps_males_processes1[1:20, ]

#investigating the unique hgps female genes
#first we need to only save the list of genes - we don't want stats
hgps_females_specific_common_list1 <- hgps_fem_specific_common1$Gene

enriched_hgps_females1 <- enrichr(hgps_females_specific_common_list1, enrichrdbs)
head(enriched_hgps_females1)

head(enriched_hgps_females1[["GO_Biological_Process_2025"]])
e_hgps_females_processes1 <- enriched_hgps_females1[["GO_Biological_Process_2025"]]
e_hgps_females_processes1$Count <- as.numeric(sub("/.*", "", e_hgps_females_processes1$Overlap))
e_hgps_females_processes1 <- e_hgps_females_processes1[order(e_hgps_females_processes1$Adjusted.P.value, decreasing = FALSE), ]
e_hgps_females_processes1 <- e_hgps_females_processes1 %>%
  separate(Term, into = c("Pathway", "GO_code"), sep = " \\(", remove = TRUE) %>%
  mutate(GO_code = sub("\\)", "", GO_code))

#disease effect
hgps_males_specific_common_list2 <- hgps_male_specific_common2$Gene
hgps_males_specific_list2 <- hgps_male_genes_cleaned2$Gene

enriched_hgps_males2 <- enrichr(hgps_males_specific_common_list2, enrichrdbs)
head(enriched_hgps_males2)

e_hgps_males_processes2 <- enriched_hgps_males2[["GO_Biological_Process_2025"]]
e_hgps_males_processes2$Count <- as.numeric(sub("/.*", "", e_hgps_males_processes2$Overlap))
e_hgps_males_processes2 <- e_hgps_males_processes2[order(e_hgps_males_processes2$Adjusted.P.value, decreasing = FALSE), ]
View(e_hgps_males_processes2)

e_hgps_males_processes2 <- e_hgps_males_processes2 %>%
  separate(Term, into = c("Pathway", "GO_code"), sep = " \\(", remove = TRUE) %>%
  mutate(GO_code = sub("\\)", "", GO_code))

#investigating the unique hgps female genes
hgps_females_specific_common_list2 <- hgps_fem_specific_common2$Gene

enriched_hgps_females2 <- enrichr(hgps_females_specific_common_list2, enrichrdbs)
head(enriched_hgps_females2)

head(enriched_hgps_females2[["GO_Biological_Process_2025"]])
e_hgps_females_processes2 <- enriched_hgps_females2[["GO_Biological_Process_2025"]]
e_hgps_females_processes2$Count <- as.numeric(sub("/.*", "", e_hgps_females_processes2$Overlap))
e_hgps_females_processes2 <- e_hgps_females_processes2[order(e_hgps_females_processes2$Adjusted.P.value, decreasing = FALSE), ]
e_hgps_females_processes2 <- e_hgps_females_processes2 %>%
  separate(Term, into = c("Pathway", "GO_code"), sep = " \\(", remove = TRUE) %>%
  mutate(GO_code = sub("\\)", "", GO_code))


#plotting the data
#males
# hgps_male_enricher_1 <- ggplot(top_20_hgps_male_enrichr, aes(x = reorder(Term, -Adjusted.P.value), 
#                                                              y = Adjusted.P.value,
#                                                              fill = Adjusted.P.value)) +
#   geom_col() +
#   coord_flip() +
#   scale_fill_viridis_c(name = "p-value", option = "viridis") +
#   labs(x = NULL, 
#        y = "p-value",
#        title = "Top 20 Enriched Pathways in HGPS Males in Order of Decreasing Significance",
#        subtitle = "System Analysed: enrichR") +
#   theme_minimal(base_size = 12) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5),
#         axis.text.y = element_text(size = 10))
# 
# hgps_male_enricher_1
# 
# 
# #females
# hgps_female_enricher_1 <- ggplot(top_20_hgps_female_enrichr, aes(x = reorder(Term, -Adjusted.P.value), 
#                                                                  y = Adjusted.P.value,
#                                                                  fill = Adjusted.P.value)) +
#   geom_col() +
#   coord_flip() +
#   scale_fill_viridis_c(name = "p-value", option = "viridis") +
#   labs(x = NULL, 
#        y = "p-value",
#        title = "Top 20 Enriched Pathways in HGPS Females in Order of Decreasing Significance",
#        subtitle = "System Analysed: enrichR") +
#   theme_minimal(base_size = 12) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5),
#         axis.text.y = element_text(size = 10))
# 
# hgps_female_enricher_1

################################################################################
#                                 2.E: gprofiler                                #
################################################################################
#as another enrichment pathway comparison we are using g profiler. 

#first investigating genes only expressed in females 1
hgps_fem_prof1 <- gost(hgps_fem_specific_common1$Gene, organism = "hsapiens", significant = TRUE, correction_method = "fdr")
names(hgps_fem_prof1)
hgps_fem_prof_res1 <- hgps_fem_prof1$result
hgps_fem_prof_res1 <- hgps_fem_prof_res1[order(hgps_fem_prof_res1$source, hgps_fem_prof_res1$p_value, decreasing = FALSE), ]

hgps_fem_prof_meta1 <- hgps_fem_prof1$meta

gostplot(hgps_fem_prof_res1, capped = TRUE, interactive = TRUE)

hgps_fem_go1 <- hgps_fem_prof_res1 %>% filter(grepl("^GO:", term_id))



#first investigating genes only expressed in males
hgps_male_prof1 <- gost(hgps_male_specific_common1$Gene, organism = "hsapiens", significant = TRUE, correction_method = "fdr")
colnames(hgps_male_prof1$result)
hgps_male_prof_res1 <- hgps_male_prof1$result
hgps_male_prof_res1 <- hgps_male_prof_res1[order(hgps_male_prof_res1$source, hgps_male_prof_res1$p_value, decreasing = FALSE), ]

hgps_male_prof_meta1 <- hgps_male_prof1$meta

gostplot(hgps_male_prof1, capped = TRUE, interactive = TRUE)

hgps_male_go1 <- hgps_male_prof_res1 %>% filter(grepl("^GO:", term_id))


#finding the pathways unique to hgps females
hgps_females_prof_unique <- hgps_fem_prof_res[!hgps_fem_prof_res$term_name %in% hgps_male_prof_res$term_name, ]
nrow(hgps_females_prof_unique)
hgps_females_prof_unique<- hgps_females_prof_unique[order(hgps_females_prof_unique$p_value, decreasing = FALSE), ]
top_20_hgps_female_prof <- hgps_females_prof_unique[1:20, ]



#finding the pathways unique to hgps males
hgps_males_prof_unique <- hgps_male_prof_res[!hgps_male_prof_res$term_name %in% hgps_fem_prof_res$term_name, ]
nrow(hgps_males_prof_unique)
hgps_males_prof_unique<- hgps_males_prof_unique[order(hgps_males_prof_unique$p_value, decreasing = FALSE), ]
top_20_hgps_male_prof <- hgps_males_prof_unique[1:20, ]

#first investigating genes only expressed in females 1
hgps_fem_prof2 <- gost(hgps_fem_specific_common2$Gene, organism = "hsapiens", significant = TRUE, correction_method = "fdr")
names(hgps_fem_prof2)
hgps_fem_prof_res2 <- hgps_fem_prof2$result
hgps_fem_prof_res2 <- hgps_fem_prof_res2[order(hgps_fem_prof_res2$source, hgps_fem_prof_res2$p_value, decreasing = FALSE), ]

hgps_fem_prof_meta2 <- hgps_fem_prof2$meta

gostplot(hgps_fem_prof_res2, capped = TRUE, interactive = TRUE)

hgps_fem_go2 <- hgps_fem_prof_res2 %>% filter(grepl("^GO:", term_id))

#first investigating genes only expressed in males
hgps_male_prof2 <- gost(hgps_male_specific_common2$Gene, organism = "hsapiens", significant = TRUE, correction_method = "fdr")
colnames(hgps_male_prof2$result)
hgps_male_prof_res2 <- hgps_male_prof2$result
hgps_male_prof_res2 <- hgps_male_prof_res2[order(hgps_male_prof_res2$source, hgps_male_prof_res2$p_value, decreasing = FALSE), ]

hgps_male_prof_meta2 <- hgps_male_prof2$meta

gostplot(hgps_male_prof2, capped = TRUE, interactive = TRUE)

hgps_male_go2 <- hgps_male_prof_res2 %>% filter(grepl("^GO:", term_id))


#plotting the data
# hgps_male_prof_1 <- ggplot(top_20_hgps_male_prof, aes(x = reorder(term_name, p_value), 
#                                                       y = p_value,
#                                                       fill = p_value)) +
#   geom_col() +
#   coord_flip() +
#   scale_fill_viridis_c(name = "p-value", option = "viridis") +
#   labs(x = NULL, 
#        y = "p-value",
#        title = "Top 20 Enriched Pathways in HGPS Males in Order of Decreasing Significance",
#        subtitle = "System analysed: g:profiler") +
#   theme_minimal(base_size = 12) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5),
#         axis.text.y = element_text(size = 10))
# 
# hgps_male_prof_1
# 
# 
# hgps_female_prof_1 <- ggplot(top_20_hgps_female_prof, aes(x = reorder(term_name, p_value), 
#                                                           y = p_value,
#                                                           fill = p_value)) +
#   geom_col() +
#   coord_flip() +
#   scale_fill_viridis_c(name = "p-value", option = "viridis") +
#   labs(x = NULL, 
#        y = "p-value",
#        title = "Top 20 Enriched Pathways in HGPS Females in Order of Decreasing p-value",
#        subtitle = "System Analysed: g:profiler") +
#   theme_minimal(base_size = 12) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5),
#         axis.text.y = element_text(size = 10))
# hgps_female_prof_1 


###
# comparing pathways
#hgps females (sex difference)
hgps_fem_com_path1 <- inner_join(hgps_fem_go1, e_hgps_females_processes1, by = c("term_id" = "GO_code"))
hgps_fem_com_path1 <- hgps_fem_com_path1[order(hgps_fem_com_path1$Adjusted.P.value, decreasing = FALSE), ]
hgps_fem_com_path1 <- hgps_fem_com_path1[!duplicated(hgps_fem_com_path1$term_name), ]
nrow(hgps_fem_com_path1)
hgps_fem_com_path1$Pathway[1:16]

hgps_fem_com_path1_plot <- hgps_fem_com_path1 %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 10) %>%
  mutate(log10_padj = -log10(Adjusted.P.value),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_fem_com_path1_plot, aes(x = term_name, y = log10_padj, fill = log10_padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#B35E00"))(10),  # 20-step gradient
    name = expression(-log[10]("adjusted p-value"))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("adjusted p-value")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title  = element_text(size = 20),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.text = element_text(size = 20),
    legend.title= element_text(size = 20))

#hpgs males (sex difference
hgps_male_com_path1 <- inner_join(hgps_male_go1, e_hgps_males_processes1, by = c("term_id" = "GO_code"))
hgps_male_com_path1 <- hgps_male_com_path1[order(hgps_male_com_path1$Adjusted.P.value, decreasing = FALSE), ]
hgps_male_com_path1 <- hgps_male_com_path1[!duplicated(hgps_male_com_path1$term_name), ]
nrow(hgps_male_com_path1)

hgps_male_com_path1_plot <- hgps_male_com_path1 %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(Adjusted.P.value),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_male_com_path1_plot, aes(x = term_name, y = log10_padj, fill = log10_padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#030e1e", "#1F497D"))(20),  # 20-step gradient
    name = expression(-log[10]("adjusted p-value"))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("adjusted p-value")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, family = "Serif"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#disease effect
hgps_fem_com_path2 <- inner_join(hgps_fem_go2, e_hgps_females_processes2, by = c("term_id" = "GO_code"))
hgps_fem_com_path2 <- hgps_fem_com_path2[order(hgps_fem_com_path2$Adjusted.P.value, decreasing = FALSE), ]
hgps_fem_com_path2 <- hgps_fem_com_path2[!duplicated(hgps_fem_com_path2$term_name), ]
nrow(hgps_fem_com_path2)

hgps_fem_com_path2_plot <- hgps_fem_com_path2 %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(Adjusted.P.value),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_fem_com_path2_plot, aes(x = term_name, y = log10_padj, fill = log10_padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#FF8C42", "#B35E00"))(2),  # 20-step gradient
    name = expression(-log[10]("padj"))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title  = element_text(size = 12),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12))

#hpgs males (sex difference
hgps_male_com_path2 <- inner_join(hgps_male_go2, e_hgps_males_processes2, by = c("term_id" = "GO_code"))
hgps_male_com_path2 <- hgps_male_com_path2[order(hgps_male_com_path2$Adjusted.P.value, decreasing = FALSE), ]
hgps_male_com_path2 <- hgps_male_com_path2[!duplicated(hgps_male_com_path2$term_name), ]
nrow(hgps_male_com_path2)


hgps_male_com_path2_plot <- hgps_male_com_path2 %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(Adjusted.P.value),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(hgps_male_com_path2_plot, aes(x = term_name, y = log10_padj, fill = log10_padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#030e1e", "#1F497D"))(2),  # 20-step gradient
    name = expression(-log[10]("padj"))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, family = "Serif"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

################################################################################
# Gene Network Analysis
################################################################################
GSE113957_dseq_cpm #matrix

hgps_male_specific_common1 #stats
hgps_fem_specific_common1 #stats
hgps_male_specific_common2 #stats
hgps_fem_specific_common2 #stats

lmna_col_hgps <- hgps_common["LMNA", ]   
hgps_male_specific_lmna1 <- rbind(hgps_male_specific_common1, lmna_col_hgps)
hgps_fem_specific_lmna1 <- rbind(hgps_fem_specific_common1, lmna_col_hgps)

PearsCorr <- sapply(setdiff(colnames(), "LMNA"),
                    function(x) cor(tra_df[, x], tra_df$LMNA, method = "pearson"))

lmna_col_fem <- female_common["LMNA", ] 
hgps_fem_specific_lmna2 <- rbind(hgps_fem_specific_common2, lmna_col_hgps)

lmna_col_male <- male_common["LMNA", ] 
hgps_male_specific_lmna2 <- rbind(hgps_male_specific_common2, lmna_col_hgps)


#### string ####
################################################################################
#                                 3.A: STRINGdb                                #
################################################################################
#code adapted from: https://www.rdocumentation.org/packages/STRINGdb/versions/1.8.1/topics/STRINGdb

#creating a string object
string <- STRINGdb$new(version="12.0", 
                       species=9606, #NCBI taxonomic ID for homo sapiens
                       score_threshold=400, 
                       input_directory="") 

#hgps females sex effect)
hgps_females_mapped1 <- string$map(hgps_fem_specific_lmna1, "Gene", removeUnmappedRows =  TRUE)

#getting the hits
#getting the hits
hgps_female_hits1 <- hgps_females_mapped1$STRING_id
#plotting the hits
string$plot_network(hgps_female_hits1)

#because there's a lot of genes here I'm going to focus on cluster interactions to see which genes LMNA is properly interacting with
#obtaining the individual clusters
hgps_female_clusters1 <- string$get_clusters(hgps_female_hits1)
#investigating the clusters
sapply(hgps_female_clusters1, length)
#from this we can see we have 13 proper clusters


#hgps_males (sex effect)
#First I'm going to add LMNA to the end of this list because I want to see if LMNA interacts with any of the gens
#changing the data frame to just hold the gene column
hgps_males_mapped1 <- string$map(hgps_male_specific_lmna1, "Gene", removeUnmappedRows =  TRUE)

#getting the hits
#getting the hits
hgps_male_hits1 <- hgps_males_mapped1$STRING_id
#plotting the hits
string$plot_network(hgps_male_hits1)

#because there's a lot of genes here I'm going to focus on cluster interactions to see which genes LMNA is properly interacting with
#obtaining the individual clusters
hgps_male_clusters1 <- string$get_clusters(hgps_male_hits1)
#investigating the clusters
sapply(hgps_male_clusters1, length)
#from this we can see we have 13 proper clusters

#converting the string ids to gene names
hgps_malegene_map1 <- hgps_males_mapped1[, c("STRING_id","Gene")]
hgps_males_ids2symb <- function(v) {m <- hgps_malegene_map1$Gene[match(v, hgps_malegene_map1$STRING_id)]
ifelse(is.na(m), v, m)
}

#saving cluster 2 to a new variable because it contains lmna
male_hgps_lmna_cluster_a <- hgps_males_ids2symb(hgps_male_clusters1[[1]])
male_hgps_lmna_cluster_b <- hgps_males_ids2symb(hgps_male_clusters1[[2]])
male_hgps_lmna_cluster_c <- hgps_males_ids2symb(hgps_male_clusters1[[3]])

string$plot_network(male_hgps_lmna_cluster_a)
string$plot_network(male_hgps_lmna_cluster_b)
string$plot_network(male_hgps_lmna_cluster_c)

#hgps females (disease effect)
hgps_fem_mapped2 <- string$map(hgps_fem_specific_lmna2, "Gene", removeUnmappedRows =  TRUE)

#getting the hits
hgps_fem_hits2 <- hgps_fem_mapped2$STRING_id
#plotting the hits
string$plot_network(hgps_fem_hits2)

#because there's a lot of genes here I'm going to focus on cluster interactions to see which genes LMNA is properly interacting with
#obtaining the individual clusters
hgps_fem_clusters2 <- string$get_clusters(hgps_fem_hits2)
#investigating the clusters
sapply(hgps_fem_clusters2, length)
#from this we can see we have 13 proper clusters

#converting the string ids to gene names
hgps_femgene_map2 <- hgps_fem_mapped2[, c("STRING_id","Gene")]
hgps_fem_ids2symb <- function(v) {m <- hgps_femgene_map2$Gene[match(v, hgps_femgene_map2$STRING_id)]
ifelse(is.na(m), v, m)
}

#investigating the individual clusters for lmna
if ("LMNA" %in% hgps_fem_ids2symb(hgps_fem_clusters2[[2]])){
  print("True")
} else {
  print("False")
}

#saving cluster 2 to a new variable because it contains lmna
fem_hgps_lmna_cluster_a <- hgps_males_ids2symb(hgps_fem_clusters2[[2]])

string$plot_network(fem_hgps_lmna_cluster_a)

hgps_fem2_lmna_interactions <- string$get_interactions(hgps_fem_clusters2[[2]])
hgps_fem2_lmna_interactions

#converting the string ids to gene names again
hgps_fem_ids2symb(hgps_fem2_lmna_interactions)


#hgps males (disease effect)
hgps_males_mapped2 <- string$map(hgps_male_specific_lmna2, "Gene", removeUnmappedRows =  TRUE)

#getting the hits
#getting the hits
hgps_male_hits2 <- hgps_males_mapped2$STRING_id
#plotting the hits
string$plot_network(hgps_male_hits2)

#because there's a lot of genes here I'm going to focus on cluster interactions to see which genes LMNA is properly interacting with
#obtaining the individual clusters
hgps_male_clusters2 <- string$get_clusters(hgps_male_hits2)
#investigating the clusters
sapply(hgps_male_clusters2, length)
#from this we can see we have 13 proper clusters

#converting the string ids to gene names
hgps_malegene_map2 <- hgps_males_mapped2[, c("STRING_id","Gene")]
hgps_males2_ids2symb <- function(v) {m <- hgps_malegene_map2$Gene[match(v, hgps_malegene_map2$STRING_id)]
ifelse(is.na(m), v, m)
}

#investigating the individual clusters for lmna
if ("LMNA" %in% hgps_males2_ids2symb(hgps_male_clusters2[[5]])){
  print("True")
} else {
  print("False")
}

#saving cluster 2 to a new variable because it contains lmna
male_hgps_lmna_cluster_a <- hgps_males2_ids2symb(hgps_male_clusters2[[5]])

string$plot_network(male_hgps_lmna_cluster_a)



################################################################################
# DCM Data Analysis
################################################################################
male_data_loc <- "redacted"
female_data_loc <- "redacted"
dcm_comb_meta_loc <- "redacted"

#cleaning the meta data 
dcm_comb_meta <- as.data.frame(read_csv(dcm_comb_meta_loc, col_names = TRUE))
rownames(dcm_comb_meta) <- dcm_comb_meta$Sample
dcm_comb_meta <- subset(dcm_comb_meta, select = -Sample)
View(dcm_comb_meta)

#cleaning the data
fem_dcm_data <- as.data.frame(read_csv(female_data_loc, col_names = TRUE))
rownames(fem_dcm_data) <- fem_dcm_data$...1
#removing the gene column
fem_dcm_data <- subset(fem_dcm_data, select = -...1)
nrow(fem_dcm_data)

male_dcm_data <- as.data.frame(read_csv(male_data_loc, col_names = TRUE))
#cleaning up the data frame
rownames(male_dcm_data) <- male_dcm_data$...1
male_dcm_data <- subset(male_dcm_data, select = -...1)
nrow(male_dcm_data)

#checking if the rows between males and females are identical
length(setdiff(rownames(fem_dcm_data), rownames(male_dcm_data)))
#combining the two data frames
dcm_data_comb <- cbind(fem_dcm_data, male_dcm_data)
dcm_data_comb <- as.matrix(dcm_data_comb)

#joining the matrix and meta data
dcm_idx <- match(colnames(dcm_data_comb), rownames(dcm_comb_meta))
#ordering the meta data to match the 
dcm_comb_meta <- dcm_comb_meta[dcm_idx,]

#Creating a new group column
dcm_comb_meta$Group <- with(dcm_comb_meta, paste(Condition, Sex, sep = "_"))

#assigning colours to the groups
dcm_colours <- c(Control_Male = "#4A90E2", 
                 DCM_Male = "#008B8B",
                 Control_Female = "#FF8C42",
                 DCM_Female = "#7B4173") 
dcm_comb_meta$Group_Colour <- dcm_colours[dcm_comb_meta$Group]
#rounding the dcm data to full numbers
dcm_data_comb <- round(dcm_data_comb)


################################################################################
#                                                                              #
#                            PART 1: DEG ANALYSIS                              #
#                                                                              #
################################################################################

################################################################################
#                            1.A: DEseq2 Analysis                              #
################################################################################
dcm_dsq <- DESeqDataSetFromMatrix(countData = dcm_data_comb,
                                  colData = dcm_comb_meta, 
                                  design = ~ Group)

#cleaning the data so that it has sufficient counts
dcm_dsq <-dcm_dsq[rowSums(counts(dcm_dsq) >= 10) >=3, ]
nrow(dcm_dsq)

dcm_dsq <- estimateSizeFactors(dcm_dsq)
sizeFactors(dcm_dsq)

dcm_dsq_log2 <- log2(1+counts(dcm_dsq, normalized=TRUE))
#hgps_log2
#fpm is similar to cpm and discussed here https://rdrr.io/bioc/DESeq2/man/fpm.html.  Here the data is log transformed
dcm_dsq_fpm <-log(fpm(dcm_dsq)+1)#this is also normalised by the sizeFactors
#rlog ad vsd are discussed here https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
dcm_dsq_vsd <- vst(dcm_dsq, blind = T)
#head(hgps_vsd, 3)
dcm_dsq_rld <- rlog(dcm_dsq, blind = T)
dcm_dsq_cpm <- fpm(dcm_dsq, robust = TRUE)
dcm_dsq_cpm_cols <- colnames(dcm_dsq_cpm)

#making pca plots
dcm_pca_dsq <- prcomp(t(assay(dcm_dsq_vsd)), scale.=TRUE)

dcm_pca_dsq_3d <- data.frame(Sample = rownames(dcm_pca_dsq$x),
                             PC1    = dcm_pca_dsq$x[,1],
                             PC2    = dcm_pca_dsq$x[,2],
                             PC3    = dcm_pca_dsq$x[,3],
                             Group  = dcm_comb_meta$Group)

plot_ly(data = dcm_pca_dsq_3d,
        x = ~PC1,
        y = ~PC2,
        z = ~PC3,
        color = ~Group,
        colors = dcm_colours,
        text = ~paste("Sample:", Sample),
        hoverinfo = "text",
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 4)) %>% layout(title = "3D PCA of DCM Patients and Healthy Controls",
                                            scene = list(xaxis = list(title = "PC1"),
                                                         yaxis = list(title = "PC2"),
                                                         zaxis = list(title = "PC3")))
#from this pca there are some prettu drastic batch effects that need to be corrected

#creating the DESeq2 contrasts
dcm_effects <- DESeq(dcm_dsq)
resultsNames(dcm_effects)

# Creating contrasts for the analysis
dcm_control_res <- results(dcm_effects, name = "Group_Control_Male_vs_Control_Female")
dcm_res <- results(dcm_effects, contrast = c("Group", "DCM_Male", "DCM_Female"))
dcm_male_res <- results(dcm_effects, contrast = c("Group", "DCM_Male", "Control_Male"))
dcm_female_res <- results(dcm_effects, name = "Group_DCM_Female_vs_Control_Female")


# Analysis between control males and control females #
#Filter for significant results
cntrl_male_vs_fem <- as.data.frame(subset(dcm_control_res, !is.na(padj)))
cntrl_male_vs_fem$Gene <- rownames(cntrl_male_vs_fem )

ks.test(cntrl_male_vs_fem$log2FoldChange, "pnorm")
#because the log2fold change distribution is not normal, I am going to use rank z_scores because that transforms the data into a normal distribution

cntrl_male_vs_fem$rank_zscore <- qnorm(rank(cntrl_male_vs_fem$log2FoldChange, na.last = "keep") / 
                                         (sum(!is.na(cntrl_male_vs_fem$log2FoldChange)) + 1))

cntrl_male_vs_fem_idx <- colData(dcm_dsq)$Group[c("Control_male", "Control_female")]

cntrl_male_vs_fem_sig_dsq <- subset(cntrl_male_vs_fem, padj < 0.05)

#isolating the genes upregu dcm_dsq#isolating the genes upregulated in males
cntrl_male_dsq_up1_dsq <- subset(cntrl_male_vs_fem, padj < 0.05 & log2FoldChange >= 0.58)
dcm_cntrl_male_up_deg1 <- nrow(cntrl_male_dsq_up1_dsq)

#isolating the genes upregulated in females
cntrl_fem_dsq_up1 <- subset(cntrl_male_vs_fem, padj < 0.05 & log2FoldChange < -0.58)
head(cntrl_fem_dsq_up1, 10)
dcm_cntrl_fem_up_deg1 <- nrow(cntrl_fem_dsq_up1)

dcm_cntrl_male_fem_tot <- dcm_cntrl_male_up_deg1 + dcm_cntrl_fem_up_deg1 

dcm_cntrl_deg_counts_dsq <- data.frame(Category = c("Total", "Control Male-Biased", "Control Female-Biased"),
                                         Count = c(dcm_cntrl_male_fem_tot, dcm_cntrl_male_up_deg1, dcm_cntrl_fem_up_deg1))
dcm_cntrl_deg_counts_dsq$Category <- factor(dcm_cntrl_deg_counts_dsq$Category, levels = c("Total", "Control Male-Biased", "Control Female-Biased"))


ggplot(dcm_cntrl_deg_counts_dsq, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7500)) +
  scale_y_continuous(breaks = seq(0, 7500, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between Control Males and Control Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "Control Male-Biased" = "#4A90E2", "Control Female-Biased" = "#FF8C42", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

################################################
# Analysis between dcm males and dcm females   #
################################################

#Filter for significant results
dcm_male_vs_female <- as.data.frame(subset(dcm_res, !is.na(padj)))
dcm_male_vs_female$Gene <- rownames(dcm_male_vs_female)

#checking the normality of the log2fold values
ks.test(rank(dcm_male_vs_female$log2FoldChange, na.last = "keep"), "pnorm")
#using ranked z_scores because the log2fold change is nor normally distributed
dcm_male_vs_female$rank_zscore <- qnorm(rank(dcm_male_vs_female$log2FoldChange, na.last = "keep") / 
                                          (sum(!is.na(dcm_male_vs_female$log2FoldChange)) + 1))

hist(dcm_male_vs_female$rank_zscore, breaks = 50, main = "Histogram of log2FC Distribution of the DCM Male vs Female", xlab = "log2 Fold Change")


dcm_male_vs_female_idx <- colData(dcm_dsq)$Group[c("DCM_Male", "DCM_Female")]
dcm_male_vs_female_idx <- colData(dcm_dsq)$Group %in% c("DCM_Male", "DCM_Female")

dcm_male_vs_female_sig <- subset(dcm_male_vs_female, padj < 0.05)

#isolating the genes upregulated in males
dcm_male_dsq_up1 <- subset(dcm_male_vs_female, padj < 0.05 & log2FoldChange >= 0.58)
head(dcm_male_dsq_up1, 10)
dcm_male_deg1_dsq <- nrow(dcm_male_dsq_up1)

#isolating the genes upregulated in dcm females
dcm_fem_dsq_up1 <- subset(dcm_male_vs_female, padj < 0.05 & log2FoldChange < -0.58)
head(dcm_fem_dsq_up1, 10)
dcm_fem_deg1_dsq <- nrow(dcm_fem_dsq_up1)

dcm_male_fem_deg_tot_dsq <- dcm_male_deg1_dsq  + dcm_fem_deg1_dsq

dcm_deg_counts_dsq <- data.frame(Category = c("Total", "DCM Male-Biased", "DCM Female-Biased"),
                                       Count = c(dcm_male_fem_deg_tot_dsq , dcm_male_deg1_dsq, dcm_fem_deg1_dsq))
dcm_deg_counts_dsq$Category <- factor(dcm_deg_counts_dsq$Category, levels = c("Total", "DCM Male-Biased", "DCM Female-Biased"))


ggplot(dcm_deg_counts_dsq, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7500)) +
  scale_y_continuous(breaks = seq(0, 7500, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between DCM Males and DCM Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "DCM Male-Biased" = "#008B8B", "DCM Female-Biased" = "#7B4173", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

#############################################
# Analysis between dcm and control females #
#############################################
#Filter for significant results
dcm_fem <- as.data.frame(subset(dcm_female_res, !is.na(padj)))
dcm_fem$Gene <- rownames(dcm_fem)

#checking the normality of the log2fold values
ks.test(rank(dcm_fem$log2FoldChange, na.last = "keep"), "pnorm")
#using ranked z_scores because the log2fold change is nor normally distributed
dcm_fem$rank_zscore <- qnorm(rank(dcm_fem$log2FoldChange, na.last = "keep") / 
                               (sum(!is.na(dcm_fem$log2FoldChange)) + 1))

hist(dcm_fem$rank_zscore, breaks = 50, main = "Histogram of log2FC Distribution of the DCM Female vs Control Female", xlab = "log2 Fold Change")

dcm_female_idx <- colData(dcm_dsq)$Group[c("DCM_Female", "Control_Female")]
dcm_female_idx <- colData(dcm_dsq)$Group %in% c("DCM_Female", "Control_Female")

#isolating the genes upregulated in dcm females
dcm_fem_dsq_up2 <- subset(dcm_fem, padj < 0.05 & log2FoldChange >= 0.58)
head(dcm_fem_dsq_up2, 10)
dcm_fem_deg2_desq <- nrow(dcm_fem_dsq_up2)

#isolating the genes upregulated in females
cntrl_fem_dsq_up2 <- subset(dcm_fem, padj < 0.05 & log2FoldChange < -0.58)
head(cntrl_fem_dsq_up2, 10)
cntrl_fem_deg2_desq <- nrow(cntrl_fem_dsq_up2)

fem_dcm_vs_cntrl_deg_dsq <- dcm_fem_deg2_desq + cntrl_fem_deg2_desq 

fem_deg_counts_dsq <- data.frame(Category = c("Total", "DCM Female-Biased", "Control Female-Biased"),
                                 Count = c(fem_dcm_vs_cntrl_deg_dsq , dcm_fem_deg2_desq, cntrl_fem_deg2_desq))
fem_deg_counts_dsq$Category <- factor(fem_deg_counts_dsq$Category, levels = c("Total", "DCM Female-Biased", "Control Female-Biased"))


ggplot(fem_deg_counts_dsq, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7500)) +
  scale_y_continuous(breaks = seq(0, 7500, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between DCM Males and DCM Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "Control Female-Biased" = "#FF8C42", "DCM Female-Biased" = "#7B4173", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))


###########################################
# Analysis between hgps and control males #
###########################################
#Filter for significant results
dcm_male <- as.data.frame(subset(dcm_male_res, !is.na(padj)))
dcm_male$Gene <- rownames(dcm_male)

#checking the normality of the log2fold values
ks.test(rank(dcm_male$log2FoldChange, na.last = "keep"), "pnorm")
#using ranked z_scores because the log2fold change is nor normally distributed
dcm_male$rank_zscore <- qnorm(rank(dcm_male$log2FoldChange, na.last = "keep") / 
                                (sum(!is.na(dcm_male$log2FoldChange)) + 1))

hist(dcm_male$rank_zscore, breaks = 50, main = "Histogram of rank z-score Distribution of the DCM Males vs Control Males", xlab = "ranked z-score")

dcm_male_idx <- colData(dcm_dsq)$Group[c("DCM_Male", "Control_Male")]
dcm_male_idx <- colData(dcm_dsq)$Group %in% c("DCM_Male", "Control_Male")


#isolating the genes upregulated in dcm males
dcm_male_dsq_up2 <- subset(dcm_male, padj < 0.05 & log2FoldChange >= 0.58)
head(dcm_male_dsq_up2, 10)
dcm_male_dsq2_deg <- nrow(dcm_male_dsq_up2)

#isolating the genes upregulated in control males
cntrl_male_dsq_up2 <- subset(dcm_male, padj < 0.05 & rank_zscore < -0.58)
head(cntrl_male_dsq_up2)
cntrl_male_dsq2_deg <- nrow(cntrl_male_dsq_up2)

dcm_male_dsq_sig_deg <- dcm_male_dsq2_deg  + cntrl_male_dsq2_deg


male_deg_counts_dsq <- data.frame(Category = c("Total", "DCM Male-Biased", "Control Male-Biased"),
                                 Count = c(dcm_male_dsq_sig_deg, dcm_male_dsq2_deg, cntrl_male_dsq2_deg))
male_deg_counts_dsq$Category <- factor(male_deg_counts_dsq$Category, levels = c("Total", "DCM Male-Biased", "Control Male-Biased"))


ggplot(male_deg_counts_dsq, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7500)) +
  scale_y_continuous(breaks = seq(0, 7500, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between DCM Males and Control Males",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "DCM Male-Biased" = "#008B8B", "Control Male-Biased" = "#4A90E2", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

#################################################################
# Finding the genes that are uniquely expressed in DCM females  #
#################################################################
#getting unique dcm female genes
dcm_fem_specific_dsq1 <- dcm_fem_dsq_up1[!dcm_fem_dsq_up1$Gene %in% cntrl_fem_dsq_up1$Gene, ]
nrow(dcm_fem_specific_dsq1)

dcm_fem_specific_dsq2 <- dcm_fem_dsq_up2[!dcm_fem_dsq_up2$Gene %in% cntrl_fem_dsq_up1$Gene, ]
nrow(dcm_fem_specific_dsq2)
dcm_fem_specific_dsq2 <- dcm_fem_specific_dsq2[!dcm_fem_specific_dsq2$Gene %in% dcm_male_dsq_up2$Gene, ]
nrow(dcm_fem_specific_dsq2)

#organise by significance
dcm_fem_specific_dsq2<- dcm_fem_specific_dsq2[order(dcm_fem_specific_dsq2$padj, decreasing = FALSE), ]
head(dcm_fem_specific_dsq2, 10)

#################################################################
# Finding the genes that are uniquely expressed in DCM Males    #
#################################################################
#getting unique dcm male genes
dcm_male_specific_dsq1 <- dcm_male_dsq_up1[!dcm_male_dsq_up1$Gene %in% cntrl_male_dsq_up1_dsq$Gene, ]
nrow(dcm_male_specific_dsq1)

dcm_male_specific_dsq2 <- dcm_male_dsq_up2[!dcm_male_dsq_up2$Gene %in% cntrl_male_dsq_up1_dsq$Gene, ]
nrow(dcm_male_specific_dsq2)
dcm_male_specific_dsq2 <- dcm_male_specific_dsq2 [!dcm_male_specific_dsq2 $Gene %in% dcm_fem_dsq_up2$Gene, ]
nrow(dcm_male_specific_dsq2)

#organise by significance
dcm_male_specific_dsq2<- dcm_male_specific_dsq2[order(dcm_male_specific_dsq2$padj, decreasing = FALSE), ]
head(dcm_male_specific_dsq2, 10)

################################################################################
#                              1.B: edgeR Analysis                             #
################################################################################
dcm_group <- factor(dcm_comb_meta$Group)

dcm_ecounts <- DGEList(counts = dcm_data_comb, group = dcm_group)

#filtering genes based on count
dcm_ecounts_keep <- rowSums(dcm_ecounts$counts >= 10) >= 3
dcm_ecounts_cleaned <- dcm_ecounts[dcm_ecounts_keep, , keep.lib.sizes = FALSE]

#normalization
dcm_ecounts_cleaned <- calcNormFactors(dcm_ecounts_cleaned)

#creating the design matrix
dcm_investigation <- model.matrix(~0 + Group, data = dcm_comb_meta)
colnames(dcm_investigation) <- levels(dcm_group)

#estimating dispersion
dcm_ecounts_cleaned <- estimateDisp(dcm_ecounts_cleaned, dcm_investigation)

#fitting the model
dcm_fit <- glmFit(dcm_ecounts_cleaned, dcm_investigation)


# Contrasts for the analysis
dcm_contrasts <- makeContrasts(DCM_female_vs_control_female = DCM_Female - Control_Female,
                               DCM_male_vs_Control_male = DCM_Male - Control_Male,
                               DCM_female_vs_DCM_male = DCM_Female - DCM_Male,
                               Control_female_vs_Control_male = Control_Female - Control_Male,
                               levels = dcm_investigation)
#DCM female Control female
dcm_fem_v_c_fem <- glmLRT(dcm_fit, contrast = dcm_contrasts[, "DCM_female_vs_control_female"])
dcm_fem_v_c_fem_res <- topTags(dcm_fem_v_c_fem, n = Inf)$table
dcm_fem_v_c_fem_res <- as.data.frame(subset(dcm_fem_v_c_fem_res, !is.na(FDR)))
dcm_fem_v_c_fem_res$Gene <- rownames(dcm_fem_v_c_fem_res)

#checking distribution of log2foldchange using jarque bera
jarque.bera.test(dcm_fem_v_c_fem_res$logFC)
#using ranked z_scores because the log2fold change is not normally distributed
dcm_fem_v_c_fem_res$rank_zscore <- qnorm(rank(dcm_fem_v_c_fem_res$logFC, na.last = "keep") / 
                                              (sum(!is.na(dcm_fem_v_c_fem_res$logFC)) + 1))

#filtering for genes upregulated in dcm females
dcm_fem_e_up2 <- subset(dcm_fem_v_c_fem_res, FDR < 0.05 & logFC >= 0.58)
dcm_fem_deg2_e <- nrow(dcm_fem_e_up2)
#filtering for genes upregulated in control females
cntrl_fem_e_up2 <- subset(dcm_fem_v_c_fem_res, FDR < 0.05 & logFC < -0.58)
cntrl_fem_deg2_e <- nrow(cntrl_fem_e_up2)

fem_dcm_v_cntrl_deq_e <- dcm_fem_deg2_e + cntrl_fem_deg2_e

fem_deg_counts_e <- data.frame(Category = c("Total", "DCM Female-Biased", "Control Female-Biased"),
                                 Count = c(fem_dcm_v_cntrl_deq_e, dcm_fem_deg2_e, cntrl_fem_deg2_e))
fem_deg_counts_e$Category <- factor(fem_deg_counts_e$Category, levels = c("Total", "DCM Female-Biased", "Control Female-Biased"))


ggplot(fem_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7200)) +
  scale_y_continuous(breaks = seq(0, 7200, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between DCM Females and Control Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "Control Female-Biased" = "#FF8C42", "DCM Female-Biased" = "#7B4173", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))



#DCM male Control male
dcm_male_v_c_male <- glmLRT(dcm_fit, contrast = dcm_contrasts[, "DCM_male_vs_Control_male"])
dcm_male_v_c_male_res <- topTags(dcm_male_v_c_male, n = Inf)$table
dcm_male_v_c_male_res <- as.data.frame(subset(dcm_male_v_c_male_res, !is.na(FDR)))
dcm_male_v_c_male_res$Gene <- rownames(dcm_male_v_c_male_res)

#checking distribution of log2foldchange using jarque bera
jarque.bera.test(dcm_male_v_c_male_res$logFC)
#using ranked z_scores because the log2fold change is not normally distributed
dcm_male_v_c_male_res$rank_zscore <- qnorm(rank(dcm_male_v_c_male_res$logFC, na.last = "keep") / 
                                             (sum(!is.na(dcm_male_v_c_male_res$logFC)) + 1))


#filtering for genes upregulated in dcm males
dcm_male_e_up2 <- subset(dcm_male_v_c_male_res, FDR < 0.05 & logFC >= 0.58)
dcm_male_deg2_e <- nrow(dcm_male_e_up2) 
#filtering for genes upregulated in control males
cntrl_male_e_up2 <- subset(dcm_male_v_c_male_res, FDR < 0.05 & logFC < -0.58)
cntrl_male_deg2_e <- nrow(cntrl_male_e_up2)

male_dcm_v_cntrl_tot_deg_e <- dcm_male_deg2_e + cntrl_male_deg2_e

male_deg_counts_e <- data.frame(Category = c("Total", "DCM Male-Biased", "Control Male-Biased"),
                                  Count = c(male_dcm_v_cntrl_tot_deg_e, dcm_male_deg2_e, cntrl_male_deg2_e))
male_deg_counts_e$Category <- factor(male_deg_counts_e$Category, levels = c("Total", "DCM Male-Biased", "Control Male-Biased"))


ggplot(male_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7200)) +
  scale_y_continuous(breaks = seq(0, 7200, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between DCM Males and Control Males",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "DCM Male-Biased" = "#008B8B", "Control Male-Biased" = "#4A90E2", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

#DCM female DCM male
dcm_fem_v_dcm_male <- glmLRT(dcm_fit, contrast = dcm_contrasts[, "DCM_female_vs_DCM_male"])
dcm_fem_v_dcm_male_res <- topTags(dcm_fem_v_dcm_male, n = Inf)$table
dcm_fem_v_dcm_male_res <- as.data.frame(subset(dcm_fem_v_dcm_male_res, !is.na(FDR)))
dcm_fem_v_dcm_male_res$Gene <- rownames(dcm_fem_v_dcm_male_res)

#checking distribution of log2foldchange using jarque bera
jarque.bera.test(dcm_fem_v_dcm_male_res$logFC)
#using ranked z_scores because the log2fold change is not normally distributed
dcm_fem_v_dcm_male_res$rank_zscore <- qnorm(rank(dcm_fem_v_dcm_male_res$logFC, na.last = "keep") / 
                                              (sum(!is.na(dcm_fem_v_dcm_male_res$logFC)) + 1))

#filtering for genes upregulated in dcm females
dcm_fem_e_up1 <- subset(dcm_fem_v_dcm_male_res, FDR < 0.05 & logFC >= 0.58)
dcm_fem_deg1_e <- nrow(dcm_fem_e_up1)
#filtering for genes upregulated in control females
dcm_male_e_up1 <- subset(dcm_fem_v_dcm_male_res, FDR < 0.05 & logFC < -0.58)
dcm_male_deg1_e <- nrow(dcm_male_e_up1 ) 

dcm_male_vs_female_deg_e <- dcm_fem_deg1_e + dcm_male_deg1_e


dcm_deg_counts_e <- data.frame(Category = c("Total", "DCM Male-Biased", "DCM Female-Biased"),
                                 Count = c(dcm_male_vs_female_deg_e, dcm_male_deg1_e, dcm_fem_deg1_e))
dcm_deg_counts_e$Category <- factor(dcm_deg_counts_e$Category, levels = c("Total", "DCM Male-Biased", "DCM Female-Biased"))


ggplot(dcm_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7200)) +
  scale_y_continuous(breaks = seq(0, 7200, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between DCM Males and DCM Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "DCM Male-Biased" = "#008B8B", "DCM Female-Biased" = "#7B4173", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

# creating the matrix
n1 <- dcm_male_vs_female_deg_e;  x1 <- dcm_male_deg1_e
n2 <- dcm_male_vs_female_deg_e;  x2 <- dcm_fem_deg1_e     

tab <- matrix(c(x1, n1 - x1,
                x2, n2 - x2),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("Male", "Female"),
                              c("Significant", "Not_significant")))

#testing for significance
prop.test(tab, correct = FALSE)


#Control female Control male
cntrl_fem_v_cntrl_male <- glmLRT(dcm_fit, contrast = dcm_contrasts[, "Control_female_vs_Control_male"])
cntrl_fem_v_cntrl_male_res <- topTags(cntrl_fem_v_cntrl_male, n = Inf)$table
cntrl_fem_v_cntrl_male_res <- as.data.frame(subset(cntrl_fem_v_cntrl_male_res, !is.na(FDR)))
cntrl_fem_v_cntrl_male_res$Gene <- rownames(cntrl_fem_v_cntrl_male_res)

#checking distribution of log2foldchange using jarque bera
jarque.bera.test(cntrl_fem_v_cntrl_male_res$logFC)
#using ranked z_scores because the log2fold change is not normally distributed
cntrl_fem_v_cntrl_male_res$rank_zscore <- qnorm(rank(cntrl_fem_v_cntrl_male_res$logFC, na.last = "keep") / 
                                                  (sum(!is.na(cntrl_fem_v_cntrl_male_res$logFC)) + 1))

#filtering for genes upregulated in dcm females
cntrl_fem_e_up1 <- subset(cntrl_fem_v_cntrl_male_res, FDR < 0.05 & logFC >= 0.58)
cntrl_fem_deg1_e <- nrow(cntrl_fem_e_up1)

#filtering for genes upregulated in control females
cntrl_male_e_up1 <- subset(cntrl_fem_v_cntrl_male_res, FDR < 0.05 & logFC < -0.58)
cntrl_male_deg1_e <- nrow(cntrl_male_e_up1)

cntrl_male_vs_female_deg_e <- cntrl_fem_deg1_e + cntrl_male_deg1_e

cntrl_male_fem_deg_counts_e <- data.frame(Category = c("Total", "Control Male-Biased", "Control Female-Biased"),
                                       Count = c(cntrl_male_vs_female_deg_e, cntrl_fem_deg1_e, cntrl_male_deg1_e))
cntrl_male_fem_deg_counts_e$Category <- factor(cntrl_male_fem_deg_counts_e$Category, levels = c("Total", "Control Male-Biased", "Control Female-Biased"))


ggplot(cntrl_male_fem_deg_counts_e, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 7200)) +
  scale_y_continuous(breaks = seq(0, 7200, by = 300)) +
  labs(title = "Number of Differentially Expressed Genes Between Control Males and Control Females",
       x = "Category", y = "Number of DEGs") +
  scale_fill_manual(values = c("Total" = "black", "Control Male-Biased" = "#4A90E2", "Control Female-Biased" = "#FF8C42", "Shared" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30, hjust = 0.5, family = "sans"),
        axis.title.x = element_text(size = 27, family = "sans"),
        axis.title.y = element_text(size = 27, family = "sans"),
        axis.text.x = element_text(size = 21, family = "sans"),
        axis.text.y = element_text(size = 21, family = "sans"))

# creating the matrix
n1 <- cntrl_male_vs_female_deg_e;  x1 <- cntrl_fem_deg1_e
n2 <- cntrl_male_vs_female_deg_e;  x2 <- cntrl_male_deg1_e  

tab <- matrix(c(x1, n1 - x1,
                x2, n2 - x2),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("Control Female", "Control Male"),
                              c("Significant", "Not_significant")))

#testing for significance
prop.test(tab, correct = FALSE)



################################################################################
#Finding the genes exclusive to DCM females
dcm_fem_specific_e1 <-dcm_fem_e_up1[!dcm_fem_e_up1$Gene %in% cntrl_fem_e_up1$Gene, ]
nrow(dcm_fem_specific_e1)
dcm_fem_specific_e1 <- dcm_fem_specific_e1[order(dcm_fem_specific_e1$FDR, decreasing = FALSE), ]
head(dcm_fem_specific_e1, 10)

dcm_fem_specific_e2 <-dcm_fem_e_up2[!dcm_fem_e_up2$Gene %in% cntrl_fem_e_up1$Gene, ]
nrow(dcm_fem_specific_e2)
dcm_fem_specific_e2 <-dcm_fem_specific_e2[!dcm_fem_specific_e2$Gene %in% dcm_male_e_up2$Gene, ]
nrow(dcm_fem_specific_e2)
dcm_fem_specific_e2 <- dcm_fem_specific_e2[order(dcm_fem_specific_e2$FDR, decreasing = FALSE), ]
head(dcm_fem_specific_e2, 10)

#Finding the genes exclusive to DCM males
dcm_male_specific_e1 <- dcm_male_e_up1[!dcm_male_e_up1$Gene %in% cntrl_male_e_up1$Gene, ]
nrow(dcm_male_specific_e1)
dcm_male_specific_e1 <- dcm_male_specific_e1[order(dcm_male_specific_e1$FDR, decreasing = FALSE), ]
head(dcm_male_specific_e1, 10)

dcm_male_specific_e2 <- dcm_male_e_up2[!dcm_male_e_up2$Gene %in% cntrl_male_e_up1$Gene, ]
nrow(dcm_male_specific_e2)
dcm_male_specific_e2 <- dcm_male_specific_e2 [!dcm_male_specific_e2 $Gene %in% dcm_fem_e_up2$Gene, ]
nrow(dcm_male_specific_e2)
dcm_male_specific_e2 <- dcm_male_specific_e2[order(dcm_male_specific_e2$FDR, decreasing = FALSE), ]
head(dcm_male_specific_e2, 10)

################################################################################
# Combining the DESeq2 and edgeR data
################################################################################
#Getting the intersect of the uniquely expressed genes in HGPS females in the deseq2 and edgeR analyses
dcm_fem_specific_common1 <- inner_join(dcm_fem_specific_dsq1, dcm_fem_specific_e1, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(dcm_fem_specific_common1)

dcm_fem_specific_common2 <- inner_join(dcm_fem_specific_dsq2, dcm_fem_specific_e2, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(dcm_fem_specific_common2)

#reorganising the dataframe
dcm_fem_specific_common1 <- dcm_fem_specific_common1[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                           "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]
dcm_fem_specific_common2 <- dcm_fem_specific_common2[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                           "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]


dcm_fem_specific_common1 <- dcm_fem_specific_common1[!duplicated(dcm_fem_specific_common1$Gene), ]
rownames(dcm_fem_specific_common1) <- dcm_fem_specific_common1$Gene

dcm_fem_specific_common2 <- dcm_fem_specific_common2[!duplicated(dcm_fem_specific_common2$Gene), ]
rownames(dcm_fem_specific_common2) <- dcm_fem_specific_common2$Gene

#organise by padj
nrow(dcm_fem_specific_common1)
dcm_fem_specific_common1 <- dcm_fem_specific_common1[order(dcm_fem_specific_common1$padj, decreasing = FALSE), ]
head(dcm_fem_specific_common1, 10)

nrow(dcm_fem_specific_common2)
dcm_fem_specific_common2 <- dcm_fem_specific_common2[order(dcm_fem_specific_common2$padj, decreasing = FALSE), ]
head(dcm_fem_specific_common2, 10)


dcm_male_specific_common1 <- inner_join(dcm_male_specific_dsq1, dcm_male_specific_e1, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(dcm_male_specific_common1)

dcm_male_specific_common2 <- inner_join(dcm_male_specific_dsq2, dcm_male_specific_e2, by = "Gene", suffix = c(".DESeq2", ".edgeR"))
nrow(dcm_male_specific_common2)

#reorganising the dataframe
dcm_male_specific_common1 <- dcm_male_specific_common1[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                         "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]
dcm_fem_specific_common2 <- dcm_fem_specific_common2[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                                         "padj", "logFC", "logCPM", "LR", "PValue", "FDR", "rank_zscore.DESeq2", "rank_zscore.edgeR")]


dcm_male_specific_common1 <- dcm_male_specific_common1[!duplicated(dcm_male_specific_common1$Gene), ]
rownames(dcm_male_specific_common1) <- dcm_male_specific_common1$Gene

dcm_male_specific_common2 <- dcm_male_specific_common2[!duplicated(dcm_male_specific_common2$Gene), ]
rownames(dcm_male_specific_common2) <- dcm_male_specific_common2$Gene

#organise by padj
nrow(dcm_male_specific_common1)
dcm_male_specific_common1 <- dcm_male_specific_common1[order(dcm_male_specific_common1$padj, decreasing = FALSE), ]
head(dcm_male_specific_common1, 10)

nrow(dcm_male_specific_common2)
dcm_male_specific_common2 <- dcm_male_specific_common2[order(dcm_male_specific_common2$padj, decreasing = FALSE), ]
head(dcm_male_specific_common2, 10)


#################################################################################
# FGSEA
#

#DCM Males
dcm_males_genes <- inner_join(dcm_male_specific_common2, ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
dcm_males_genes <- dcm_males_genes[order(dcm_males_genes$padj, decreasing = FALSE), ]
nrow(hgps_males_genes1)

#clean the data set
#remove NAs
dcm_male_genes_cleaned <- dcm_males_genes[!is.na(dcm_males_genes$entrezgene_id), ]
nrow(dcm_male_genes_cleaned)

#collapse duplicates
dcm_male_genes_cleaned <- dcm_male_genes_cleaned[ !duplicated(dcm_male_genes_cleaned$entrezgene_id), ]
nrow(dcm_male_genes_cleaned)

#set the row names to entrez gene id
rownames(dcm_male_genes_cleaned) <- dcm_male_genes_cleaned$entrezgene_id

dcm_male_gene_names <- rownames(dcm_male_genes_cleaned)
dcm_male_gene_stats <- setNames(dcm_male_genes_cleaned$stat, dcm_male_genes_cleaned$entrezgene_id)
dcm_male_gene_stats <- dcm_male_gene_stats[is.finite(dcm_male_gene_stats)]

#creating the fgsea table
dcm_males_fgsea <- fgsea(pathways = gene_sets, 
                           stats = dcm_male_gene_stats, 
                           minSize = 5, 
                           maxSize = 500)

dcm_males_fgsea <- dcm_males_fgsea[order(dcm_males_fgsea$padj, decreasing = FALSE), ]
head(dcm_males_fgsea, 20)

#DCM Females
dcm_fem_genes <- inner_join(dcm_fem_specific_common2, ens2entrez, by = join_by("Gene"=="hgnc_symbol"))
dcm_fem_genes <- dcm_fem_genes[order(dcm_fem_genes$padj, decreasing = FALSE), ]
nrow(dcm_fem_genes)

#clean the data set
#remove NAs
dcm_fem_genes_cleaned <- dcm_fem_genes[!is.na(dcm_fem_genes$entrezgene_id), ]
nrow(dcm_fem_genes_cleaned)

#collapse duplicates
dcm_fem_genes_cleaned <- dcm_fem_genes_cleaned[ !duplicated(dcm_fem_genes_cleaned$entrezgene_id), ]
nrow(dcm_fem_genes_cleaned)

#set the row names to entrez gene id
rownames(dcm_fem_genes_cleaned) <- dcm_fem_genes_cleaned$entrezgene_id

dcm_fem_gene_names <- rownames(dcm_fem_genes_cleaned)
dcm_fem_gene_stats <- setNames(dcm_fem_genes_cleaned$stat, dcm_fem_genes_cleaned$entrezgene_id)
dcm_fem_gene_stats <- dcm_fem_gene_stats[is.finite(dcm_fem_gene_stats)]

#creating the fgsea table
dcm_fem_fgsea <- fgsea(pathways = gene_sets, 
                           stats = dcm_fem_gene_stats, 
                           minSize = 5, 
                           maxSize = 500)

nrow(dcm_fem_fgsea)
head(dcm_fem_fgsea, 20)
dcm_fem_fgsea <- dcm_fem_fgsea[order(dcm_fem_fgsea$padj, decreasing = FALSE), ]

################################################################################
#                                  2.D: enrichR                                #
################################################################################
enrichrdbs <- listEnrichrDbs()
colnames(enrichrdbs)
enrichrdbs <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025", "GO_Biological_Process_2025")

#investigating the unique hgps male genes
#first we need to only save the list of genes - we don't want stats
dcm_males_specific_common_list2 <- dcm_male_specific_common2$Gene

enriched_dcm_males2 <- enrichr(dcm_males_specific_common_list2, enrichrdbs)
head(enriched_dcm_males2)

e_dcm_males_processes2 <- enriched_dcm_males2[["GO_Biological_Process_2025"]]
e_dcm_males_processes2$Count <- as.numeric(sub("/.*", "", e_dcm_males_processes2$Overlap))
e_dcm_males_processes2 <- e_dcm_males_processes2[order(e_dcm_males_processes2$Adjusted.P.value, decreasing = FALSE), ]
View(e_dcm_males_processes2)

e_dcm_males_processes2 <- e_dcm_males_processes2 %>%
  separate(Term, into = c("Pathway", "GO_code"), sep = " \\(", remove = TRUE) %>%
  mutate(GO_code = sub("\\)", "", GO_code))
top_20_dcm_male_enrichr2 <- e_dcm_males_processes2[1:20, ]

#investigating the unique hgps female genes
#first we need to only save the list of genes - we don't want stats
dcm_females_specific_common_list2 <- dcm_fem_specific_common2$Gene

enriched_dcm_females2 <- enrichr(dcm_females_specific_common_list2, enrichrdbs)
head(enriched_dcm_females2)

head(enriched_dcm_females2[["GO_Biological_Process_2025"]])
e_dcm_females_processes2 <- enriched_dcm_females2[["GO_Biological_Process_2025"]]
e_dcm_females_processes2$Count <- as.numeric(sub("/.*", "", e_dcm_females_processes2$Overlap))
e_dcm_females_processes2 <- e_dcm_females_processes2[order(e_dcm_females_processes2$Adjusted.P.value, decreasing = FALSE), ]
e_dcm_females_processes2 <- e_dcm_females_processes2 %>%
  separate(Term, into = c("Pathway", "GO_code"), sep = " \\(", remove = TRUE) %>%
  mutate(GO_code = sub("\\)", "", GO_code))


################################################################################
#                                 2.E: gprofiler                                #
################################################################################
#as another enrichment pathway comparison we are using g profiler. 

#first investigating genes only expressed in females 1
dcm_fem_prof2 <- gost(dcm_fem_specific_common2$Gene, organism = "hsapiens", significant = TRUE, correction_method = "fdr")
names(dcm_fem_prof2)
dcm_fem_prof_res2 <- dcm_fem_prof2$result
dcm_fem_prof_res2 <- dcm_fem_prof_res2[order(dcm_fem_prof_res2$source, dcm_fem_prof_res2$p_value, decreasing = FALSE), ]

dcm_fem_prof_meta2 <- hgps_fem_prof2$meta

#gostplot(dcm_fem_prof_res2, capped = TRUE, interactive = TRUE)

dcm_fem_go2 <- dcm_fem_prof_res2 %>% filter(grepl("^GO:", term_id))

#first investigating genes only expressed in males
dcm_male_prof2 <- gost(hgps_male_specific_common2$Gene, organism = "hsapiens", significant = TRUE, correction_method = "fdr")
colnames(dcm_male_prof2$result)
dcm_male_prof_res2 <- dcm_male_prof2$result
dcm_male_prof_res2 <- dcm_male_prof_res2[order(dcm_male_prof_res2$source, dcm_male_prof_res2$p_value, decreasing = FALSE), ]

dcm_male_prof_meta2 <- dcm_male_prof2$meta

gostplot(dcm_male_prof2, capped = TRUE, interactive = TRUE)

dcm_male_go2 <- dcm_male_prof_res2 %>% filter(grepl("^GO:", term_id))


###
# comparing pathways
#disease effect
dcm_fem_com_path2 <- inner_join(dcm_fem_go2, e_dcm_females_processes2, by = c("term_id" = "GO_code"))
dcm_fem_com_path2 <- dcm_fem_com_path2[order(dcm_fem_com_path2$Adjusted.P.value, decreasing = FALSE), ]
dcm_fem_com_path2 <- dcm_fem_com_path2[!duplicated(dcm_fem_com_path2$term_name), ]
nrow(dcm_fem_com_path2)

dcm_fem_com_path2_plot <- dcm_fem_com_path2 %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(Adjusted.P.value),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(dcm_fem_com_path2_plot, aes(x = term_name, y = log10_padj, fill = log10_padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#FF8C42", "#B35E00"))(2),  # 20-step gradient
    name = expression(-log[10]("padj"))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title  = element_text(size = 12),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12))

#hpgs males (sex difference
dcm_male_com_path2 <- inner_join(dcm_male_go2, e_dcm_males_processes2, by = c("term_id" = "GO_code"))
dcm_male_com_path2 <- dcm_male_com_path2[order(dcm_male_com_path2$Adjusted.P.value, decreasing = FALSE), ]
dcm_male_com_path2 <- dcm_male_com_path2[!duplicated(dcm_male_com_path2$term_name), ]
nrow(dcm_male_com_path2)


dcm_male_com_path2_plot <- dcm_male_com_path2 %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20) %>%
  mutate(log10_padj = -log10(Adjusted.P.value),                
         term_name = factor(term_name, levels = rev(term_name)))

ggplot(dcm_male_com_path2_plot, aes(x = term_name, y = log10_padj, fill = log10_padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#030e1e", "#1F497D"))(2),  # 20-step gradient
    name = expression(-log[10]("padj"))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10]("padj")),
    title = "Top Enriched Pathways"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, family = "Serif"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

################################################################################
# Gene Network Analysis
################################################################################
GSE113957_dseq_cpm #matrix

hgps_male_specific_common1 #stats
hgps_fem_specific_common1 #stats
hgps_male_specific_common2 #stats
hgps_fem_specific_common2 #stats

lmna_col_hgps <- hgps_common["LMNA", ]   
hgps_male_specific_lmna1 <- rbind(hgps_male_specific_common1, lmna_col_hgps)
hgps_fem_specific_lmna1 <- rbind(hgps_fem_specific_common1, lmna_col_hgps)

PearsCorr <- sapply(setdiff(colnames(), "LMNA"),
                    function(x) cor(tra_df[, x], tra_df$LMNA, method = "pearson"))

lmna_col_fem <- female_common["LMNA", ] 
hgps_fem_specific_lmna2 <- rbind(hgps_fem_specific_common2, lmna_col_hgps)

lmna_col_male <- male_common["LMNA", ] 
hgps_male_specific_lmna2 <- rbind(hgps_male_specific_common2, lmna_col_hgps)



################################################################################
#                   2.C: Combining the dcm and hgps analysis                  #
################################################################################
## what genes are prevalent in dcm females and hgps females? from the affected vs control females analysis. 
lamin_fem <- intersect(rownames(hgps_fem_specific_common2), rownames(dcm_fem_specific_common2))
length(lamin_fem)
head(lamin_fem)

lamin_fem_hgps <- hgps_fem_specific_common2[c(lamin_fem), , drop = FALSE ]
lamin_fem_hgps_samples <- colnames(GSE113957_dseq_cpm)[GSE113957_female_hgps_vs_cntrl_idx]
lamin_fem_hgps_annot <- data.frame(Group = GSE113957_meta$Group[match(hgps_male_fem_dseq_samples, rownames(GSE113957_meta))],
                                       row.names = hgps_male_fem_dseq_samples)
lamin_fem_hgps_annot_cols <- list(Group = c("HGPS_female" = "#B35E00", "Control_female" = "#1F497D"))

lamin_fem_hgps_names <- rownames(lamin_fem_hgps)

lamin_fem_hgps_annot <- data.frame(Group = GSE113957_meta$Group[match(lamin_fem_hgps_samples, rownames(GSE113957_meta))],
                                              row.names = lamin_fem_hgps_samples)
pheatmap(GSE113957_dseq_cpm[lamin_fem_hgps_names, GSE113957_female_hgps_vs_cntrl_idx],
         scale="row",
         show_rownames=T,
         main="Top 20 Genes Upregulated in HGPS Females Relative to Control Females",
         fontfamily = "Arial", 
         fontsize = 20,
         annotation_col= fem_hgps_cntrl_dseq_annot,
         annotation_colors = fem_hgps_cntrl_dseq_annot_cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE)

lamin_fem_dcm <- dcm_fem_specific_common2[c(lamin_fem), , drop = FALSE ]
  
lamin_male <- intersect(rownames(hgps_male_specific_common2), rownames(dcm_male_specific_common2))
length(lamin_male)
lamin_male

lamin_male_hgps <- hgps_male_specific_common2[c(lamin_male), , drop = FALSE ]
lamin_male_dcm <- dcm_male_specific_common2[c(lamin_male), , drop = FALSE ]



lmna_fem_com_pathways2 <- intersect(hgps_fem_com_path2$Pathway, dcm_fem_com_path2$Pathway)
lmna_male_com_pathways2 <- intersect(hgps_male_com_path2$Pathway, dcm_male_com_path2$Pathway)

install.packages("circlize")   # once
library(circlize)

female_paths <- lmna_fem_com_pathways2

edges_f <- rbind(
  data.frame(from = "HGPS_female", to = female_paths, weight = 1, stringsAsFactors = FALSE),
  data.frame(from = "DCM_female",  to = female_paths, weight = 1, stringsAsFactors = FALSE)
)

circos.clear()

# Sector order: put the set nodes at opposite sides with pathway nodes between
order_f <- c("HGPS_female", female_paths, "DCM_female")

# Optional: bigger gaps around the two set nodes to make them stand out
gap_vec <- c(12, rep(1, length(female_paths)), 12)

circos.par(gap.after = gap_vec, start.degree = 90)  # rotate start if you like

# Pre-allocate a track to add labels
chordDiagram(
  x = edges_f[, c("from", "to", "weight")],
  order = order_f,
  grid.col = c(
    "HGPS_female" = "#B35E00",   # pick your colours
    "DCM_female"  = "#7B4173",
    setNames(rep("#97503A", length(female_paths)), female_paths)
  ),
  transparency = 0.6,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08)
)

# Add readable sector labels (wrap long pathway names if needed)
circos.trackPlotRegion(
  track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim), y = ylim[1] + 0.1,
      labels = sector.name, facing = "clockwise",
      niceFacing = TRUE, adj = c(0, 0.5),
      cex = ifelse(sector.name %in% c("HGPS_female","DCM_female"), 0.9, 0.55)
    )
  },
  bg.border = NA
)

male_paths <- lmna_male_com_pathways2

edges_m <- rbind(
  data.frame(from = "HGPS_male", to = male_paths, weight = 1, stringsAsFactors = FALSE),
  data.frame(from = "DCM_male",  to = male_paths, weight = 1, stringsAsFactors = FALSE)
)

circos.clear()
order_m <- c("HGPS_male", male_paths, "DCM_male")
gap_vec_m <- c(12, rep(1, length(male_paths)), 12)
circos.par(gap.after = gap_vec_m, start.degree = 90)

chordDiagram(
  x = edges_m[, c("from", "to", "weight")],
  order = order_m,
  grid.col = c(
    "HGPS_male" = "#1F497D",
    "DCM_male"  = "#008B8B",
    setNames(rep("#0243DF", length(male_paths)), male_paths)
  ),
  transparency = 0.6,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08)
)

circos.trackPlotRegion(
  track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim), y = ylim[1] + 0.1,
      labels = sector.name, facing = "clockwise",
      niceFacing = TRUE, adj = c(0, 0.5),
      cex = ifelse(sector.name %in% c("HGPS_male","DCM_male"), 0.9, 0.55)
    )
  },
  bg.border = NA
)


