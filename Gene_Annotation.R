getwd()
setwd("C:/Users/amani/OneDrive/Desktop/Intern_IITM/Project/Files")

if("org.Hs.eg.db" %in% rownames(installed.packages())) {
  print("Package is installed")
} else {
  print("Package is not installed")
}
BiocManager::install("clusterProfiler")
BiocManager::install("GO.db")
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)

df <- read.csv("C:/Users/amani/OneDrive/Desktop/Intern_IITM/Project/Files/gene_1.csv", stringsAsFactors = FALSE)
l1 <- as.list(df[c(2)])  

# Step 1: Flatten the list into a character vector
gene_char <- as.character(unlist(l1))

# Step 2: Split on '///' and trim whitespace
split_genes <- strsplit(gene_char, "\\s*///\\s*|\\s*-\\s*|\\s*\\|\\s*|\\s*_\\s*")
split_genes <- lapply(split_genes, trimws)

# Step 3: Flatten and get unique genes
clean_genes <- unique(unlist(split_genes))

#step 4: remove microRNA
# Remove miRNAs (optional)
clean_genes <- clean_genes[!grepl("^MIR|^miR", clean_genes)]

# Remove empty strings and NAs
clean_genes <- clean_genes[clean_genes != "" & !is.na(clean_genes)]

#validation-> valid genes are there 
valid_genes <- bitr(clean_genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

#Group Genes by GO Biological Process
go_groups <- groupGO(gene     = valid_genes$ENTREZID,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "BP",       # Biological Process
                     level    =2,          # as level increases the depth of info increases
                     readable = TRUE)
#exploring results
(go_groups@result[, c("Description", "Count")])

#it tells you which biological processes are statistically enriched (overrepresented) in your gene list.

ego <- enrichGO(gene         = valid_genes$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",        # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable     = TRUE)


head(ego@result[, c("Description", "geneID", "p.adjust")])

#Metabolic Genes Extraction

metabolic_terms <- ego@result[grepl("metaboli", ego@result$Description, ignore.case = TRUE), ]

metabolic_genes <- unique(unlist(strsplit(metabolic_terms$geneID, "/")))

library(dplyr)
library(tidyr)

# Step 1: Create a vector of metabolic keywords (you already did this)
metabolic_keywords <- c(
  "metabolic process", "biosynthetic process", "catabolic process", "primary metabolic process",
  "glycolysis", "gluconeogenesis", "TCA cycle", "citric acid cycle", "oxidative phosphorylation",
  "fatty acid beta-oxidation", "amino acid metabolism", "lipid metabolism", "nucleotide metabolism",
  "carbohydrate metabolism", "cholesterol biosynthetic process", "steroid biosynthetic process",
  "pyruvate metabolism", "pentose phosphate pathway", "urea cycle", "beta-oxidation",
  "glyoxylate cycle", "purine metabolism", "pyrimidine metabolism",
  "heme biosynthetic process", "acetyl-CoA metabolism", "ketone body metabolism"
)

# Step 2: Collapse keywords into a regex pattern
metabolic_pattern <- paste(metabolic_keywords, collapse = "|")

# Step 3: Filter GO terms related to metabolism
metabolic_terms_df <- ego@result %>%
  filter(grepl(metabolic_pattern, Description, ignore.case = TRUE))

# Step 4: Expand genes into one per row
metabolic_gene_df <- metabolic_terms_df %>%
  select(Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(GeneSymbol = geneID, Function = Description) %>%
  distinct()

# Optional: sort by gene
metabolic_gene_df <- metabolic_gene_df %>% arrange(GeneSymbol)

# View the final table
head(metabolic_gene_df)
metabolic_pattern <- paste(metabolic_keywords, collapse = "|")
metabolic_terms_2 <- ego@result[grepl(metabolic_pattern, ego@result$Description, ignore.case = TRUE), ]
metabolic_genes_1 <- unique(unlist(strsplit(metabolic_terms_2$geneID, "/")))

write.csv(metabolic_gene_df, "metabolic_gene_list_function.csv", row.names = FALSE)

write.csv(metabolic_genes, "metabolic_gene_list.csv", row.names = FALSE)

# For regulatory genes extraction

regulatory_keywords <- c(
  # General regulation
  "regulation of transcription",
  "positive regulation of transcription",
  "negative regulation of transcription",
  "regulation of gene expression",
  "epigenetic regulation",
  "regula",
  # Cell cycle and signal regulation
  "regulation of cell cycle",
  "regulation of cell proliferation",
  "regulation of kinase activity",
  "regulation of protein phosphorylation",
  "regulation of apoptotic process",
  
  # Transcriptional/post-transcriptional control
  "transcription factor activity",
  "RNA polymerase II transcription regulator activity",
  "chromatin remodeling",
  "histone modification",
  "DNA methylation"
)

# Regulatory
regulatory_pattern <- paste(regulatory_keywords, collapse = "|")
regulatory_terms <- ego@result[grepl(regulatory_pattern, ego@result$Description, ignore.case = TRUE), ]
regulatory_genes <- unique(unlist(strsplit(regulatory_terms$geneID, "/")))

write.csv(regulatory_genes, "regulatory_gene_list.csv", row.names = FALSE)

#Immune Response
immune_keywords <- c(
  "immune response",
  "immun",
  "inflammatory response",
  "adaptive immune response",
  "innate immune response",
  "T cell activation",
  "B cell activation",
  "antigen processing",
  "cytokine production",
  "leukocyte activation",
  "macrophage activation",
  "complement activation",
  "interferon signaling"
)
immune_pattern <- paste(immune_keywords, collapse = "|")
immune_terms <- ego@result[grepl(immune_pattern, ego@result$Description, ignore.case = TRUE), ]
immune_genes <- unique(unlist(strsplit(immune_terms$geneID, "/")))

write.csv(immune_genes, "immune_gene_list.csv", row.names = FALSE)

#apoptosis
apop_keywords <- c(
  "apoptotic process", "programmed cell death", "cell death", "regulation of apoptosis",
  "negative regulation of apoptosis", "positive regulation of apoptosis", "caspase activation",
  "intrinsic apoptotic signaling pathway", "extrinsic apoptotic signaling pathway","apop"
)
apop_pattern <- paste(apop_keywords, collapse = "|")
apop_terms <- ego@result[grepl(apop_pattern, ego@result$Description, ignore.case = TRUE), ]
apop_genes <- unique(unlist(strsplit(apop_terms$geneID, "/")))

write.csv(apop_genes, "apop_gene_list.csv", row.names = FALSE)

#signalling
SignalTransduction <- c(
  "signal transduction", "cell communication", "MAPK cascade", "Wnt signaling pathway",
  "Notch signaling pathway", "Hedgehog signaling pathway", "calcium signaling pathway",
  "receptor signaling", "second messenger mediated signaling", "cAMP signaling", "JAK-STAT cascade"
)
signal_pattern <- paste(SignalTransduction, collapse = "|")
signal_terms <- ego@result[grepl(signal_pattern, ego@result$Description, ignore.case = TRUE), ]
signal_genes <- unique(unlist(strsplit(signal_terms$geneID, "/")))

write.csv(signal_genes, "signal_gene_list.csv", row.names = FALSE)

#transport
Transport <- c(
  "transport", "transmembrane transport", "vesicle-mediated transport", "protein transport",
  "ion transport", "lipid transport", "glucose transport", "ABC transport",
  "endocytosis", "exocytosis"
)
transport_pattern <- paste(Transport, collapse = "|")
transport_terms <- ego@result[grepl(transport_pattern, ego@result$Description, ignore.case = TRUE), ]
transport_genes <- unique(unlist(strsplit(transport_terms$geneID, "/")))

write.csv(transport_genes, "transport_gene_list.csv", row.names = FALSE)

Development <- c(
  "developmental process", "embryonic development", "organ development", "cell differentiation",
  "neuron differentiation", "axon guidance", "morphogenesis", "organ morphogenesis",
  "pattern specification process", "limb development","developement"
)
dev_pattern <- paste(Development, collapse = "|")
dev_terms <- ego@result[grepl(dev_pattern, ego@result$Description, ignore.case = TRUE), ]
dev_genes <- unique(unlist(strsplit(dev_terms$geneID, "/")))

write.csv(dev_genes, "dev_gene_list.csv", row.names = FALSE)


# Load required libraries
library(dplyr)
library(tidyr)

# STEP 1: Define keyword dictionaries for each category
category_keywords <- list(
  Metabolic = c(
    "metabolic process", "biosynthetic process", "catabolic process", "primary metabolic process",
    "glycolysis", "gluconeogenesis", "TCA cycle", "citric acid cycle", "oxidative phosphorylation",
    "fatty acid beta-oxidation", "amino acid metabolism", "lipid metabolism", "nucleotide metabolism",
    "carbohydrate metabolism", "cholesterol biosynthetic process", "steroid biosynthetic process",
    "pyruvate metabolism", "pentose phosphate pathway", "urea cycle", "beta-oxidation",
    "glyoxylate cycle", "purine metabolism", "pyrimidine metabolism",
    "heme biosynthetic process", "acetyl-CoA metabolism", "ketone body metabolism","metaboli"
  ),
  
  Immune = c(
    "immune response", "inflammatory response", "adaptive immune response", "innate immune response",
    "T cell activation", "B cell activation", "antigen processing", "cytokine production",
    "leukocyte activation", "macrophage activation", "complement activation", "interferon signaling"
  ),
  
  Apoptosis = c(
    "apoptotic process", "programmed cell death", "cell death", "regulation of apoptosis",
    "negative regulation of apoptosis", "positive regulation of apoptosis", "caspase activation",
    "intrinsic apoptotic signaling pathway", "extrinsic apoptotic signaling pathway"
  ),
  
  SignalTransduction = c(
    "signal transduction", "cell communication", "MAPK cascade", "Wnt signaling pathway",
    "Notch signaling pathway", "Hedgehog signaling pathway", "calcium signaling pathway",
    "receptor signaling", "second messenger mediated signaling", "cAMP signaling", "JAK-STAT cascade"
  ),
  
  Transport = c(
    "transport", "transmembrane transport", "vesicle-mediated transport", "protein transport",
    "ion transport", "lipid transport", "glucose transport", "ABC transport",
    "endocytosis", "exocytosis"
  ),
  
  Development = c(
    "developmental process", "embryonic development", "organ development", "cell differentiation",
    "neuron differentiation", "axon guidance", "morphogenesis", "organ morphogenesis",
    "pattern specification process", "limb development"
  ),
  
  regulatory =c(
    # General regulation
    "regulation of transcription",
    "positive regulation of transcription",
    "negative regulation of transcription",
    "regulation of gene expression",
    "epigenetic regulation",
    "regula",
    # Cell cycle and signal regulation
    "regulation of cell cycle",
    "regulation of cell proliferation",
    "regulation of kinase activity",
    "regulation of protein phosphorylation",
    "regulation of apoptotic process",
    
    # Transcriptional/post-transcriptional control
    "transcription factor activity",
    "RNA polymerase II transcription regulator activity",
    "chromatin remodeling",
    "histone modification",
    "DNA methylation"
  )
)

# STEP 2: Create a combined data frame
combined_df <- data.frame()

for (category in names(category_keywords)) {
  pattern <- paste(category_keywords[[category]], collapse = "|")
  
  # Filter GO terms by keyword match
  cat_terms_df <- ego@result %>%
    filter(grepl(pattern, Description, ignore.case = TRUE)) %>%
    select(Description, geneID) %>%
    mutate(Functional_Category = category)
  
  # Expand gene list into separate rows
  cat_gene_df <- cat_terms_df %>%
    separate_rows(geneID, sep = "/") %>%
    rename(GeneSymbol = geneID, Function = Description)
  
  combined_df <- bind_rows(combined_df, cat_gene_df)
}

# STEP 3: Clean and remove duplicates
classified_genes <- combined_df %>%
  distinct(GeneSymbol, Function, Functional_Category) %>%
  arrange(GeneSymbol)

# STEP 4: Save to CSV (optional)
write.csv(classified_genes, "Classified_Genes_Functional_Groups.csv", row.names = FALSE)

# STEP 5: View first few rows
head(classified_genes)


#venn diagram
library(dplyr)
library(tibble)
library(VennDiagram)

# Filter to top 5 categories (example)
top_categories <- c("Metabolic", "Immune", "Apoptosis", "SignalTransduction", "Transport")

gene_sets_sub <- classified_genes %>%
  filter(Functional_Category %in% top_categories) %>%
  group_by(Functional_Category) %>%
  summarise(Genes = list(unique(GeneSymbol))) %>%
  deframe()

fill_colors <- c("skyblue", "orange", "mediumseagreen", "orchid", "plum")

venn.plot <- venn.diagram(
  x = gene_sets_sub,
  category.names = names(gene_sets_sub),
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 500,
  fill = fill_colors,
  cat.cex = 0.8,
  cex = 0.7,
  fontfamily = "sans",
  main = "Gene Overlaps Between Functional Categories (Top 5)"
)

png("Gene_Overlap_Venn_Top5.png", width = 3000, height = 3000, res = 500)
grid.draw(venn.plot)
dev.off()

library(dplyr)
library(tidyr)

# Assume `classified_genes` has columns: GeneSymbol, Function, Functional_Category

# Step 1: Remove duplicates for gene-category pairs
unique_gene_category <- classified_genes %>%
  select(GeneSymbol, Functional_Category) %>%
  distinct() %>%              # Keep unique gene-category pairs only
  arrange(GeneSymbol, Functional_Category)

# Step 2: View or save the cleaned dataframe
head(unique_gene_category)

# Optional: Save to CSV
write.csv(unique_gene_category, "Genes_Unique_Gene_Category.csv", row.names = FALSE)






