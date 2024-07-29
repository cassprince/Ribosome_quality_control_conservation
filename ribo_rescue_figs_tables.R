library(tidyverse)
library(rhmmer)
library(treeio)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(castor)
library(ggbreak)
library(ggprism)

table_for_tree = function(phylum, gene, accessions){
  table = as.data.frame(prop.table(table(phylum, gene), margin = 1)*100)
  table_new = data.frame(cbind(gene = table$Freq[table$gene == "TRUE"]))
  rownames(table_new) = accessions
  return(table_new)
}

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina")

df_GCF_nc_compressed = read.csv("ribo_rescue_df.csv") 

df_GCF_nc = separate_rows(df_GCF_nc_compressed, nuccore, sep = ";")

# Collapsed tree
tree = read.newick("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//16S_fasttree.tre")
tree = get_subtree_with_tips(tree, only_tips = df_GCF_nc$nuccore)$subtree
tree_tip_GCF = left_join(data.frame(tree$tip.label), df_GCF_nc, by = join_by("tree.tip.label" == "nuccore"), multiple = "any")
tree$tip.label = tree_tip_GCF$assembly

dfShort = data.frame(df_GCF_nc_compressed[df_GCF_nc_compressed$phylum %in% names(which(table(df_GCF_nc_compressed$phylum)>10)),])
n_vals = data.frame(table(dfShort$phylum))

random = dfShort %>%
  group_by(phylum) %>%
  sample_n(1) %>%
  select(-X) %>%
  ungroup()

rownames(random) = random$assembly
#If the row names and the first column don't match, the tip labels don't get written... So weird.
#write.csv(random, "random_genomes_for_tree_new.csv")

random = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//random_genomes_for_tree_new.csv")
random = column_to_rownames(random, 'X')

subtree_bar = get_subtree_with_tips(tree, only_tips = rownames(random))$subtree
tree_mid = midpoint(subtree_bar)

table_smpB = table_for_tree(dfShort$phylum, dfShort$smpB, random$assembly)
table_ssrA = table_for_tree(dfShort$phylum, dfShort$ssrA, random$assembly)
table_arfA = table_for_tree(dfShort$phylum, dfShort$arfA, random$assembly)
table_arfB = table_for_tree(dfShort$phylum, dfShort$arfB, random$assembly)
table_rqcH = table_for_tree(dfShort$phylum, dfShort$rqcH, random$assembly)
table_mutS2 = table_for_tree(dfShort$phylum, dfShort$mutS2, random$assembly)
table_smrB = table_for_tree(dfShort$phylum, dfShort$smrB, random$assembly)

genes = data.frame(cbind(smpB = table_smpB$gene, ssrA = table_ssrA$gene, arfB = table_arfB$gene, arfA = table_arfA$gene, smrB = table_smrB$gene, rqcH = table_rqcH$gene, mutS2 = table_mutS2$gene))
rownames(genes) = rownames(random)



p = ggtree(tree_mid, size = 0.8) %<+% random + 
  xlim(0, 11.75) + 
  geom_tiplab(aes(label=phylum), align = TRUE, size = 4) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") + 
  new_scale_fill()

p1 = gheatmap(p, genes, offset = 6.7, width=0.7, font.size=3.5, colnames = TRUE, color=NA, colnames_angle = 90, colnames_offset_y = -0.2) + 
  scale_fill_viridis_c(option="A", direction = -1, name="Percent\nwith gene")

ggsave("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//Figures//ribo_rescue_collapsed_tree_7_24_24.png", p1, width = 8.5, height = 6.4, dpi = 600, units = "in")

####Stats

rrf_table = data.frame(cbind(smpB = table(df_GCF_nc_compressed$smpB), ssrA = table(df_GCF_nc_compressed$ssrA), arfA = table(df_GCF_nc_compressed$arfA), arfB = table(df_GCF_nc_compressed$arfB), rqcH = table(df_GCF_nc_compressed$rqcH), mutS2 = table(df_GCF_nc_compressed$mutS2), smrB = table(df_GCF_nc_compressed$smrB)))

rrf_prop_table = data.frame(cbind(smpB = prop.table(table(df_GCF_nc_compressed$smpB)), ssrA = prop.table(table(df_GCF_nc_compressed$ssrA)), arfA = prop.table(table(df_GCF_nc_compressed$arfA)), arfB = prop.table(table(df_GCF_nc_compressed$arfB)), rqcH = prop.table(table(df_GCF_nc_compressed$rqcH)), mutS2 = prop.table(table(df_GCF_nc_compressed$mutS2)), smrB = prop.table(table(df_GCF_nc_compressed$smrB))))

format(round(rrf_prop_table*100,1),nsmall=1)
format(rrf_table,nsmall=1)
df_GCF_nc_compressed %>%
  summarise(sum = sum(ssrA == TRUE|smpB == TRUE), n = n(), perc = sum/n)

#rqcH and mutS2
rqcH = factor(c("without", "with", "without","with"))
mutS2 = factor(c("without", "without", "with", "with"))
vals = data.frame(table(df_GCF_nc_compressed$rqcH, df_GCF_nc_compressed$mutS2))$Freq
conf = data.frame(rqcH, mutS2, vals)

plot = ggplot(data =  conf, mapping = aes(x = rqcH, y = mutS2)) +
  geom_tile(aes(fill = vals), colour = "black") +
  geom_text(aes(label = sprintf("%1.0f", vals)), vjust = 1, size = 5) + 
  scale_fill_gradient(low = "white", high = "#F3AC71") +
  theme_bw() + theme(legend.position = "none") +
  theme(text = element_text(size = 20)) +
  theme(axis.title=element_text(face="italic"))


ggsave("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//Figures//rqcH_mutS2_confmatrix_7_26_24.png", plot, width = 5, height = 3.5, dpi = 600, units = "in")


#ssrA and smpB
ssrA = factor(c("without", "with", "without","with"))
smpB = factor(c("without", "without", "with", "with"))
vals = data.frame(table(df_GCF_nc_compressed$ssrA, df_GCF_nc_compressed$smpB))$Freq
conf = data.frame(ssrA, smpB, vals)

plot = ggplot(data =  conf, mapping = aes(x = ssrA, y = smpB)) +
  geom_tile(aes(fill = vals), colour = "black") +
  geom_text(aes(label = sprintf("%1.0f", vals)), vjust = 1, size = 5) +
  scale_fill_gradient(low = "white", high = "#F3AC71") +
  theme_bw() + theme(legend.position = "none") +
  theme(text = element_text(size = 20)) +
  theme(axis.title=element_text(face="italic"))

ggsave("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//Figures//ssrA_smrB_confmatrix_7_26_24.png", plot, width = 5, height = 3.5, dpi = 600, units = "in")

# Tables
random_IDs = random %>%
  select(phylum, assembly)

sup_table = df_GCF_nc_compressed %>%
  group_by(phylum) %>%
  summarise(total_genomes = n(), smpB = sum(smpB == TRUE), ssrA = sum(ssrA == TRUE), arfB = sum(arfB), arfA = sum(arfA), smrB = sum(smrB), rqcH = sum(rqcH), mutS2 = sum(mutS2)) 


sup_table_final = sup_table %>%
  left_join(random_IDs, by = join_by(phylum)) %>%
  rename(accession_for_tree = assembly)


write_csv(sup_table_final, "C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//ribo_rescue_gene_counts.csv")

# Supp table 1

df_lineage_string = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//ref_lineage.txt", col.names = "organism.taxId", header = FALSE) %>% 
  separate(organism.taxId, into = c("organism.taxId", "taxonomy"), sep = "\\s", extra = "merge")
df_lineage_string$organism.taxId = as.numeric(df_lineage_string$organism.taxId) 

df_lineage_string = df_lineage_string %>%
  drop_na()

sup_table_2 = inner_join(df_GCF_nc_compressed, df_lineage_string, unmatched = "drop", multiple = "any", by = "organism.taxId") %>% 
  select(assembly, nuccore, smpB, ssrA, arfB, arfA, smrB, rqcH, mutS2, organism.taxId, checkmInfo.completeness, checkmInfo.contamination, taxonomy)


write.csv(sup_table_2, file = "C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//Table_S1_spp_and_genes_7_26_24.csv")


#### Big tree
dfShort = dfShort %>%
  distinct(assembly, .keep_all = TRUE)

smpB = data.frame(dfShort$smpB)
rownames(smpB) = dfShort$assembly
ssrA = data.frame(dfShort$ssrA)
rownames(smpB) = dfShort$assembly
phylum = data.frame(dfShort$phylum)
rownames(phylum) = dfShort$assembly

df_big_tree = dfShort %>% select(smpB, ssrA, phylum) 
rownames(df_big_tree) = dfShort$assembly

subtree = get_subtree_with_tips(tree, only_tips = rownames(df_big_tree))$subtree
big_tree_mid = midpoint(subtree)

p = ggtree(big_tree_mid, layout='circular', size=0.2)

p1 = gheatmap(p, df_big_tree, offset = 0, width=0.2, font.size=1.5, colnames = FALSE, color=NA) + 
  scale_fill_manual(values=c("TRUE" = "white", "FALSE" = "black"), na.value = "white") + 
  new_scale_fill()

p2 = gheatmap(p1, phylum, offset = 0.8, width=0.05, font.size=1, colnames = FALSE, color=NA) 



ggsave(".//Figures//smpB_ssrA_tree.png", p2, width = 15, height = 12, dpi = 600, units = "in")
