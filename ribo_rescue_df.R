library(tidyverse)
library(rhmmer)
library(treeio)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(castor)
library(ggbreak)
library(ggprism)
library(jsonlite)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//HMMER_6_6_24")

read_nhmmer = function(file){
  df = read.table(file, skip = 2, sep = "", colClasses = c("character", rep("NULL", 3), rep("character", 11), rep("NULL", 8)), fill = TRUE, row.names = NULL)
  df = df[,1:12]
  df = replace(df, df=='', NA)
  colnames(df) = c("accessions", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "strand", "eval", "score", "bias")
  df = drop_na(df)
  df = df %>% mutate_at(c("hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "eval", "score", "bias"), as.numeric)
  
  return(df)
}

# Upload lineages
lineages = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//ref_lineage.txt", col.names = "taxID", header = FALSE) %>% 
                      separate(taxID, into = c("organism", "domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "merge") %>% 
                      separate(organism, into = c("taxID", "organism"), sep = "\\s", extra = "merge") %>%
  drop_na()

# Upload GCF assembly accessions w corresponding NC nuccore accessions.
df_GCF_nc = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//GCF_nuccore_reps_clean.csv")

lines = readLines("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//prfB//Data//Bioinformatics//assembly_data_report.jsonl")
lines = lapply(lines, fromJSON)
lines = lapply(lines, unlist)
df_GCF_tax = bind_rows(lines)

df_GCF_tax_select = df_GCF_tax %>%
  select(accession, assemblyInfo.assemblyLevel, organism.organismName, assemblyInfo.bioprojectLineage.bioprojects.title, organism.taxId, checkmInfo.checkmMarkerSet, checkmInfo.completeness, checkmInfo.contamination, assemblyStats.totalSequenceLength, assemblyStats.gcPercent, assemblyStats.numberOfComponentSequences)

df_GCF_nc_tax = left_join(df_GCF_nc, df_GCF_tax_select, by = join_by(assembly == accession)) %>%
  mutate(checkmInfo.completeness = as.numeric(checkmInfo.completeness)) %>%
  mutate(checkmInfo.contamination = as.numeric(checkmInfo.contamination)) %>%
  mutate(assemblyStats.totalSequenceLength = as.numeric(assemblyStats.totalSequenceLength)) %>%
  mutate(assemblyStats.gcPercent = as.numeric(assemblyStats.gcPercent)) %>%
  mutate(assemblyStats.numberOfComponentSequences = as.numeric(assemblyStats.numberOfComponentSequences))

df_smpB = read_nhmmer("smpB_HMMER_5e2.csv")
df_ssrA_most = read_nhmmer("ssrA_HMMER_5e2.csv")
df_ssrA_alphas = read_nhmmer("ssrA_HMMER_alpha_5e2.csv")
df_ssrA = bind_rows(df_ssrA_most, df_ssrA_alphas)
df_rqcH = read_nhmmer("rqcH_HMMER_5e2.csv")
df_arfA = read_nhmmer("arfA_HMMER_5e2.csv")
df_arfB = read_nhmmer("arfB_HMMER_5e2.csv")
df_mutS2 = read_nhmmer("mutS2_HMMER_5e2.csv")
df_smrB = read_nhmmer("smrB_HMMER_5e2.csv")

df_ssrA = filter(df_ssrA, abs(as.numeric(df_ssrA$hmmfrom) - df_ssrA$hmmto) > 200)
df_arfB = filter(df_arfB, abs(as.numeric(df_arfB$hmmfrom) - df_arfB$hmmto) > 150)
df_mutS2 = filter(df_mutS2, abs(as.numeric(df_mutS2$hmmfrom) - df_mutS2$hmmto) > 1250)
df_smrB = filter(df_smrB, (abs(as.numeric(df_smrB$hmmfrom) - df_smrB$hmmto) > 400))

df = data.frame(cbind(df_GCF_nc_tax , smpB = (df_GCF_nc_tax$nuccore %in% df_smpB$accessions), ssrA = (df_GCF_nc_tax$nuccore %in% df_ssrA$accessions), rqcH = (df_GCF_nc_tax$nuccore %in% df_rqcH$accessions), arfA = (df_GCF_nc_tax$nuccore %in% df_arfA$accessions), arfB = (df_GCF_nc_tax$nuccore %in% df_arfB$accessions), mutS2 = (df_GCF_nc_tax$nuccore %in% df_mutS2$accessions), smrB = (df_GCF_nc_tax$nuccore %in% df_smrB$accessions)))


df_presence = df %>%
  group_by(assembly) %>%
  summarise(across(smpB:smrB, ~sum(.))) %>%
  mutate_if(is.numeric, ~1 * (. != 0))

bools = ifelse(df_presence[-1] == 1,"TRUE","FALSE")
bools = data.frame(cbind(assembly = df_presence$assembly, bools))

names =  df %>%
  distinct(nuccore, .keep_all = TRUE) %>%
  group_by(assembly) %>%
  summarize(nuccore=paste(nuccore, collapse=";"))

# Final dataframe w assembly and nuccore accessions, 
df_distinct = inner_join(names, bools, by = "assembly") %>%
  left_join(df[,1:11], by = "assembly", multiple = "any") %>%
  select(-X, -nuccore.y) %>%
  rename(nuccore = nuccore.x) %>% 
  inner_join(lineages, by = join_by("organism.taxId" == "taxID"), multiple = "all") %>%
  filter(!grepl("Archaea", domain)) %>%
  filter(checkmInfo.completeness > 80) %>%
  filter(checkmInfo.contamination < 10) %>%
  distinct(assembly, .keep_all = TRUE)

df_distinct$phylum[df_distinct$phylum == "delta/epsilon subdivisions"] = "Pseudomonadota"
df_distinct$phylum[df_distinct$phylum == "Pseudomonadota"] = "Proteobacteria"

df_distinct$phylum[df_distinct$phylum == "Terrabacteria group"] = df_distinct$class[df_distinct$phylum == "Terrabacteria group"]
df_distinct$phylum[df_distinct$phylum == "Bacillota"] = "Firmicutes"
df_distinct$phylum[df_distinct$phylum == "Actinomycetota"] = "Actinobacteriota"
df_distinct$phylum[df_distinct$phylum == "Abditibacteriota"] = "Armatimonadota"

df_distinct$phylum[df_distinct$phylum == "Aquificota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Campylobacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Deferribacterota"] = "Aquificota + Campylobacterota + Deferribacterota"

df_distinct$phylum[df_distinct$phylum == "Thermodesulfobacteriota"] = "Desulfuromonadota + Desulfobacterota"
df_distinct$phylum[df_distinct$phylum == "Calditrichota"] = "FCB group"
df_distinct$phylum[df_distinct$phylum == "Thermomicrobiota"] = "Chloroflexota"

df_distinct$phylum[df_distinct$phylum == "Proteobacteria"] = df_distinct$class[df_distinct$phylum == "Proteobacteria"]


df_no_tmrna = df_distinct %>%
  filter(smpB == "FALSE" & ssrA == "FALSE")

df_nothing = df_distinct %>%
  filter(smpB == "FALSE" & ssrA == "FALSE" & rqcH == "FALSE" & smrB == "FALSE" & mutS2 == "FALSE" & arfA == "FALSE" & arfB == "FALSE")

df_nothing_complete = df_nothing %>%
  filter(assemblyInfo.assemblyLevel == "Complete Genome" | assemblyInfo.assemblyLevel == "Chromosome")

min(df_distinct$assemblyStats.totalSequenceLength)

ggplot(df_no_tmrna, aes(x = checkmInfo.completeness)) + 
  geom_histogram(bins = 10) +
  xlab("CheckM Completeness")


write.csv(df_distinct, file = "C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//ribo_rescue_df.csv")
