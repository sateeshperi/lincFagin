library(tidyr)
library(dplyr)
library(fagin)
library(ggplot2)

get_tag(m)

rbind_with_name <- function(xs, grpname){
  for(name in names(xs)){
    xs[[name]][[grpname]] <- name
  }
  out <- do.call(rbind, xs)
  rownames(out) <- NULL
  out
}


# convert the group labels corresponding to the fagin paper
name_conversion <-  c(O1="A_gen", O2="A_trn", O3="A_orf",
                      N1="N_cds", N2="N_exo", N3="N_rna", N4="N_dna",
                      U2="U_ind", U5="U_scr", U6="U_unk", U1="U_una", U3="U_nst", U7="U_tec")

# function to clean the species name
clean_name <- function(x, suffix){
  x$target_species <- sub(suffix, "", x$target_species)
  x
}

# determine the homolog feature for each gene as figure 3 in fagin paper 
feats <- get_value(m, tag='feature_table') %>%
  rbind_with_name("target_species") %>%
  clean_name("feature_table/query/") %>% 
  clean_name("feature_table/control/")

# load phylostrata levels for query and control gene list
# see phylostratr pacakge for the format
# strata <- readr::read_tsv("strata.tab")

# classify homolog group for query gene list

id_query <- get_value(m, tag='query_genes')[[1]]

#query <- subset(strata, seqid %in% id_query) %>%
#  dplyr::select(seqid, std_strata_name = mrca, std_strata_level = ps)

vs_query <- get_value(m, tag='query_labels')[[1]]

labels_query <- rbind_with_name(vs_query$labels, "target_species") %>%
  dplyr::select(seqid, homology_class=secondary, target_species) %>%
  dplyr::mutate(homology_class = name_conversion[homology_class])

query <- merge(feats,labels_query, by=c("seqid", "target_species"))

query <- merge(query, feats, by=c("seqid", "target_species"))

ori_query <- get_value(m, tag="query_origins")[[1]]$backbone
ori_query <- as.matrix(ori_query)
ori_query[ori_query == "O"] <- "A"
ori_query <- as.data.frame(ori_query)
ori_query$seqid <- rownames(ori_query)

query <- merge(query, ori_query, by="seqid")
query$group = "query"

# classify homolog group for control gene list
id_control <- get_value(m, tag='control_genes')[[1]]

#control <- subset(strata, seqid %in% id_control) %>%
# dplyr::select(seqid, std_strata_name = mrca, std_strata_level = ps)

vs_control <- get_value(m, tag='control_labels')[[1]]


labels_control <- rbind_with_name(vs_control$labels, "target_species") %>%
  dplyr::select(seqid, homology_class=secondary, target_species) %>%
  dplyr::mutate(homology_class = name_conversion[homology_class])

control <- merge(feats, labels_control, by=c("seqid", "target_species"))

control <- merge(control, feats, by=c("seqid", "target_species"))

ori_control <- get_value(m, tag="control_origins")[[1]]$backbone
ori_control <- as.matrix(ori_control)
ori_control[ori_control == "O"] <- "A"
ori_control <- as.data.frame(ori_control)
ori_control$seqid <- rownames(ori_control)

control <- merge(control, ori_control, by='seqid')
control$group = "control"

# combine query gene and control gene result into a dataframe
fagin_result <- rbind(query, control)

# barplot for homolog class as figure 5 in fagin paper
ggplot(query, aes(x=homology_class, fill=target_species)) +
  geom_bar(stat="count", position=position_dodge()) +
  theme_classic() 

ggplot(control, aes(x=homology_class, fill=target_species)) +
  geom_bar(stat="count", position=position_dodge()) +
  theme_classic() 





get_tag(m)

test <- get_value(m, tag='feature_table')
names(test)
