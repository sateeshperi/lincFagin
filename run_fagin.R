library(fagin)                                                                     
library(rmonad)                                                                                                 
library(knitr)                                                                                                  
library(magrittr)                                                                                               
library(readr)                                                                                                  

get_brassicaceae_config <- function(){                                                                          
  con <- fagin::config()                                                                                        
  con@archive = "brass-test"                                                                                    
  con@synder@offsets = c(0L,1L)
  con@synder@trans = "d" # proportion transform satsuma
  con@alignment@dna2dna_maxspace = 1e8L
  con@input@focal_species = "Arabidopsis_thaliana"
  con@input@gff <- list(
    Arabidopsis_thaliana = "gff/Arabidopsis_thaliana.gff"
    , Eutrema_salsugineum  = "gff/Eutrema_salsugineum.gff"
  )
  con@input@fna <- list(
    Arabidopsis_thaliana = "fna/Arabidopsis_thaliana.fna"
    , Eutrema_salsugineum  = "fna/Eutrema_salsugineum.fna"
  )
  con@input@syn <- list(
    Eutrema_salsugineum = "syn/Arabidopsis_thaliana.vs.Eutrema_salsugineum.syn"
  )
  con@input@tree <- "tree.newick"
  con@input@query_gene_list <- "query.txt"
  con@input@control_gene_list <- "control.txt"
  fagin::validate_config(con)
  con
}

con <- get_brassicaceae_config()

m <- run_fagin(con)
