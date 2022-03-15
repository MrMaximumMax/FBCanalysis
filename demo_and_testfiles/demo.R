#Demo to test and let run the functionalities e.g. on a local machine with the provided test csv files form this folder

library(imputeTS)
library(lubridate)
library(arsenal)
library(FCPS)
library(dplyr)
library(readr)
library(emdist)
library(cluster)
library(RankAggreg)
library(glmnet)
library(clValid)
library(mclust)

#get the FBCanalysis package
install.packages("devtools")
devtools::install_github("MrMaximumMax/FBCanalysis")
#Provided demo files can be found in the demo directory
library(FBCanalysis)

#A) Data preparation and processing
list <- patient_list('.../ts_demofiles1') #no specific csv.-file, just insert
#folder where csv.-files they are stored
#Sample frequency is daily
patient_ts_plot(list,"testpat_1","PEF")
patient_boxplot(list,c("ID_2","testpat_1","testpat_2","a301"), "FEV1")
patient_hist(list,"testpat_1","PEF")

#B) Generate EMD data
#In case EMD function cannot handle max_iter, proceed according to following
#instructions: https://github.com/s-u/emdist/issues/2
matrix <- emd_matrix(list, "FEV1")
emd_heatmap(matrix)
max_fluc(list, "PEF")

#C) Time series data clustering and enrichment analysis
clustering <- clust_matrix(matrix, method = "kmeans", nclust = 3)
enr <- add_enrich(list, '.../enrichment_dat.csv') #check directory for demo file
enr <- add_clust2enrich(enr, clustering)
ts <- add_clust2ts(list, clustering)
enr_obs_clust(ts, enr, 1)
path <- ".../enrichment_dat.csv"
test <- sim_sample_enr(list,path,clustering,1,100)

#D) Cluster validation measure analysis
list <- patient_list('.../ts_demofiles2') #sample frequency is twice daily
distmat <- emd_matrix(list, "PEF", maxIter = 5000)
parameters <- init_clValid()
output <- clValid_flow(distmat, parameters)
clustdat <- clust_matrix(distmat, output$method, as.numeric(output$clust_num))
enr <- add_enrich(list, '.../enrichment_dat.csv')
enr <- add_clust2enrich(enr, clustdat, output$method)
ts <- add_clust2ts(list, clustdat, output$method)
enr_obs_clust(ts, enr, 1)

#E) Random data removal/Jaccard
testlist <- patient_list('.../ts_demofiles1')
testlist_rm <- rnd_dat_rm(testlist, 0.95)
output <- sim_jaccard_emd(testlist, "PEF", 0.05, 10, "hierarchical", 2)
output <- jaccard_run_emd(testlist,"PEF",10,"hierarchical",1,3,c(0.001,0.005,0.01,0.05))

output <- sim_jaccard_cognate(testlist,"PEF",0.1,10,"hierarchical",3)
output <- jaccard_run_cognate(testlist,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
