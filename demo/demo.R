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
list <- FBCanalysis::patient_list("/Users/maximiliankohler/Documents/MasterBiomedEng/Thesis/Package/FBCanalysis/demo/dummy") #no specific csv.-file, just insert
#folder where csv.-files they are stored
#Sample frequency is daily
patient_ts_plot(list,"aPSOW","PEF")
patient_boxplot(list,c("aPSOW","0s03n"), "FEV1")
patient_hist(list,"0s03n","PEF")

#B) Generate EMD data
#In case EMD function cannot handle max_iter, proceed according to following
#instructions: https://github.com/s-u/emdist/issues/2
matrix <- emd_matrix(list, "PEF", maxIter = 5000)
emd_heatmap(matrix)
max_fluc(list, "PEF", maxIter = 5000)

#C) Time series data clustering and enrichment analysis
clustering <- clust_matrix(matrix, method = "kmeans", nclust = 2)
enr <- add_enrich(list, '/Users/maximiliankohler/Documents/MasterBiomedEng/Thesis/Package/FBCanalysis/demo/enrichment_dummy/enrichment.csv') #check directory for demo file
enr <- add_clust2enrich(enr, clustering)
ts <- add_clust2ts(list, clustering)
enr_obs_clust(ts, enr, 1)
path <- '/Users/maximiliankohler/Documents/MasterBiomedEng/Thesis/Package/FBCanalysis/demo/enrichment_clin/enrichment.csv'
test <- sim_sample_enr(list,path,clustering,1,100)

#D) Cluster validation measure analysis
parameters <- init_clValid()
output <- clValid_flow(matrix, parameters)
enr <- add_enrich(list, '/Users/maximiliankohler/Documents/MasterBiomedEng/Thesis/Package/FBCanalysis/demo/enrichment_dummy/enrichment.csv')
enr <- add_clust2enrich(enr, output)
ts <- add_clust2ts(list, output)
enr_obs_clust(ts, enr, 1)

#E) Random data removal/Jaccard

list_rm <- rnd_dat_rm(list, 0.95)
output <- sim_jaccard_emd(list, "PEF", 0.05, 3, "hierarchical", 2, Iter = 5000)
output <- jaccard_run_emd(list,"PEF",10,"hierarchical",1,3,c(0.001,0.005,0.01,0.05))

output <- sim_jaccard_cognate(list,"PEF",0.1,10,"hierarchical",3)
output <- jaccard_run_cognate(list,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
