#Packages required to run the functions
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
install.packages(c("devtools","processx"))
devtools::install_github("MrMaximumMax/FBCanalysis")
library(FBCanalysis)

#A) Data preparation and processing
list <- patient_list("https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv", GitHub = TRUE)
#Sample frequency is daily (twice a day would work likewise); Day:Month:Year time format
patient_ts_plot(list,"a301","PEF")
patient_boxplot(list,c("a301","ID_2"), "FEV1")
patient_hist(list,"ID_2","PEF")

#B) Generate EMD data
#In case EMD function cannot handle max_iter, proceed according to following
#instructions: https://github.com/s-u/emdist/issues/2
matrix <- emd_matrix(list, "PEF", maxIter = 5000)
emd_heatmap(matrix)
max_fluc(list, "PEF", maxIter = 5000)

#C) Time series data clustering and enrichment analysis
clustering <- clust_matrix(matrix, method = "kmeans", nclust = 2)
enr <- add_enrich(list, 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv')
enr <- add_clust2enrich(enr, clustering)
ts <- add_clust2ts(list, clustering)
enr_obs_clust(ts, enr, 1)
path <- 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv'
test <- sim_sample_enr(list,path,clustering,1,100)

#D) Cluster validation measure analysis
parameters <- init_clValid()
output <- clValid_flow(matrix, parameters)
enr <- add_enrich(list, 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv')
enr <- add_clust2enrich(enr, output)
ts <- add_clust2ts(list, output)
enr_obs_clust(ts, enr, 1)
test <- sim_sample_enr(list,path,output,1,100)

#E) Random data removal/Jaccard
list_rm <- rnd_dat_rm(list, 0.5) #Removes 50% of time series data
output <- sim_jaccard_cognate(list,"PEF",0.1,10,"hierarchical",3)
output <- jaccard_run_cognate(list,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
output <- sim_jaccard_emd(list, "PEF", 0.05, 3, "hierarchical", 2, Iter = 5000)
output <- jaccard_run_emd(list,"PEF",10,"kmeans",1,2,c(0.001,0.005,0.01,0.05))

