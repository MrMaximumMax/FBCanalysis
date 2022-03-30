#Packages required to run the functions from the FBCanalysis package
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
#Installation via GitHub required since no official release so fat
install.packages(c("devtools","processx")) #In case they are not installed yet
devtools::install_github("MrMaximumMax/FBCanalysis")
#Load the library
library(FBCanalysis)
#Get overview of the package and its functions
#?function_name can also be called
help("FBCanalysis")
?sim_sample_enr

#A) Data preparation, processing and visualization
list <- patient_list("https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv", GitHub = TRUE)
#Sample frequency is daily (twice a day would work likewise);
#alternatively work with following real world beta data: https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/beta/phys/data.csv
#In real world beta data, FEV1 is normalized, PEF is not normalized -> the further functions
#can handle both cases!
#Day:Month:Year time format for demo; Year:Month:Day for real world beta

#Visualize the data stored in the list
patient_ts_plot(list,"a301","PEF")
patient_boxplot(list,c("a301","ID_2"), "FEV1", normalize = TRUE) #normalize TRUE/FALSE can be used here likewise
#By default, normalize = TRUE; This indicates that data does not have to be
#normalized before calling patient_list(); In case data has been normalized
#before, set normalize = FALSE
patient_hist(list,"ID_2","PEF")

#B) Generate EMD data
#In case EMD function cannot handle max_iter, proceed according to following
#instructions: https://github.com/s-u/emdist/issues/2
#The user can indicate if the data needs to be normalized or not
matrix <- emd_matrix(list, "PEF", maxIter = 5000, normalize = TRUE)
emd_heatmap(matrix)
#this would also work:
emd_heatmap(list,"FEV1",normalize = FALSE)
?emd_matrix
#Find Patient_ID pair data distribution with largest EMD and visualize boxplot
max_fluc(list, "PEF", maxIter = 5000)
?max_fluc

#C) Time series data clustering and enrichment analysis (with data preprocessing)
clustering <- clust_matrix(matrix, method = "kmeans", nclust = 2)
enr <- add_enrich(list, 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv')
#alternatively work with following real world beta data: https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/beta/enrich/enrichment.csv
enr <- add_clust2enrich(enr, clustering)
ts <- add_clust2ts(list, clustering)
enr_obs_clust(ts, enr, 1, numeric = "anova", categorical = "fe")
path <- 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv'
test <- sim_sample_enr(list,path,clustering,1,100, numeric = "kwt", categorical = "chisq")

#D) Cluster validation measure analysis
#Initialize the clustering techniques and number of clusters of interest
parameters <- init_clValid()
#Perfrom analysis with the defined techniques and cluster numbers
output <- clValid_flow(matrix, parameters)
#Now transfer the output to the enrichmebnt analysis (also random dat. removal possible)
enr <- add_enrich(list, 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv')
enr <- add_clust2enrich(enr, output)
ts <- add_clust2ts(list, output)
enr_obs_clust(ts, enr, 1, numeric = "anova", categorical = "chisq")
test <- sim_sample_enr(list,path,output,1,100, numeric = "anova", categorical = "fe")

#E) Random data removal/Jaccard
#Removes 50% of time series data
list_rm <- rnd_dat_rm(list, 0.5)
#Stability tested for specific amount of data removal; output generated on all clusters
?sim_jaccard_global
removal_out <- sim_jaccard_global(list, parameter =  "PEF", removal =  0.1,
                                  n_simu =  10, method = "hierarchical", n_clust = 2,
                                  normalize = TRUE)
#Stability tested for multiple removal amounts; output generated for specific cluster
removal_out <- jaccard_run_global(list, "PEF", 10, "hierarchical",clust_num = 1, n_clust = 3,
                                  c(0.005,0.01,0.05,0.1,0.8), normalize = TRUE, maxIter = 5000)

#Analogous data removal approach via Earth Mover's distances (Version 1)
removal_out <- sim_jaccard_emd(list, "PEF", 0.05, 10, "hierarchical", 2, normalize = TRUE,
                               maxIter = 5000)
removal_out <- jaccard_run_emd(list,"PEF",10,"kmeans",1,2,c(0.001,0.005,0.01,0.05))

#Analogous data removal approach via Earth Mover's distances (Version 2)
removal_out <- sim_jaccard_emd_2(list, "FEV1", 0.1, 10, "pam", 2, normalize = TRUE,
                               maxIter = 5000)
removal_out <- jaccard_run_emd_2(list,"FEV1",10,"pam",1,2,c(0.01,0.05,0.1,0.2,0.3))

#Check reappearance frequency
removal_reap <- reap_freq(list,"PEF",0.01,10,"hierarchical",2,5000)
removal_reap <- reap_freq_run(list,"PEF",10,"hierarchical",3,c(0.01,0.05,0.1,0.2,0.3,0.4))

