% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jaccard.R
\name{reap_freq}
\alias{reap_freq}
\title{Determine reappearance frequency}
\usage{
reap_freq(plist, parameter, removal, n_simu, method, n_clust, Iter, normalize)
}
\arguments{
\item{plist}{Object of type list storing patient time series data (also see function: \link{patient_list})}

\item{parameter}{Parameter of interest in time series data list}

\item{removal}{Amount of random data removal to determine Jaccard index}

\item{n_simu}{Number of simulations}

\item{method}{Clustering method (also see function: \link{clust_matrix})}

\item{n_clust}{Number of clusters (also see function: \link{clust_matrix})}

\item{Iter}{Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
default is 5,000 for this function}

\item{normalize}{Indicates if parameter indicated needs to be normalized or not (TRUE by default)}
}
\value{
Vector of length of n_simu where reappearance frequency is stored for each simulation run
}
\description{
Determine the reappearance frequency for clustered elements in perturbed data
for specific amount of random data removal
}
\details{
To begin, each participant's collection of measured values is randomly
depleted of a certain proportion of measurements. The clustering technique is
then performed using the perturbed data, and the resulting clusters are compared
to the original clusters created with the unperturbed gold standard.
This technique is done iteratively in order to provide statistics indicating
the original clusters' stability following the elimination of random data.
These stability statistics are calculated using two cluster similarity metrics:
Jaccard's index, a measure of global similarity that quantifies the extent to
which the original and modified clusters overlap. Additionally, it is utilized
to decide which cluster is considered the cognate cluster. Then it was transformed
into a local measure, meaning the frequency with which each member of the original
clusters reappeared between iterations.
}
\references{
Edgar Delgado-Eckert, Oliver Fuchs, Nitin Kumar, Juha Pekkanen,
Jean-Charles Dalphin, Josef Riedler, Roger Lauener, Michael Kabesch, Maciej Kupczyk,
Sven-Erik Dahlen, et al. Functional phenotypes determined by fluctuation-based
clustering of lung function measurements in healthy and asthmatic cohort participants.
Thorax, 73(2):107–115, 2018.
}
