% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_prep_and_process.R
\name{patient_list}
\alias{patient_list}
\title{Process patient time series data}
\usage{
patient_list(path, GitHub)
}
\arguments{
\item{path}{Path where csv file(s) are stored (only folder, not specific file(s))}

\item{GitHub}{Set TRUE when csv file comes form GitHub (FALSE by default); only in demo needed}
}
\value{
Object of type list storing patient time series data
}
\description{
Process patient time series data by interpolation options and store data in
an object of type list.
}
\details{
It should be outlined that no prior z-normalization is required!
Within FBCanalysis' functions, the user can indicate if normalization is
required or not.

Prior to undertaking an analysis using one of the FBC procedures,
it is necessary to adequately process and prepare the relevant time series data.
The function then creates an interactive flow using the console in R Studio.
To begin, the method retrieves all csv files in the provided folder, indicating
that it is capable of handling multiple files. The function extracts all csv files
from the given directory and merges them into a single raw data frame. The user
then indicates which column represents Patient ID and time for adequate processing.
The csv files are merged, columns are selected where the Patient ID column will
be renamed ”Patient_ID” and the time column will be titled ”Time”. This
standardization approach is critical for subsequent features because it enables
the easy detection of time series data and the consistent computation and
processing of data, for example z-normalization.

The user should also indicate in the interactive console the time formate which
will be standardized with the help of \link{lubridate}. This is crucial because
the technique can now filter the raw data by Patient ID, extract the start and
end timestamps for each Patient ID, and then align the data if any records are
missing while maintaining the indicated sample frequency.

The user may choose between seven approaches: L1 Regularization/Least absolute
shrinkage and selection operator (LASSO) Regression, L2 Regularization/Ridge
Regression, Elastic Net Regularization, Linear interpolation, Cubic C2 interpolation
or, according to recent articles, fill in missing values using the highest or
lowest quartile of measurements in the given time series data distribution.

The Regression and Regularization techniques generate adequate polynomials for
each possible degree n-1 (where n is the total number of data points). Afterwards,
cross-validation (from \link{glmnet}) is applied to determine the lambda value
for the lowest MSE of the model. Afterwards, the model with polynomial degree
for the lowest MSE is chosen and the missing data is interpolated with the
regularized model.

It is also possible to apply a simple linear interpolation in between missing
time series data points. It may be the easiest option to employ straight lines
between neighboring points (also see \link{na_interpolation}). Nevertheless,
these basic spline polynomials may be notoriously inexact. Cubic spline polynomials
mostly provide better results.

Another option for the user is to apply the interpolation by using a cubic C spline.
It implies that the composite function S must be twice continuously differentiable
from all boundaries or subintervals (also see \link{na_interpolation}).

Without regression, regularization or interpolation, the user may opt to sample
missing values within time series data by randomly choosing a value from the
greatest or lowest quartile readings from each patient distribution. The R
function then loops over each NA element in the time series data distribution
of a patient for the specified parameter and randomly samples a value for the
chosen quartile until the data frame is complete.
}
\examples{
list <- patient_list(
"https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
GitHub = TRUE)
#Sampling frequency is supposed to be daily

}
\references{
Jerome Friedman, Trevor Hastie, Rob Tibshirani, Balasubramanian
Narasimhan, Ken- neth Tay, Noah Simon, and Junyang Qian. Package ‘glmnet’.
Journal of Statistical Software. 2010a, 33(1), 2021.

Hui Zou and Trevor Hastie. Regularization and variable selection via the elastic net.
Journal of the royal statistical society: series B (statistical methodology),
67(2):301– 320, 2005.

Steffen Moritz and Thomas Bartz-Beielstein. imputets: time series missing value
imputation in r. R J., 9(1):207, 2017.
}
