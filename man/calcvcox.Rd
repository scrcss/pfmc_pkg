\name{calcvcox}
\alias{calcvcox}
\title{
Cross-validation for the penalized finite-mixture Cox PH model
}
\description{
  Cross-validation for the choice of tuning parameters in the penalized
  finite-mixture Cox regression method
}
\usage{
calcvcox(Time, Delta, X, ncore = 2, K = 2, nfolds = 10, 
         nopenaltyresult = NULL, alpha = 1, scad = FALSE, 
		 adpcoef = NULL, nlambda = 100, seed = 1)
}
\arguments{
  \item{Time}{
observed time
}
  \item{Delta}{
survival status
}
  \item{X}{
a data matrix of explanatory variables, where each colomn
correponds to one predictor and each row indicates one sample.   
}
  \item{ncore}{
the number of CPU cores occupied for parallele implementation
}
  \item{K}{
number of components in the finite-mixture Cox model
}
  \item{nfolds}{
number of folds split in CV partitioning
}
\item{nopenaltyresult}{
  fitting of the full finite-mixture Cox PH model, which is a list with
  \code{U}, a matrix of the posterior probabilities of each
  sample belonging to each component, and \code{pi}, the estimate
  for the mixing probability
}
\item{alpha}{
  the elastic net mixing parameter
}
  \item{scad}{
is the SCAD penalty applied?
}
\item{adpcoef}{
  the adaptive weight of adaptive LASSO method
}
\item{nlambda}{
the number of \code{lambda} values
}
  \item{seed}{
random seeding for CV split
}
}
\value{
  a data matrix, where the first row \code{lambda.min} presents the
  value of tuning parameter that gives the minimum CV error in each component
  of the finite-mixutre Cox PH model, and the second row
  \code{lambda.1se} provides the largest value of tuning parameter such
  that CV error is within 1 standard error of the minimum in each component.
}
\references{
Subtype classification and heterogeneous prognosis model construction in
precision medicine. Na You, Shun He, Xueqin Wang, Junxian Zhu and Heping Zhang
}
\author{
Shijie Quan, Shun He
}
\seealso{
pmixcox
}
\examples{
## see pmixcox ##
}

