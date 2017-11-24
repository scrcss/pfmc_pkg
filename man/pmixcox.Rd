\name{pmixcox}
\alias{pmixcox}
\title{
IRLS-EM algorithm
}
\description{
  IRLS-EM algorithm for the optimization of penalized maximum likelihood
of finite-mixture Cox PH model using the penalty of LASSO, adaptive LASSO, SCAD or elastic net.
}
\usage{
pmixcox(Time, Delta, X, K = 3, iter.max = 1000, u.init = NULL, 
        tparm = 0, alpha = 1, scad = FALSE, adpcoef = NULL, 
		abstol = 0.01, reltol = 1e-04, seed = 1)
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
\item{K}{
  number of components in the finite-mixture Cox model
}
\item{iter.max}{
  maximum number of EM iterations
}
  \item{u.init}{
initial value of U, a data matrix that gives the probability that each sample belongs to each component. Each row corresponds to one sample and each column indicates one component.
}
\item{tparm}{
value of tuning parameter for variable selection
}
\item{alpha}{
  the elastic net mixing parameter
}
\item{scad}{
  is the SCAD penalty applied?
}
\item{adpcoef}{
  the adaptive weights in the adaptive LASSO method
}
\item{abstol}{
  absolute tolerance of EM algorithm
}
\item{reltol}{
relative tolerance of EM algorithm
}
\item{seed}{
  random seed for initialing U if it is not given
}
}
\value{
  a list, where
  \item{U}{posterior probability matrix}
  \item{fit}{a list with \code{K} components, each of which is an object with class \code{glmnet} or \code{ncvsurv}, output by the last EM iteration}
  \item{pi}{estimate for the mixing probability}
  \item{class}{subtype classification for each sample}
  \item{ploglik}{partial log-likelihood}
  \item{mixloglik}{observed log-likelihood}
  \item{iter}{number of EM iterations}
  \item{convergence}{indicator to show the algorithm converged or not}
}
\references{
  Subtype classification and heterogeneous prognosis model construction in
precision medicine. Na You, Shun He, Xueqin Wang, Junxian Zhu and Heping Zhang
}
\author{
  Shijie Quan, Shun He
}
\examples{
require(PFMC)
data(example.data)
nopenalty.fit = pmixcox(Time = example.data$time,
  Delta = example.data$status,
  X = example.data[,-1:-2],
  K = 2,
  tparm = 0,
  seed = 1)

adpcoef = lapply(nopenalty.fit$fit, function(x) 1/abs(x))
alasso.cv = calcvcox(Time = example.data$time,
  Delta = example.data$status,
  X = example.data[,-1:-2],
  K = 2,
  nopenaltyresult = nopenalty.fit,
  adpcoef = adpcoef)

tparmchoice = alasso.cv[1,]
alasso.fit = pmixcox(Time = example.data$time,
  Delta = example.data$status,
  X = example.data[,-1:-2],
  K = 2,
  u.init = nopenalty.fit$U,
  tparm = tparmchoice,
  adpcoef = adpcoef)

}