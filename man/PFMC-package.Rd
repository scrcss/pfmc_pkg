\name{PFMC-package}
\alias{PFMC-package}
\alias{PFMC}
\docType{package}
\title{
  IRLS-EM algorithm for the penalized finite-mixture Cox PH model
}
\description{
  A variable selection method for the finite-mixture Cox PH model using
  the penalty of LASSO, adaptive LASSO, SCAD or elastic net
}
\details{
  The internal funtions \code{convexMin}, \code{cv.ncvsurv}, \code{cvf.surv},
  \code{lamNames}, \code{ncvgetmin}, \code{ncvsurv},
  \code{setupLambdaCox} and \code{std} are imported from R package
  \code{ncvreg} and adapted to allow the argument \code{weights}.
}
\author{
  Shijie Quan, Shun He
  
  Maintainer: Na You (youn@mail.sysu.edu.cn)
}
\references{
  R package \code{ncvreg}, https://cran.r-project.org/web/packages/ncvreg/
  
  Breheny, P. and Huang, J. (2011) Coordinate descent algorithms for
  nonconvex penalized regression, with applications to biological
  feature selection. Ann. Appl. Statist., 5: 232-253.
  
  Subtype classification and heterogeneous prognosis model construction
  in precision medicine. Na You, Shun He, Xueqin Wang, Junxian Zhu and
  Heping Zhang. 
}
\keyword{PFMC}
