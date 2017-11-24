\name{PFMC-internal}
\title{Internal PFMC functions}
\alias{cv.ncvsurv}
\alias{ncvsurv}
\alias{convexMin}
\alias{getmin}
\alias{lamNames}
\alias{setupLambdaCox}
\alias{std}
\description{Internal PFMC functions}
\usage{
cv.ncvsurv(X, y, ..., cluster, nfolds=10, seed, returnY=FALSE, trace=FALSE)
ncvsurv(X, y, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, 3),
        alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, 
		eps=1e-4, max.iter=10000, convex=TRUE, dfmax=p, 
		penalty.factor=rep(1, ncol(X)),weights=rep(1,length(y[,1])), 
		warn=TRUE, returnX=FALSE, ...)
convexMin(b, X, penalty, gamma, l2, penalty.factor, a, Delta=NULL, weights)
ncvgetmin(lambda,cvm,cvsd)
lamNames(l)
setupLambdaCox(X, y, Delta, alpha, lambda.min, nlambda, 
           penalty.factor, weights=rep(1,length(Delta))) 
std(X, weights)
}
\references{
  R package \code{ncvreg}, https://cran.r-project.org/web/packages/ncvreg/
  
  Breheny, P. and Huang, J. (2011) Coordinate descent algorithms for
  nonconvex penalized regression, with applications to biological
  feature selection. Ann. Appl. Statist., 5: 232-253.
  
  Subtype classification and heterogeneous prognosis model construction
  in precision medicine. Na You, Shun He, Xueqin Wang, Junxian Zhu and Heping Zhang
}
\author{
  Shijie Quan, Shun He
}

\details{
Functions from R package \code{ncvreg}, which are adapted to allow the
argument \code{weights} on samples.
}
\keyword{internal}
