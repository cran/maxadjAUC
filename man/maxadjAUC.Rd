\name{maxadjAUC}
\alias{maxadjAUC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximizing the Covariate-Adjusted AUC
}
\description{
Often there is interest in combining several predictors or biomarkers into a linear combination for diagnosis, prognosis or screening. This can be done by targeting measures of predictive capacity. In the presence of a discrete covariate, such as batch or study center, an appropriate summary of discriminatory performance is the covariate-adjusted area under the receiver operating characteristic curve (AUC), or aAUC. This function estimates a linear combination of predictors by maximizing a smooth approximation to the aAUC.
}
\usage{
maxadjAUC(outcome, predictors, covariate, initialval="rGLM", approxh = 1/3,
conditional=FALSE, tolval = 1e-6, stepsz = 1e-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outcome}{
    A vector of outcome (disease) indicators for each observation (1 for diseased, 0 for non-diseased). Missing values are not allowed.
}
\item{predictors}{
    A numeric matrix with one row for each observation and one column for each candidate predictor. Missing values are not allowed. The columns of the matrix will be (re)named "V1", "V2", ....
}
\item{covariate}{
    A numeric vector of covariate values for each observation. The covariate should have a limited number of values (i.e., it should be a discrete covariate). Missing values are not allowed.
}
\item{initialval}{
    Starting values of the predictor combination for the SaAUC algorithm. Default value is \code{"rGLM"}, which means that estimates from robust logistic regression, specifically the method of Bianco and Yohai (implemented via the \code{aucm} package), are used as starting values. If any other value of \code{initialval} is given, or if robust logistic regression fails to converge, estimates from standard logistic regression are used as starting values. For both robust and standard logistic regression, the covariate will be included as a stratifying variable.
}
\item{approxh}{
    The tuning parameter for the smooth approximation to the covariate-specific AUC is the ratio of the standard deviation of the linear combination (based on the starting values) to \eqn{n_c^{approxh}}, where \eqn{n_c} is the number of observations with covariate value \eqn{c}. In particular, larger values of \code{approxh} will provide a better approximation to the AUC, though estimation may become unstable if \code{approxh} is too large. Default 1/3.
}
\item{conditional}{
    A logical value indicating whether standard logistic regression should be conditional if TRUE (i.e., \code{survival::clogit}) or unconditional if FALSE (\code{stats::glm}). Default is FALSE.
}
\item{tolval}{
    Controls the tolerance on feasibility and optimality for the optimization procedure (performed by \code{solnp} in the \code{Rsolnp} package). Default 1e-6.
}
\item{stepsz}{
    Controls the step size for the optimization procedure (performed by \code{solnp} in the \code{Rsolnp} package). Default 1e-5.
}
}
\details{
The function seeks to optimize a smooth approximation to the covariate-adjusted AUC,
\eqn{SaAUC = \sum_{c=1}^m w_c SAUC_c}
where \eqn{SAUC_c} are the smooth approximations to the covariate-specific AUC and \eqn{w_c} are covariate-specific weights for a covariate with \eqn{m} values in the data.
}

\value{
A list will be returned with the following components:
\item{NumCov }{The number of covariate strata used after removing concordant strata.}
\item{FittedCombs }{A list containing four fitted combinations: \code{InitialVal} (either robust logistic regression, if \code{initialval="rGLM"}, standard unconditional logistic regression, if \code{initialval} is not \code{"rGLM"} and \code{conditional=FALSE}, or standard conditional logistic regression, if \code{initialval} is not \code{"rGLM"} and \code{conditional=TRUE}), \code{NormGLM} (standard unconditional or conditional logistic regression, depending on \code{conditional}), \code{NormrGLM} (robust logistic regression), \code{MaxSaAUC} (SaAUC approach). All fitted combination vectors are normalized. If robust logistic regression fails to converge, standard unconditional or conditional logistic regression (depending on \code{conditional}) is used instead, and a warning is given.}
\item{aAUCTR }{A vector of the aAUC in the training data for the four fitted combinations.}
\item{varTR }{A vector of the variability in the covariate-specific AUCs around the aAUC in the training data for the four fitted combinations.}
}
\references{
Bianco, A.M. and Yohai, V.J. (1996) Robust estimation in the logistic regression model. In Robust statistics, data analysis, and computer intensive methods (ed H. Rieder), pp 17-34. Springer.

Janes, H. and Pepe, M.S. (2009) Adjusting for covariate effects on classification accuracy using the covariate-adjusted
receiver operating characteristic curve. Biometrika, pages 1-12.

Meisner, A., Parikh, C.R., and Kerr, K.F. (2017). Developing biomarker combinations in multicenter studies via direct maximization and penalization. UW Biostatistics Working Paper Series, Working Paper 421.
}
% \author{
%%  ~~who you are~~
%}
\note{
The function automatically removes any covariate strata that are concordant on the outcome (i.e., all 0 or all 1).

Warnings are issued if the SaAUC algorithm does not converge or if robust logistic regression fails to converge.

The standard unconditional or conditional logistic regression algorithm may not converge, producing a warning. If such a convergence failure occurs, the "GLM" results will be affected, as will the "rGLM" results if the robust logistic model also fails to converge, and the "SaAUC" results if \code{initialval} is not "rGLM" or if the robust logistic model fails to converge. Thus, users should be alert to any convergence failures.
}

\seealso{
 \code{\link[aucm]{rlogit}}, \code{\link[Rsolnp]{solnp}}
}
\examples{
  \dontshow{
  expit <- function(x){
    exp(x)/(1+exp(x))
  }

  set.seed(1)

  covar <- rep(c(1:2),each=30)

  x1 <- rnorm(60,0,rep(runif(2,0.8,1.2),each=30))
  x2 <- rnorm(60,0,rep(runif(2,0.8,1.2),each=30))

  covint <- rep(runif(2,-1.5,1.5), each=30)

  y <- rbinom(60,1,expit(covint + 1*x1 - 1*x2))
  X <- cbind(x1,x2)

  output <- maxadjAUC(outcome=y, predictors=X, covariate=covar, initialval="rGLM",
                      approxh = 1/3, conditional=FALSE, tolval = 1e-2, stepsz = 1e-2)
  output
  }
  \donttest{
  expit <- function(x){
    exp(x)/(1+exp(x))
  }

  set.seed(1)

  covar <- rep(c(1:4),each=100)

  x1 <- rnorm(400,0,rep(runif(4,0.8,1.2),each=100))
  x2 <- rnorm(400,0,rep(runif(4,0.8,1.2),each=100))
  x3 <- rnorm(400,0,rep(runif(4,0.8,1.2),each=100))
  x4 <- rnorm(400,0,rep(runif(4,0.8,1.2),each=100))

  covint <- rep(runif(4,-1.5,1.5), each=100)

  y <- rbinom(400,1,expit(covint + 1*x1 - 1*x2 + 1*x3 - 1*x4))
  X <- cbind(x1,x2,x3,x4)

  output <- maxadjAUC(outcome=y, predictors=X, covariate=covar, initialval="rGLM",
                      approxh = 1/3, conditional=FALSE, tolval = 1e-6, stepsz = 1e-5)
  output
  }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
