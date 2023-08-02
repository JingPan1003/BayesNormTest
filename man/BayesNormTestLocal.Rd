\name{BayesNormTestLocal}
\alias{BayesNormTestLocal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{BayesNormTestLocal
%%  ~~function to do ... ~~
}
\description{Bayesian hypothesis test for normality with local prior of gamma (skewness parameter). The null hypothesis H0 is that the data set follows normal distribution. The alternative hypothesis H1 is that the data set follows two-piece normal distribution (asymmetric distribution).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BayesNormTestLocal(data, m0, s0, a0, b0, m0l,s0l,initN,init2L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{the data set used for testing}
  \item{m0}{the mean of normal distribution for prior of mu}
  \item{s0}{the standard deviation of normal distribution for prior of mu}
  \item{a0}{the shape parameter of gamma distribution for prior of sigma}
  \item{b0}{the rate parameter of gamma distribution for prior of sigma}
  \item{m0l}{the mean of normal distribution for local prior of gamma}
  \item{s0l}{the standard deviation of normal distribution for local prior of gamma}
  \item{initN}{numeric vector, initial values for the parameters to be optimize for normal model}
  \item{init2L}{numeric vector, initial values for the parameters to be optimized for two-piece model}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{probNorm }{The probability that data set follows a normal distribution.}
  \item{probTPN }{The probability that data set follows a two-piece normal distribution.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
df <- swiss
swiss_dat <- df$Infant.Mortality
BayesNormTestLocal(dat=swiss_dat, a0 = 0.005, b0 = 0.005,
                   m0 = 0, s0 = 100, m0l =0,s0l=3,initN = c(0,0),init2L=c(0,0,0))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
