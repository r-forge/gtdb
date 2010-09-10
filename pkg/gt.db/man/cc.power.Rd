%
% Copyright (C) 2010, 23andMe, Inc.
%
% Written by David A. Hinds <dhinds@sonic.net>
%
% This is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the license, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>
%
\name{cc.power}
\alias{cc.power}
\title{Power Calculation for Case-Control Association Studies}
\description{
  Compute power, or sample size, for a case-control association study,
  with adjustment for linkage disequilibrium.
}
\usage{
cc.power(p, prevalence, low.risk, rel.risk, odds.ratio,
         N, case.fraction=0.5, alpha=0.05, power,
         pop.controls=FALSE, m, dprime, rsqr,
         method=c('trend','binomial','schork','simulate'), ...)
}
\arguments{
  \item{p}{frequency of the high risk allele.}
  \item{prevalence}{population prevalence of case status.}
  \item{low.risk}{baseline risk for the AA genotype.}
  \item{rel.risk}{either the per-allele (log-additive) relative risk,
    or a vector of two relative risks for AB vs AA, and BB vs AA.}
  \item{odds.ratio}{either the per-allele (log-additive) odds ratio,
    or a vector of two odds ratios for AB vs AA, and BB vs AA.}
  \item{N}{either the total number of cases and controls, or (if length
    = 2) counts of cases and controls.}
  \item{case.fraction}{if \code{N} is missing or gives the total study
    size, the fraction of the total samples that are cases.}
  \item{alpha}{the desired significance level.}
  \item{power}{if specified, the desired study power.}
  \item{pop.controls}{logical: indicates if controls are unverified and
    are assumed to include cases at the population prevalence.}
  \item{m}{if specified, the marker allele frequency.}
  \item{dprime,rsqr}{alternate metrics for specifying linkage
    disequilibrium between the marker and the causal variant.}
  \item{method}{the scoring criterion used to determine power.  See
    details.}
  \item{\dots}{additional method-specific arguments.  See details.}
}
\details{
  This calculates study power for a given sample size (if \code{N} is
  specified), or sample size for a desired power (if \code{power} is
  specified).

  Several methods for calculating power are implemented:
  \describe{
    \item{trend}{power for the Cochran Armitage trend test, based on
      Freidlin \emph{et al.}, 2002.}
    \item{binomial}{power for a two-sample binomial test.}
    \item{delta.p}{power for a test on the allele frequency difference
      between cases and controls.}
    \item{simulate}{simulated power for an arbitrary scoring function.}
  }

  The \code{simulate} method has the following additional arguments:
  \describe{
    \item{score.fn}{the scoring function to use: defaults to
      \code{\link{score.trend}}.}
    \item{tries}{the number of iterations to perform.}
    \item{progress}{logical: indicates whether to show a progress bar.}
    \item{\dots}{additional arguments passed to \code{score.fn}}
  }
}
\value{
  A list with four elements:
  \item{model}{a list as returned by \code{\link{cc.model}}.}
  \item{N}{the study size (numbers of cases and controls).}
  \item{alpha}{the desired false positive rate (probability of
    incorrectly rejecting the null hypothesis when it is true).}
  \item{power}{the study power (probability of rejecting the null
    hypothesis).}
}
\references{
  Freidlin, B., Zheng, G., Li, Z., & Gastwirth, J.L. (2002) Trend tests
  for case-control studies of genetic markers: power, sample size, and
  robustness. \emph{Hum. Hered.} \bold{53}: 146-152.

  Altshuler, D., Daly, M. J., & Lander, E.S. (2008) Genetic mapping in
  human disease. \emph{Science} \bold{332}: 881-888.
}
\seealso{
  \code{\link{cc.model}}, \code{\link{score.trend}}.
}
\examples{
# Fig. 2 from Altshuler et al. (2008)

study.size <- function(p, odds, alpha=1e-8, power=0.9)
    sum(cc.power(p, prevalence=0.05, odds=odds,
                 alpha=alpha, power=power)$N)
odds <- c(10,5,3,2,1.5,1.3,1.2,1.1)
d <- expand.grid(p=c(0.003,0.01,0.03,0.10,0.30), odds=odds)
d$N <- mapply(study.size, d$p, d$odds)
fn <- function(x) -log(x-1)
panel.fn <- function(...)
{
    panel.grid(-1,0)
    ref <- trellis.par.get('reference.line')
    panel.abline(v=fn(unique(d$odds)), col.line=ref$col)
    panel.xyplot(...)
}
xyplot(N~fn(odds), d, groups=p, xlab='Odds ratio', ylab='Sample size',
       scales=list(x=list(at=fn(odds), labels=format(odds)),
                   y=list(log=TRUE, at=c(1e5,3e4,1e4,3e3,1e3,3e2,1e2))),
       type=c('p','l'), panel=panel.fn, xlim=fn(c(10,1.1)), ylim=c(100,1e5))

# sample sizes in Fig. 2, Altshuler et al. (2008)

N <- c(1e5,3e4,1e4,3e3,1e3,3e2,1e2)
pwr.fn <- function(OR, N)
    cc.power(0.1, 0.05, odds.ratio=OR, N=N, alpha=1e-8)$power - 0.90
OR <- sapply(N, function(N) uniroot(pwr.fn, c(1.1,25), N)$root)
equiv.size <- function(alpha, power)
    sapply(OR, function(x) sum(cc.power(0.1, 0.05, odds.ratio=x,
                                        alpha=alpha, power=power)$N))
data.frame(gwas_90=N,
           gwas_50=equiv.size(1e-8,0.5),
           gwas_10=equiv.size(1e-8,0.1),
           rare_90=equiv.size(1e-5,0.9),
           nomn_90=equiv.size(0.05,0.9))
}
\keyword{manip}
