#
# Jonckheere-Terpstra test implementation
#
# As posted by Chris Andrews <candrews@buffalo.edu>
# on the R-Help mailing list, Fri 30 Jun 2006
# http://tolstoy.newcastle.edu.au/R/help/06/06/30112.html
#

jt.test <- function(x, y,
   alternative = c("two.sided", "decreasing", "increasing"),
   asymp=TRUE, correct=FALSE, perm=0, na.action=c("omit", "fail"),
   permgraph=FALSE, permreps=FALSE)
{

# Jonckheere-Terpstra test
#  x = response
#  y = group membership
#  alternative hypothesis
#  asymptotic = asymptotic formula for variance or don't bother
#  correct = continuity correction
#  perm = number of repetitions for a permutation test
#  na.action = what to do with missing data
#  permgraph = draw a histogram of the permutations?
#  permreps = return the permutations?

   alternative <- match.arg(alternative)
   complete <- complete.cases(x,y)
   nn <- sum(complete)
   #cat(nn, "complete cases.\n")
   if (!all(complete) && (na.action=="fail")) {
     stop("option na.action is 'fail' and",
       sum(!complete), "cases are not complete.\n")
   }
   response <- x[complete]
   groupvec <- y[complete, drop=TRUE]

   if(!is.ordered(groupvec)) {
     groupvec <- as.ordered(groupvec)
     warning("y was not an ordered factor.  Redefined to be one.\n")
   }

   lev <- levels(groupvec)
   nlev <- length(lev)
   if(nlev <= 2) {
     if (nlev < 2) {
       stop("Not enough groups (k =", nlev, ") for Jonckheere.\n")
     } else {
       #warning("Two groups.  Could use wilcox.test.\n")
     }
   }

   computestat <- function(response, groupvec, n.groups) {
     pairwise <- function(x, y) {
       lx <- length(x)
       lx*(length(y) + (1+lx)/2) - sum(rank(c(x,y))[seq(along=x)])
     }
     H<-0
     for (i in seq(n.groups-1))
       for (j in seq(i+1, n.groups))
         H <- H+pairwise(response[groupvec == lev[i]],
                         response[groupvec == lev[j]])
     return(H)
   }

   retval <- list(
     statistic=computestat(response, groupvec, nlev),
     alternative = paste(alternative,
       paste(levels(groupvec),collapse=switch(alternative,
         two.sided = ", ", decreasing= " > ", increasing = " < ")), sep=": "),
     method = "Jonckheere-Terpstra test",
     data.name = paste(deparse(substitute(x)), "by", deparse(substitute(y))))

   if (asymp) {
     ns <- tabulate(as.numeric(groupvec))
     retval$EH <- (nn * nn - sum(ns * ns))/4
     retval$VH <- if (!any(duplicated(response))) {
       (nn * nn * (2 * nn + 3) - sum(ns * ns * (2 * ns + 3)))/72
     } else {
       ds <- as.vector(table(response))
       ((nn*(nn-1)*(2*nn+5) - sum(ns*(ns-1)*(2*ns+5)) - 
       sum(ds*(ds-1)*(2*ds+5)))/72 +
       sum(ns*(ns-1)*(ns-2))*sum(ds*(ds-1)*(ds-2))/(36*nn*(nn-1)*(nn-2)) +
       sum(ns*(ns-1))*sum(ds*(ds-1))/(8*nn*(nn - 1)))
     }
     retval$N <- nn
     pp <- if (!correct) {
       pnorm(retval$statistic, retval$EH, sqrt(retval$VH))
     } else if (retval$statistic >= retval$EH + 1) {
       pnorm(retval$statistic - 1, retval$EH, sqrt(retval$VH))
     } else if (retval$statistic <= retval$EH - 1) {
       pnorm(retval$statistic + 1, retval$EH, sqrt(retval$VH))
     } else {
       .5
     }
     retval$p.value <- switch(alternative,
       two.sided = 2*min(pp,1-pp), decreasing = pp, increasing = 1-pp)
   }
   if (perm>0) {
     reps <- numeric(perm)
     for (i in seq(perm)) {
       reps[i] <- computestat(response, sample(groupvec), nlev)
     }
     retval$EH.perm <- mean(reps)
     retval$VH.perm <- var(reps)
     ppp <- if (!correct) {
       pnorm(retval$statistic, retval$EH.perm, sqrt(retval$VH.perm))
     } else if (retval$statistic >= retval$EH.perm + 1) {
       pnorm(retval$statistic - 1, retval$EH.perm, sqrt(retval$VH.perm))
     } else if (retval$statistic <= retval$EH.perm - 1) {
       pnorm(retval$statistic + 1, retval$EH.perm, sqrt(retval$VH.perm))
     } else {
       .5
     }
     retval$p.value.perm <- switch(alternative,
       two.sided = 2*min(ppp,1-ppp), decreasing = ppp, increasing = 1-ppp)
     rr <- rank(c(retval$statistic, reps))[1]/(perm+2)
     retval$p.value.rank <- switch(alternative,
       two.sided = 2*min(rr,1-rr), decreasing = rr, increasing = 1-rr)
     if (permreps)
       retval$reps <- reps
     if (permgraph) {
       hist(reps, xlab="Jonckheere Statistic", prob=TRUE,
         main="Histogram of Permutations", sub=paste(perm, "permutations"))
       abline(v=retval$statistic, col="red")
       points((-3:3)*sqrt(retval$VH.perm)+retval$EH.perm, rep(0,7), col="blue",
         pch=as.character(c(3:1,0:3)))
     }
   }
   class(retval) <- "htest"
   names(retval$statistic) <- "J"
   return(retval)

# returns list (of class htest) with the following components
#
# Always:
#  statistic = value of J
#  alternative = same as input
#  method = "Jonckheere"
#  data.name = source of data based on call
#
# Most of the time (i.e., when asymp=TRUE):
#  EH = expected value based on sample size
#  VH = variance (adjusted for ties if necessary)
#  p.value = asymptotic p value
#
# Sometimes (i.e., when perm>0):
#  EH.perm = expected value based on permutations
#  VH.perm = variance based on permutations
#  p.value.perm = p value based on normal approximation to permutations
#  p.value.rank = p value based on rank within permutation
#  perms = vector of permutation.

}
