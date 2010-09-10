#
# Copyright (C) 2010, 23andMe, Inc.
#
# Written by David A. Hinds <dhinds@sonic.net>
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the license, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#-----------------------------------------------------------------------
#
# power for tests on quantitative traits
#
#-----------------------------------------------------------------------

.additive.effect <- function(p, H, d)
{
    # note: d here is really d/a
    x = (1-d*(2*p-1))^2 + 2*p*(1-p)*d^2
    sqrt((H/(1-H)) / (2*p*(1-p)*x))
}

.heritability <- function(p, a, d)
{
    # here, we assume Ve = 1
    q  <- 1-p
    Va <- 2*p*q*a^2*(1+d*(q-p))^2
    Vd <- (2*p*q*a*d)^2
    (Va+Vd)/(Va+Vd+1)
}

qtl.model <- function(p, add, dom=0, H)
{
    q <- 1-p
    gt.freq <- c(q^2, 2*p*q, p^2)
    if (missing(add))
        add <- .additive.effect(p, H, dom)
    else
        H <- .heritability(p, add, dom)
    gn <- c('AA','AB','BB')
    list(allele.freq=c(A=q,B=p),
         gt.freq=structure(gt.freq,names=gn),
         add.effect=add, dom.effect=add*dom,
         heritability=H)
}

.qtl.simulate <- function(model, N, rdist)
{
    gt <- .rgt(N, model$gt.freq)
    pt <- rdist(N) + gt*model$add.effect + (gt==1)*model$dom.effect
    data.frame(trait=pt, genotype=gt)
}

.qtl.mpower <- function(model, N, alpha=0.05, power)
{
    if (missing(power)) {
        power <- pnorm(sqrt(N * model$heritability) + qnorm(alpha/2))
    } else {
        N <- (qnorm(power) - qnorm(alpha/2))^2 / model$heritability
        N <- ceiling(N)
    }
    list(N=N, alpha=alpha, power=power)
}

.qtl.spower <- function(model, N, alpha=0.05, score.fn=score.lm,
                        tries=1000, progress=TRUE, rdist=rnorm, ...)
{
    # power for an arbitrary scoring function
    fn <- function(x) {
        if (progress) progress.bar(x, tries)
        data <- .qtl.simulate(model, N, rdist)
        score.fn(trait~genotype, data, ...)$pvalue[1]
    }
    if (progress) progress.bar(0, tries)
    pwr <- mean(sapply(1:tries, fn) < alpha)
    list(N=N, alpha=alpha, power=pwr)
}

qtl.power <- function(p, add, dom=0, H, N, alpha=0.05,
                      method=c('model','simulate'), ...)
{
    model <- qtl.model(p, add, dom, H)
    fn <- switch(match.arg(method),
                 model=.qtl.mpower,
                 simulate=.qtl.spower)
    c(list(model=model), fn(model, N, alpha, ...))
}

#-----------------------------------------------------------------------
#
# power calculations for simple case-control studies
#
#-----------------------------------------------------------------------

cc.model <-
    function(p, prevalence, low.risk, rel.risk,
             odds.ratio, pop.controls=FALSE)
{
    q <- 1-p
    gt.freq <- c(q^2, 2*p*q, p^2)
    if (!missing(rel.risk)) {
        if (length(rel.risk) == 1)
            rel.risk <- rel.risk^(1:2)
        if (missing(low.risk))
            low.risk <- prevalence / sum(gt.freq * c(1, rel.risk))
        odds.ratio <- (1-low.risk)*rel.risk/(1-rel.risk*low.risk)
    } else {
        if (length(odds.ratio) == 1)
            odds.ratio <- odds.ratio^(1:2)
        fn <- function(low.risk) prevalence - low.risk *
            sum(gt.freq * c(1, odds.ratio/(1+(odds.ratio-1)*low.risk)))
        if (missing(low.risk))
            low.risk <- uniroot(fn, c(0,1), tol=0.0001*prevalence)$root
        rel.risk <- odds.ratio/(1+(odds.ratio-1)*low.risk)
    }
    if (missing(prevalence))
        prevalence <- low.risk * sum(gt.freq * c(1,rel.risk))
    prevalence <- c(case=prevalence,
                    control=ifelse(pop.controls,1,1-prevalence))
    gn <- c('AA','AB','BB')
    rn <- c('AB vs AA', 'BB vs AA')
    penetrance <- low.risk * c(1,rel.risk)
    penetrance <-
        rbind(case=structure(penetrance, names=gn),
              control=if (pop.controls) c(1,1,1) else 1-penetrance)
    m <- list(prevalence=prevalence,
              allele.freq=c(A=q,B=p),
              gt.freq=structure(gt.freq,names=gn),
              penetrance=penetrance,
              odds.ratio=structure(odds.ratio,names=rn),
              rel.risk=structure(rel.risk,names=rn))
    m$exp.gt.freq <- with(m, rep(gt.freq,each=2)*penetrance/prevalence)
    m$exp.freq <- m$exp.gt.freq %*% matrix(c(2:0,0:2),3,2) * 0.5
    dimnames(m$exp.freq)[[2]] <- c('A','B')
    m
}

adjust.ld <- function(model, m, dprime, rsqr)
{
    # adjust a cc.model to account for incomplete LD
    stopifnot(missing(dprime)+missing(rsqr) == 1)
    model <- c(list(orig.model=model), model)
    p <- model$allele.freq[2]
    if (!missing(dprime))
        D <- dprime * min(m*(1-p), p*(1-m))
    else
        D <- sqrt(rsqr * p*(1-p)*m*(1-m))

    # adapted from Risch and Teng, 1998
    q <- 1-p
    wd <- D/m; pd <- D/(1-m)
    pp <- p+wd ; pm <- p-pd ; qp <- q+pd ; qm <- q-wd
    mm <- matrix(c(qp^2, 2*pm*qp, pm^2,
                   qm*qp, pp*qp+qm*pm, pp*pm,
                   qm^2, 2*pp*qm, pp^2),3,3)
    model$penetrance <- model$penetrance %*% mm

    model$allele.freq <- c(A=1-m, B=m)
    model$gt.freq <- c(AA=(1-m)^2,AB=2*m*(1-m),BB=m^2)
    model$exp.gt.freq <-
        with(model, rep(gt.freq,each=2)*penetrance/prevalence)
    model$exp.freq <- model$exp.gt.freq %*% matrix(c(2:0,0:2),3,2) * 0.5
    dimnames(model$exp.freq)[[2]] <- c('A','B')
    model
}

.rgt <- function(N, gt.freq)
{
    # generate N random genotypes given AA, AB, BB probabilities
    r <- runif(N)
    ifelse(r < gt.freq[1], 0, ifelse(r > 1-gt.freq[3], 2, 1))
}

.cc.simulate <- function(model, N)
{
    # generate simulated case-control data based on a cc.model
    f <- model$exp.gt.freq
    a <- data.frame(trait=1, genotype=.rgt(N[1], f['case',]))
    b <- data.frame(trait=0, genotype=.rgt(N[2], f['control',]))
    d <- rbind(a,b)
}

.cc.tpower <-
    function(model, N, alpha, scores=0:2)
{
    # power for Cochran Armitage trend test
    # from Freidlin et al., Hum Hered 53: 146-152, 2002
    ssq <- function(x,p) sum(x^2*p)-sum(x*p)^2
    nn <- N[1]+N[2]
    p <- model$exp.gt.freq[1,]
    q <- model$exp.gt.freq[2,]
    u <- N[1]*N[2]/nn^2 * sum(scores*(p-q))
    s = (N[2]*ssq(scores,p) + N[1]*ssq(scores,q)) * N[1]*N[2]/nn^3
    t = (N[1]*ssq(scores,p) + N[2]*ssq(scores,q)) * N[1]*N[2]/nn^3
    za <- -qnorm(alpha/2)
    p <- pnorm((-za*sqrt(t+u^2) - sqrt(nn)*u)/sqrt(s), lower.tail=TRUE) +
         pnorm((za*sqrt(t+u^2) - sqrt(nn)*u)/sqrt(s), lower.tail=FALSE)
    list(N=N, alpha=alpha, power=unname(p))
}

.cc.bpower <- function(model, N, alpha, power)
{
    # power for two-sample binomial test
    require(Hmisc)
    if (missing(power)) {
        power <- bpower(model$exp.freq[1,1], model$exp.freq[2,1],
                        n1=2*N[1], n2=2*N[2], alpha=alpha)
    } else {
        N <- bsamsize(model$exp.freq[1,1], model$exp.freq[2,1],
                      N[1]/(N[1]+N[2]), alpha=alpha, power=power)
        N <- structure(ceiling(N/2), names=c('case','control'))
    }
    list(N=N, alpha=alpha, power=unname(power))
}

.cc.qpower <- function(model, N, alpha, power)
{
    # from Schork et al., 2000
    p <- model$exp.freq[,2]
    pbar <- (N[1]*p[1]+N[2]*p[2])/(N[1]+N[2])
    v <- (1+N[1]/N[2])*pbar*(1-pbar)
    if (missing(power)) {
        power <- pnorm(sqrt(2*N[1]*(p[2]-p[1])^2 / v) + qnorm(alpha/2))
    } else {
        N <- N * (qnorm(power) - qnorm(alpha/2))^2 * v / (2*(p[2]-p[1])^2)
        N <- structure(ceiling(N), names=c('case','control'))
    }
    list(N=N, alpha=alpha, power=unname(power))
}

.cc.spower <- function(model, N, alpha, score.fn=score.trend,
                      tries=1000, progress=TRUE, ...)
{
    # empirical power for an arbitrary scoring function
    fn <- function(x) {
        if (progress) progress.bar(x, tries)
        data <- .cc.simulate(model, N)
        score.fn(trait~genotype, data, ...)$pvalue[1]
    }
    if (progress) progress.bar(0, tries)
    list(N=N, alpha=alpha, power=mean(sapply(1:tries, fn) < alpha))
}

.cc.split <- function(N, case.fraction)
{
    if (length(N) == 1)
        N <- N*c(case.fraction,1-case.fraction)
    stopifnot(length(N) == 2)
    structure(N,names=c('case','control'))
}

.cc.power <-
    function(model, N, case.fraction=0.5, alpha=0.05, power,
             method=c('trend','binomial','schork','simulate'), ...)
{
    fn <- switch(match.arg(method),
                 trend=.cc.tpower, binomial=.cc.bpower,
                 simulate=.cc.spower, schork=.cc.qpower)
    if (!missing(power)) {
        N <- c(1,(1-case.fraction)/case.fraction)
        if ('power' %in% names(formals(fn)))
            ret <- fn(model, N, alpha, power=power, ...)
        else {
            nfn <- function(n) fn(model, n*N, alpha, ...)$power - power
            x <- try(uniroot(nfn, c(1,1000000), tol=0.5), TRUE)
            if (inherits(x, 'try-error'))
                x <- list(root=NA)
            N <- structure(ceiling(N * x$root), names=c('case','control'))
            ret <- list(N=N, alpha=alpha, power=power)
        }
    } else {
        N <- .cc.split(N, case.fraction)
        ret <- fn(model, N, alpha, ...)
    }
    c(list(model=model), ret)
}

cc.power <-
    function(p, prevalence, low.risk, rel.risk, odds.ratio,
             N, case.fraction=0.5, alpha=0.05, power,
             pop.controls=FALSE, m, dprime, rsqr,
             method=c('trend','binomial','schork','simulate'), ...)
{
    model <- cc.model(p, prevalence, low.risk, rel.risk,
                      odds.ratio, pop.controls)
    if (!missing(m))
        model <- adjust.ld(model, m, dprime, rsqr)
    .cc.power(model, N, case.fraction, alpha, power, method, ...)
}

#-----------------------------------------------------------------------
#
# construct a dichotomous trait from a QTL
#
# from Schork et al (2000)
#
#-----------------------------------------------------------------------

.qtl.tail <- function(p, a, d, area, pdist=pnorm)
{
    # given allele frequency, additive, and dominance effects,
    # determine the threshold for given tail area
    mixture <- function(x)
        pdist(x+a)*(1-p)^2 + pdist(x-d)*2*p*(1-p) + pdist(x-a)*p^2 - area
    uniroot(mixture, interval=c(-20,20))$root
}

qtl.to.cc <- function(model, upper.tail, lower.tail, pdist=pnorm)
{
    p <- model$allele.freq[2]
    a <- model$add.effect
    d <- model$dom.effect
    t1 <- .qtl.tail(p, a, d, 1-upper.tail, pdist)
    t2 <- .qtl.tail(p, a, d, lower.tail, pdist)
    model$prevalence <- c(case=upper.tail,control=lower.tail)
    model$penetrance <-
        rbind(case=pdist(t1-c(-a,d,a),lower.tail=FALSE),
              control=pdist(t2-c(-a,d,a),lower.tail=TRUE))
    colnames(model$penetrance) <- c('AA','AB','BB')
    model$exp.gt.freq <-
        with(model, rep(gt.freq,each=2) * penetrance / prevalence)
    model$exp.freq <- model$exp.gt.freq %*% matrix(c(2:0,0:2),3,2) * 0.5
    dimnames(model$exp.freq)[[2]] <- c('A','B')
    model
}

qtl.cc.power <-
    function(p, add, dom=0, H, upper.tail, lower.tail,
             m, dprime, rsqr, ...)
{
    model <- qtl.model(p, add, dom, H)
    model <- qtl.to.cc(model, upper.tail, lower.tail)
    if (!missing(m))
        model <- adjust.ld(model, m, dprime, rsqr)
    .cc.power(model, ...)
}

#-----------------------------------------------------------------------
#
# from Risch and Teng (1998)
#
#-----------------------------------------------------------------------

.adjust.penetrance <- function(f1, f2, p, m, delta)
{
    q <- 1-p
    w <- p*(1-m)/m
    wd <- w*delta; pd <- p*delta
    g2 <- ((f2*(p+wd)^2 + 2*f1*(p+wd)*(q-wd) + (q-wd)^2) /
	   (f2*(p-pd)^2 + 2*f1*(p-pd)*(q+pd) + (q+pd)^2))
    g1 <- (f2*(p+wd)*(p-pd) + (q-wd)*(q+pd) +
	   f1*((p+wd)*(q+pd) + (q-wd)*(p-pd))) /
	   (f2*(p-pd)^2 + 2*f1*(p-pd)*(q+pd) + (q+pd)^2)
    list(f1=g1, f2=g2)
}

.mating.dist <- function(r, p, f1, f2)
{
    # distribution of parental mating types conditional on number
    # of affected sibs and the penetrance of AB and BB genotypes
    gt <- c((1-p)^2,2*p*(1-p),p^2)
    gf <- gt %*% t(gt)
    m <- matrix(c(1, (f1+1)/2, f1,
                  (f1+1)/2, (f2+2*f1+1)/4, (f2+f1)/2,
                  f1, (f2+f1)/2, f2), 3, 3)^r * gf
    m/sum(m)
}

.sib.unr.mu.sd <- function(r, u, p, f1, f2)
{
    m <- .mating.dist(r, p, f1, f2)

    e10 <- f1 / (2*f1+2)
    e11 <- (f2+f1) / (f2+2*f1+1)
    e12 <- (2*f2+f1) / (2*f2+2*f1)
    e <- rbind(c(0.0,e10,0.5), c(e10,e11,e12), c(0.5,e12,1.0))

    t10 <- f1 / (4 * (f1+1)^2)
    t11 <- (f2*f1+2*f2+f1) / (2*(f2+2*f1+1)^2)
    t12 <- f2*f1 / (4*(f1+f2)^2)
    t <- rbind(c(0.0,t10,0.0), c(t10,t11,t12), c(0.0,t12,0.0))

    mu <- sum(e*m) - p
    sd <- sqrt(sum(t*m/r + e^2*m) - (mu+p)^2 + p*(1-p) / (2*u))
    list(mu=mu, sd=sd)
}

.sib.unr.power <- function(model, N, r=1, u=1, alpha=0.05, power)
{
    # r affected sibs, u unaffected controls
    f <- model$penetrance[1,]
    p <- unname(model$allele.freq[2])
    m <- .sib.unr.mu.sd(r, u, p, f[2]/f[1], f[3]/f[1])
    pc <- p + 2*r*m$mu/(2*r+(r+1)*u)
    a <- qnorm(1-alpha/2)*sqrt((r*u+2*r+u)*pc*(1-pc) / (4*r*u))
    if (missing(power)) {
        power <- pnorm((sqrt(N)*m$mu - a)/m$sd)
    } else {
        N <- ceiling(((m$sd*qnorm(power) + a) / m$mu)^2)
        N <- c(case=N*r, control=N*u, families=N)
    }
    list(N=N, alpha=alpha, power=power)
}

.sib.sib.mu.sd <- function(r, s, p, f1, f2)
{
    m <- .mating.dist(r, p, f1, f2)

    p10 <- 0.25 * (f1-1) / (f1+1)
    p11 <- 0.5 * (f2-1) / (1+2*f1+f2)
    p12 <- 0.25 * (f2-f1) / (f1+f2)
    p <- rbind(c(0.0,p10,0.0), c(p10,p11,p12), c(0.0,p12,0.0))

    w10 <- 0.25 * f1 / (f1+1)^2
    w11 <- 0.5 * (2*f2+f2*f1+f1) / (1+2*f1+f2)^2
    w12 <- 0.25 * f2*f1 / (f2+f1)^2
    w <- rbind(c(0.0,w10,0.0), c(w10,w11,w12), c(0.0,w12,0.0))

    mu <- sum(p*m)
    sd <- sqrt(sum(w*m/r + p^2*m) - mu^2 +
               sum(outer(c(0,1,0),c(0,1,0),'+') * m/(16*s)))
    pc <- sum(outer(0:2,0:2,'+') * m/4) + (r/(r+s))*mu
    list(mu=mu,sd=sd,pc=pc)
}

.sib.sib.power <- function(model, N, r=1, s=1, alpha=0.05, power)
{
    # r affected sibs, s unaffected sib controls
    f <- model$penetrance[1,]
    p <- unname(model$allele.freq[2])
    m <- .sib.sib.mu.sd(r, s, p, f[2]/f[1], f[3]/f[1])
    a <- qnorm(1-alpha/2)*sqrt((r+s)*m$pc*(1-m$pc) / (4*r*s))
    if (missing(power)) {
        power <- pnorm((sqrt(N)*m$mu - a)/m$sd)
    } else {
        N <- ceiling(((m$sd*qnorm(power) + a) / m$mu)^2)
        N <- c(case=N*r, control=N*s, families=N)
    }
    list(N=N, alpha=alpha, power=power)
}
