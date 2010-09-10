#
# Copyright (C) 2009, Perlegen Sciences, Inc.
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

.ld.3x3.iterate <- function(m, measure, epsilon, max.it)
{
    n.aa <- 2*m[,1,1] + m[,1,2] + m[,2,1]
    n.ab <- 2*m[,1,3] + m[,1,2] + m[,2,3]
    n.ba <- 2*m[,3,1] + m[,2,1] + m[,3,2]
    n.bb <- 2*m[,3,3] + m[,2,3] + m[,3,2]
    n.hh <- m[,2,2]

    # prepare for EM iteration
    j <- 1:length(n.hh)
    x.aa <- n.aa
    x.ab <- n.ab + n.hh
    x.ba <- n.ba + n.hh
    x.bb <- n.bb
    x.hh <- n.hh
    x.aa.bb <- n.hh * 0.5
    aa <- ab <- ba <- bb <- vector()

    # We only iterate instances that have not converged
    for (i in 1:max.it) {
        j.aa <- x.aa + x.aa.bb ; j.ab <- x.ab - x.aa.bb
        j.ba <- x.ba - x.aa.bb ; j.bb <- x.bb + x.aa.bb
        a <- j.aa*j.bb; b <- j.ab*j.ba
        b <- x.hh * a / (a+b)
        if ((i == max.it) || (i %% 10 == 0)) {
            aa[j] <- j.aa; ab[j] <- j.ab; ba[j] <- j.ba; bb[j] <- j.bb
            w <- which(abs(x.aa.bb-b)>epsilon)
            if (!length(w)) break
            j <- j[w] ; b <- b[w] ; x.hh <- x.hh[w]
            x.aa <- x.aa[w] ; x.ab <- x.ab[w]
            x.ba <- x.ba[w] ; x.bb <- x.bb[w]
        }
        x.aa.bb <- b
    }

    if (!length(w)) j <- integer(0)
    if (measure=='failed') {
        return (1:length(aa) %in% j)
    }

    nn <- aa+ab+ba+bb
    na_ <- aa+ab ; nb_ <- ba+bb
    n_a <- aa+ba ; n_b <- ab+bb
    D <- aa*bb - ab*ba

    if (measure=='dprime') {
        r <- (D / ifelse(D>0, pmin(na_*n_b, n_a*nb_),
                             -pmin(na_*n_a, nb_*n_b)))
        return(structure(r, failed=length(j)))
    }

    var <- na_ * nb_ * n_a * n_b
    zlog10 <- function(x) log10(ifelse(x<=0,1,x))
    loglik <- function(D)
       (n.aa*zlog10(na_*n_a+D) + n.ab*zlog10(na_*n_b-D) +
        n.ba*zlog10(nb_*n_a-D) + n.bb*zlog10(nb_*n_b+D) +
        n.hh*zlog10(2*var + D*(na_-nb_)*(n_a-n_b) + 2*D^2))

#    For calculating likelihood profiles; but what do they mean?
#    {
#        Dmax <-  pmin(na_*n_b, n_a*nb_)
#        Dmin <- -pmin(na_*n_a, nb_*n_b)
#        lo <- sapply(seq(1,0,-step), function(x) loglik(Dmin*x))
#        hi <- sapply(seq(0,1, step)[-1], function(x) loglik(Dmax*x))
#        lr <- cbind(lo,hi) - loglik(D)
#    }

    if (measure=='rsqr') {
        r <- (D^2 / var)
    } else if (measure=='delta') {
        r <- (D / sqrt(var))
    } else if (measure=='lod') {
        r <- loglik(D) - loglik(0)
    } else if (measure=='pvalue') {
        lr <- log(10) * (loglik(D) - loglik(0))
        r <- pchisq(2*lr, df=1, lower.tail=FALSE)
    } else {
        r <- D * NA
    }
    structure(r, failed=length(j))
}

.ld.3x3.exact <-
function (m, measure, epsilon)
{
    n.aa <- 2*m[,1,1] + m[,1,2] + m[,2,1]
    n.ab <- 2*m[,1,3] + m[,1,2] + m[,2,3]
    n.ba <- 2*m[,3,1] + m[,2,1] + m[,3,2]
    n.bb <- 2*m[,3,3] + m[,2,3] + m[,3,2]
    n.hh <- m[,2,2]

    # we solve the cubic equation using counts, not frequencies,
    # because that ends up simplifying some intermediate terms
    nn <- 2*rowSums(m)
    p <- (n.aa+n.ab+n.hh)
    q <- (n.aa+n.ba+n.hh)
    k <- nn - 2*p - 2*q
    b <- k - 2*n.aa - n.hh
    c <- p*q - n.aa*k - n.hh*(nn-p-q)
    d <- -n.aa*p*q

    # indexing a data matrix is faster than a data frame
    v <- cbind(x=-b/6, del2=(b^2-6*c)/36)
    v <- cbind(v, h2=v[,'del2']^3,
               y=(-4*v[,'x']^3 + c*v[,'x'] + d)/4)

    s <- v[,'y']^2 - v[,'h2']
    n11 <- matrix(NA,nrow(v),3)
    if (any(s>0, na.rm=TRUE)) {
        w <- which(s>0)
        t <- sqrt(s[w])
        cubrt <- function(x) sign(x)*(abs(x)^(1/3))
        n11[w,1] <- v[w,'x'] + cubrt(t-v[w,'y']) - cubrt(t+v[w,'y'])
    }
    if (any(s==0, na.rm=TRUE)) {
        w <- which(s==0)
        n11[w,] <- v[w,'x']+outer(v[w,'y']^(1/3), c(1,-2,NA))
    }
    if (any(s<0, na.rm=TRUE)) {
        w <- which(s<0)
        vw <- v[w,,drop=FALSE]
        theta <- outer(acos(-vw[,'y']/sqrt(vw[,'h2'])),c(0,2,4)*pi,'+')/3
        n11[w,] <- vw[,'x']+2*sqrt(vw[,'del2'])*cos(theta)
    }

    n12 <- p - n11
    n21 <- q - n11
    n22 <- nn - n11 - n12 - n21

    # remove biologically impossible roots
    bad <- pmin(n11,n12,n21,n22,na.rm=TRUE) < -epsilon*nn
    n11[bad] <- NA

    # compute disequilibrium for all remaining roots
    Dx <- n11*n22 - n12*n21
    nr <- rowSums(!is.na(Dx))

    # organize counts into a data matrix
    n <- cbind(aa=n.aa, ab=n.ab, ba=n.ba, bb=n.bb, hh=n.hh,
               'a_'=p, 'b_'=(nn-p), '_a'=q, '_b'=(nn-q))
    n <- cbind(n, var=n[,'a_']*n[,'b_']*n[,'_a']*n[,'_b'])

    zlog10 <- function(x) log10(ifelse(x<=0,1,x))
    loglik <- function(n,D)
       (n[,'aa']*zlog10(n[,'a_']*n[,'_a']+D) +
        n[,'ab']*zlog10(n[,'a_']*n[,'_b']-D) +
        n[,'ba']*zlog10(n[,'b_']*n[,'_a']-D) +
        n[,'bb']*zlog10(n[,'b_']*n[,'_b']+D) +
        n[,'hh']*zlog10(2*n[,'var'] + 2*D^2 +
                        D*(n[,'a_']-n[,'b_'])*(n[,'_a']-n[,'_b'])))

    # choose most-likely root, if we have more than one
    D <- ifelse(nr==1, rowSums(Dx, na.rm=TRUE), NA)
    if (any(nr>1)) {
        w <- which(nr>1)
        LL <- if.na(loglik(n[w,,drop=FALSE],Dx[w,,drop=FALSE]),0)
        D[w] <- ifelse(LL[,1]>LL[,2],
                       ifelse(LL[,1]>LL[,3], Dx[w,1], Dx[w,3]),
                       ifelse(LL[,2]>LL[,3], Dx[w,2], Dx[w,3]))
    }

    if (measure=='dprime') {
        (D / ifelse(D>0, pmin(n[,'a_']*n[,'_b'], n[,'_a']*n[,'b_']),
                        -pmin(n[,'a_']*n[,'_a'], n[,'b_']*n[,'_b'])))
    } else if (measure=='rsqr') {
        (D^2 / n[,'var'])
    } else if (measure=='delta') {
        (D / sqrt(n[,'var']))
    } else if (measure=='lod') {
        loglik(n,D) - loglik(n,0)
    } else if (measure=='pvalue') {
        lr <- log(10) * (loglik(n,D) - loglik(n,0))
        pchisq(2*lr, df=1, lower.tail=FALSE)
    } else {
        D * NA
    }
}

.ld.inner <-
function(g1, g2, measure, method, epsilon, max.it)
{
    s <- unique(c(as.character(g1$ploidy),as.character(g2$ploidy)))
    if (any(is.na(s))) {
        warning('assays with unknown ploidy will be assumed autosomal')
        s <- s[!is.na(s)]
    }
    if (length(s) > 1)
        stop('inconsistent ploidy information')
    if (measure == 'none')
        return(rep(NA,max(nrow(g1),nrow(g2))))
    if (s == 'A') {
        m <- ch.table(g1$genotype, g2$genotype, c('a','h','b'))
    } else if (s == 'X') {
        gf <- mask.gt.data(g1, gender(g1)=='F')
        gm <- mask.gt.data(g1, gender(g1)=='M')
        m <- (ch.table(gf$genotype, g2$genotype, c('a','h','b')) +
              0.5*ch.table(gm$genotype, g2$genotype, c('a','_','b')))
    } else if (s %in% c('M','Y')) {
        m <- 0.5 * ch.table(g1$genotype, g2$genotype, c('a','_','b'))
    } else {
        stop('unknown ploidy')
    }
    if (method == 'iterate')
        .ld.3x3.iterate(m, measure, epsilon, max.it)
    else
        .ld.3x3.exact(m, measure, epsilon)
}

.ld.outer <-
function(g1, g2, measure, method, epsilon, max.it)
{
    fn <- function(x,d1,d2)
        .ld.inner(d1[x,], d2, measure, method, epsilon, max.it)
    if (nrow(g1) < nrow(g2))
        l <- lapply(1:nrow(g1), fn, g1, g2)
    else
        l <- lapply(1:nrow(g2), fn, g2, g1)
    r <- do.call('cbind', l)
    if (measure != 'failed') {
        x <- unlist(sapply(l, attr, 'failed'))
        if (!is.null(x))
            attr(r,'failed') <- sum(x)
    }
    if (nrow(g1) < nrow(g2)) t(r) else r
}

ld.gt.data <-
function(g1, g2=g1, outer=TRUE,
         measure=c('rsqr','dprime','delta','pvalue','lod','none','failed'),
         method=c('iterate','exact'), epsilon=1e-6, max.it=100)
{
    measure <- match.arg(measure)
    method <- match.arg(method)
    if (outer)
        r <- .ld.outer(g1, g2, measure, method, epsilon, max.it)
    else
        r <- .ld.inner(g1, g2, measure, method, epsilon, max.it)
    if (measure != 'failed') {
        x <- attr(r,'failed')
        if (!is.null(x) && x) warning(x, ' cases did not converge')
    }
    r
}

.panel.ldplot <-
    function(x, y, z, at=pretty(z), rug=FALSE, ..., draw=TRUE,
             col.regions=regions$col, alpha.regions=regions$alpha)
{
    regions <- trellis.par.get("regions")
    x <- as.numeric(x)
    y <- as.numeric(y)
    zcol <- level.colors(z, at, col.regions, colors = TRUE)

    ux <- sort(unique(x[!is.na(x)]))
    bx <- if (length(ux) > 1)
        c(3 * ux[1] - ux[2], ux[-length(ux)] + ux[-1], 3 *
          ux[length(ux)] - ux[length(ux) - 1])/2
    else ux + c(-1,1)
    lx <- 0.5 * diff(bx)
    cx <- (bx[-1] + bx[-length(bx)])/2

    if (is.list(rug))
        do.call('panel.rug', c(list(x=ux),rug))
    else if (rug)
        panel.rug(ux)

    idx <- match(x, ux)
    idy <- match(y, ux)
    x0 <- (cx[idx]+cx[idy])/2 - (lx[idx]+lx[idy])/2
    xm <- rbind(x0, x0+lx[idy], x0+lx[idx]+lx[idy], x0+lx[idx])
    y0 <- (max(cx) - (cx[idx]-cx[idy]) + (lx[idx]-lx[idy]))/2
    ym <- rbind(y0, y0+lx[idy], y0+lx[idy]-lx[idx], y0-lx[idx])*2
    gp <- gpar(fill=zcol, lwd=1e-5, col='transparent', alpha=alpha.regions)
    grid.polygon(x=as.vector(xm), y=as.vector(ym),
                 id.lengths=rep(4,length(x0)), gp=gp,
                 default.units = "native", draw=draw)
}

ld.plot <-
function(gt.data, col=gray(seq(1,0,-0.01)), measure='rsqr',
         rotate=FALSE, equal=TRUE, colorkey=list(height=0.5),
         scales, ...)
{
    gt.data <- gt.data[order(gt.data$scaffold,gt.data$position),]
    ld <- ld.gt.data(gt.data, measure=measure[1])
    if (length(measure)==2) {
        dn <- ld.gt.data(gt.data, measure=measure[2])
        ld[lower.tri(ld)] <- dn[lower.tri(dn)]
    }
    if (missing(scales)) {
        scales <- if (equal) list(alternating=0) else list()
        if (rotate) scales <- list(y=list(draw=FALSE))
    }
    if (rotate) {
        panel.fn <- .panel.ldplot
        aspect <- 0.5
    } else {
        panel.fn <- panel.levelplot
        aspect <- 1.0
    }
    if (equal) {
        levelplot(ld, aspect=aspect,
                  col.regions=col, colorkey=colorkey,
                  xlab=NULL, ylab=NULL, scales=scales,
                  panel=panel.fn, ...)
    } else {
        dimnames(ld) <- list(gt.data$position)[c(1,1)]
        xd <- as.data.frame.table(ld, stringsAsFactors=FALSE)
        xd[,1] <- as.numeric(xd[,1])
        xd[,2] <- as.numeric(xd[,2])
        levelplot(Freq~Var1*Var2, xd, aspect=aspect,
                  col.regions=col, colorkey=colorkey,
                  xlab=NULL, ylab=NULL, scales=scales,
                  panel=panel.fn, ...)
    }
}

ld.prune <-
function(gt.data, min.maf=0.01, max.rsqr=0.2, span=20, subsets=TRUE)
{
    for (s in subsets) {
        af <- summary.gt.data(gt.data, s)$freq.a
        gt.data <- subset(gt.data, af >= min.maf & 1-af >= min.maf)
    }
    gt.data <- gt.data[order(gt.data$position),]
    gm <- lapply(subsets, function(s) mask.gt.data(gt.data, s, TRUE))

    nr <- n <- nrow(gt.data)
    k <- 1:n
    progress.bar(0, nr)
    while (n > 1) {
        w <- k[max(1,n-span):(n-1)]
        fn <- function(g) which(ld.gt.data(g[k[n],], g[w,]) > max.rsqr)
        # hide warnings about convergence
        r <- suppressWarnings(lapply(gm, fn))
        w <- n - length(w) - 1 + unique(do.call('c', r))
        if (length(w)) {
            k <- k[-w]
            n <- min(length(k), max(w)+span-length(w))
        } else {
            if (n %% 10 == 0) progress.bar(nr-k[n]+1, nr)
            n <- n - 1
        }
    }
    progress.bar(nr, nr)
    gt.data[k,]
}
