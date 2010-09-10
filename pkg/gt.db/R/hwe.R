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

#---------------------------------------------------------------------

# Efficient, vectorized HWE tests

.hwe.exact <- function(aa, ab, bb, tail)
{
    # use integer constants to avoid implicit casts
    m  <- ab %/% 2L
    xx <- aa + m + 1L
    xy <- ab %% 2L + 1L
    yy <- bb + m + 1L
    n  <- pmin(xx,yy) - 1L

    # precompute lookup table of log gamma values
    # note that xx, xy, yy are offset by 1 for indexing
    LG <- lgamma(1:(max(xx,yy,2L*n)+2))
    v <- log(2)*ab - LG[aa+1L] - LG[ab+1L] - LG[bb+1L]
    log.p <- function(xx, xy, yy, i)
        log(2)*(xy+2L*i-1L) - LG[xx-i] - LG[xy+2L*i] - LG[yy-i]

    fn.lower <- function(v, m, n, ...)
    {
        p <- exp(log.p(...,0:n)-v)
        sum(p[0:m+1L])/sum(p)
    }
    fn.upper <- function(v, m, n, ...)
    {
        p <- exp(log.p(...,0:n)-v)
        sum(p[m:n+1L])/sum(p)
    }
    fn.both <- function(v, m, n, ...)
    {
        p <- exp(log.p(...,0:n)-v)
        sum(p[p<=1])/sum(p)
    }

    mapply(switch(tail, lower=fn.lower, upper=fn.upper, fn.both),
           v, m, n, xx, xy, yy)
}

hwe.test <-
function(aa, ab, bb, test=c('lratio','chisq','exact'),
         tail=c('both','lower','upper'))
{
    test <- match.arg(test)
    tail <- match.arg(tail)

    bad <- (aa<0 | ab<0 | bb<0 | aa+ab+bb==0)
    aa[bad] <- ab[bad] <- bb[bad] <- 0

    if (test == 'exact') {
        pv <- .hwe.exact(aa, ab, bb, tail)
        return(ifelse(bad, NA, pv))
    }

    n <- aa+ab+bb
    if (test == 'lratio') {
        nln <- function(x) x*log(ifelse(x,x,1))
        a <- 2*aa+ab ; b <- 2*bb+ab
        D <- 4*n*aa - a^2
        x <- 2 * (0.25*nln(4*n) - ab*log(2) - nln(a) - nln(b)
                  + nln(aa) + nln(ab) + nln(bb))
    } else if (test == 'chisq') {
        p <- (2*aa+ab)/(2*n)
        D <- (aa/n) - p^2
        x <- n * D^2 / (p*(1-p))^2
    }

    pv <- pchisq(pmax(0,x), df=1, lower.tail=FALSE)
    pv <- switch(tail, both=pv,
                 lower=0.5*ifelse(D<0, pv, 2-pv),
                 upper=0.5*ifelse(D>0, pv, 2-pv))
    ifelse(bad, NA, pv)
}
