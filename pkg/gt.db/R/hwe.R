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
    nr <- pmin(aa,bb)*2+ab
    nn <- aa+ab+bb
    xy <- as.integer(nr*(2*nn-nr) / (2*nn))
    xy <- xy + ((nr+xy) %% 2)
    xx <- (nr-xy) %/% 2
    yy <- (nn-xy-xx)

    p.fn <- function(xx, xy, yy)
    {
        if (is.na(xx+xy+yy)) return(NA)
        i <- 1:(xy %/% 2)
        lt <- cumprod((xy-2*i+2)*(xy-2*i+1)/(4*(xx+i)*(yy+i)))
        i <- 1:min(xx,yy)
        ut <- cumprod(4*(xx-i+1)*(yy-i+1)/((xy+2*i)*(xy+2*i-1)))
        c(rev(lt),1,ut)
    }

    fn.lower <- function(brk, p) sum(p[1:brk]) / sum(p)
    fn.upper <- function(brk, p) sum(p[brk:length(p)]) / sum(p)
    fn.both  <- function(brk, p) sum(p[p <= p[brk]]) / sum(p)
    fn <- switch(tail, lower=fn.lower, upper=fn.upper, fn.both)
    mapply(function(brk,xx,xy,yy) fn(brk, p.fn(xx,xy,yy)),
           (ab %/% 2)+1, xx, xy, yy)
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
