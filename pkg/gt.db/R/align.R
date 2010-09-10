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

nw.align <- function(gt1, gt2, gap=-1, pm=1, mm=0, ends=FALSE)
{
    gt1 <- gsub('[^ACGT]','N',toupper(gt1))
    gt2 <- gsub('[^ACGT]','n',toupper(gt2))
    if (length(ends)==1) ends <- rep(ends,4)
    x <- .Call('do_nw_align', gt1, gt2, c(gap,pm,mm),
               ends, PACKAGE='gt.db')
    p <- nsubstr(x[[3]][2],'|') / nchar(c(gt1,gt2))
    list(score=x[[1]], pct.matched=100*p,
         ends=x[[2]], alignment=x[[3]])
}

nw.score <- function(gt1, gt2, gap=-1, pm=1, mm=0, ends=FALSE)
{
    gt1 <- gsub('[^ACGT]','N',toupper(gt1))
    gt2 <- gsub('[^ACGT]','n',toupper(gt2))
    if (length(ends)==1) ends <- rep(ends,4)
    param <- c(gap, pm, mm)
    mapply(.Call, name='do_nw_score', gt1=gt1, gt2=gt2,
           MoreArgs=list(param=param, ends=ends, PACKAGE='gt.db'),
           USE.NAMES=FALSE)
}

nw.orient.assay <- function(gt1, gt2, delta=1)
{
    l1 <- sub('(.*)_.*', '\\1', gt1)
    r1 <- sub('.*_(.*)', '\\1', gt1)
    l2 <- sub('(.*)_.*', '\\1', gt2)
    r2 <- sub('.*_(.*)', '\\1', gt2)
    rhs <- c(FALSE,TRUE,FALSE,TRUE)
    fwd <- (nw.score(l1,l2,ends=!rhs) + nw.score(r1,r2,ends=rhs))
    rev <- (nw.score(l1,revcomp(r2,TRUE),ends=!rhs) +
            nw.score(r1,revcomp(l2,TRUE),ends=rhs))
    f <- ifelse(fwd>rev+delta, '+', ifelse(rev>fwd+delta, '-', NA))
    factor(f, levels=c('+','-'))
}

print.align <- function(x, ..., width=50, pad=10)
{
    do.pad <- function(s)
    {
        n <- seq(1, nchar(s), pad)
        paste(mapply(substr, s, n, n+pad-1), collapse=' ')
    }
    wrap <- function(s)
    {
        n <- seq(1, nchar(s), width)
        mapply(substr, s, n, n+width-1)
    }
    a <- lapply(x$alignment, wrap)
    n1 <- cumsum(c(x$ends[1],nchar(a[[1]])-nsubstr(a[[1]], '-')))
    n2 <- cumsum(c(x$ends[3],nchar(a[[3]])-nsubstr(a[[3]], '-')))
    s1 <- sprintf("%05d ", n1)
    s2 <- sprintf("%05d ", n2)
    e1 <- sprintf(" %05d\n", n1[-1]-1)
    e2 <- sprintf(" %05d\n", n2[-1]-1)
    for (i in 1:length(a[[1]])) {
        cat(s1[i], do.pad(a[[1]][i]), e1[i], sep='')
        cat('      ', do.pad(a[[2]][i]), '\n', sep='')
        cat(s2[i], do.pad(a[[3]][i]), e2[i], sep='')
        if (i < length(a[[1]])) cat('\n')
    }
}
