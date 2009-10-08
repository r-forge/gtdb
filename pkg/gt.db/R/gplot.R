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

gplot <-
function(formula, data, aggr.fn=max, rescale=FALSE, binsz=1e6,
         subset=TRUE, col.regions=rev(heat.colors(100)[10:100]),
         scales=list(x=list(at=seq(0,250,20), draw=TRUE)),
         shrink=list(x=1,y=0.75), colorkey=list(height=0.25),
         xlab='Position, Mb', ylab='Chromosome', ...)
{
    subset <- eval(substitute(subset), data, environment(formula))
    data <- data[subset,]
    z <- latticeParseFormula(formula, data)$right
    bin <- list(factor(data$scaffold), trunc(data$position/binsz))
    m <- tapply(z, bin, aggr.fn)
    if (rescale) m <- m/mean(m,na.rm=TRUE)
    colnames(m) <- (as.numeric(colnames(m))+0.5) * binsz / 1e6

    chr <- .sort.levels(rownames(m), decreasing=TRUE)
    pos <- rep(as.numeric(colnames(m)), each=nrow(m))
    d <- data.frame(pos=pos, chr=chr, z=as.vector(m))
    p <- function(...) {
        panel.abline(v=scales$x$at,col='gray');
        panel.levelplot(...)
    }

    levelplot(z~pos+chr, d, scales=scales, aspect='fill',
              col.regions=col.regions, colorkey=colorkey,
              shrink=shrink, xlab=xlab, ylab=ylab, panel=p, ...)
}
