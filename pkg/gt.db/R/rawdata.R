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

#---------------------------------------------------------------------

.unpack.signal <- function(data)
{
    nr <- length(data)/4
    i <- readBin(data, what='int', n=2*nr, size=2,
                 signed=FALSE, endian='little')
    i <- na.if(i, 65535)
    dim(i) <- c(2,nr)
    data.frame(signal.a=i[1,], signal.b=i[2,])
}

.pack.signal <- function(data)
{
    a <- writeBin(if.na(data$signal.a,65535),
                  raw(), size=2, endian='little')
    b <- writeBin(if.na(data$signal.b,65535),
                  raw(), size=2, endian='little')
    dim(a) <- dim(b) <- c(2,nrow(data))
    as.vector(rbind(a,b))
}

#---------------------------------------------------------------------

.unpack.seqread <- function(data)
{
    i <- na.if(as.integer(data), 255)
    dim(i) <- c(4,length(i)/4)
    data.frame(fwd.a=i[1,], rev.a=i[2,], fwd.b=i[3,], rev.b=i[4,])
}

.pack.seqread <- function(data)
{
    as.vector(rbind(as.raw(if.na(data$fwd.a,255)),
                    as.raw(if.na(data$rev.a,255)),
                    as.raw(if.na(data$fwd.b,255)),
                    as.raw(if.na(data$rev.b,255))))
}

#---------------------------------------------------------------------

.unpack.chpdata <- function(data)
{
    nr <- length(data)/5
    dim(data) <- c(5,nr)
    i <- readBin(data[1:2,], what='int', n=nr, size=2,
                 signed=TRUE, endian='little')
    j <- readBin(data[3:4,], what='int', n=nr, size=2,
                 signed=FALSE, endian='little')
    k <- factor(as.integer(data[5,]), levels=1:4,
                labels=c('a','h','b','n'))
    data.frame(log.ratio=i/256, strength=j/256, forced.call=k)
}

.pack.chpdata <- function(data)
{
    x <- writeBin(as.integer(256*data$log.ratio),
                  raw(), size=2, endian='little')
    y <- writeBin(as.integer(256*data$strength),
                  raw(), size=2, endian='little')
    z <- as.raw(data$forced.call)
    dim(x) <- dim(y) <- c(2,nrow(data))
    as.vector(rbind(x,y,z))
}

#---------------------------------------------------------------------

unpack.raw.data <- function(data, raw.layout)
{
    if (raw.layout == 'signal') {
        .unpack.signal(data)
    } else if (raw.layout == 'seqread') {
        .unpack.seqread(data)
    } else if (raw.layout == 'chpdata') {
        .unpack.chpdata(data)
    } else {
        stop('unknown raw data layout')
    }
}

pack.raw.data <- function(data, raw.layout)
{
    if (raw.layout == 'signal') {
        .pack.signal(data)
    } else if (raw.layout == 'seqread') {
        .pack.seqread(data)
    } else if (raw.layout == 'chpdata') {
        .pack.chpdata(data)
    } else {
        stop('unknown raw data layout')
    }
}
