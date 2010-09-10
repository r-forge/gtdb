#
# Copyright (C) 2010, 23andMe, Inc.
#
# Written by David A. Hinds <dhinds@23andMe.com>
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

.unpack.dosage <- function(raw)
{
    i <- readBin(raw, what='int', n=length(raw)/2, size=2,
                 signed=FALSE, endian='little')
    data.frame(dosage=na.if(i,65535)/1000)
}

.pack.dosage <- function(data)
{
    writeBin(as.integer(if.na(round(1000*data$dosage),65535)),
             raw(), size=2, endian='little')
}

#---------------------------------------------------------------------

.unpack.signal <- function(raw)
{
    nr <- length(raw)/4
    i <- readBin(raw, what='int', n=2*nr, size=2,
                 signed=FALSE, endian='little')
    i <- na.if(i, 65535)
    dim(i) <- c(2,nr)
    data.frame(signal.a=i[1,], signal.b=i[2,])
}

.pack.signal <- function(data)
{
    a <- writeBin(as.integer(if.na(data$signal.a,65535)),
                  raw(), size=2, endian='little')
    b <- writeBin(as.integer(if.na(data$signal.b,65535)),
                  raw(), size=2, endian='little')
    dim(a) <- dim(b) <- c(2,nrow(data))
    as.vector(rbind(a,b))
}

#---------------------------------------------------------------------

.unpack.seqread <- function(raw)
{
    i <- na.if(as.integer(raw), 255)
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

.unpack.chpdata <- function(raw)
{
    nr <- length(raw)/5
    dim(raw) <- c(5,nr)
    i <- readBin(raw[1:2,], what='int', n=nr, size=2,
                 signed=TRUE, endian='little')
    j <- readBin(raw[3:4,], what='int', n=nr, size=2,
                 signed=FALSE, endian='little')
    k <- factor(as.integer(raw[5,]), levels=1:4,
                labels=c('a','h','b','n'))
    data.frame(log.ratio=i/256, strength=j/256, forced.call=k)
}

.pack.chpdata <- function(data)
{
    x <- writeBin(as.integer(round(256*data$log.ratio)),
                  raw(), size=2, endian='little')
    y <- writeBin(as.integer(round(256*data$strength)),
                  raw(), size=2, endian='little')
    z <- as.raw(factor(data$forced.call, levels=c('a','h','b','n')))
    dim(x) <- dim(y) <- c(2,nrow(data))
    as.vector(rbind(x,y,z))
}

#---------------------------------------------------------------------

unpack.raw.data <- function(raw, raw.layout)
{
    if (raw.layout == 'signal') {
        .unpack.signal(raw)
    } else if (raw.layout == 'seqread') {
        .unpack.seqread(raw)
    } else if (raw.layout == 'chpdata') {
        .unpack.chpdata(raw)
    } else if (raw.layout == 'dosage') {
        .unpack.dosage(raw)
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
    } else if (raw.layout == 'dosage') {
        .pack.dosage(data)
    } else {
        stop('unknown raw data layout')
    }
}
