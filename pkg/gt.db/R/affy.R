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

# Code to handle importing Affymetrix SNP annotations and genotype
# data into GT.DB

.read.affy.header <- function(file)
{
    txt <- character()
    repeat {
        x <- grep('^\\#%.*=', readLines(file, 100), value=TRUE)
        if (!length(x)) break;
        txt <- c(txt,x)
    }
    val <- sub('^..[^=]*=(.*)', '\\1', txt)
    structure(val, names=sub('^..([^=]*)=.*', '\\1', txt))
}

read.affy.anno <- function(file)
{
    cc <- c(Minor.Allele='NULL',Minor.Allele.Frequency='NULL',
            Associated.Gene='NULL',Genetic.Map='NULL',
            Microsatellite='NULL',Allele.Frequencies='NULL',
            Heterozygous.Allele.Frequencies='NULL',OMIM='NULL',
            Number.of.individuals.Number.of.chromosomes='NULL',
            dbSNP.RS.ID='character',Flank='character')
    d <- read.table(file=file, header=TRUE, quote='"', comment.char='#',
                    sep=',', na.strings='---', row.names=1, colClasses=cc)
    d$Chromosome <- .sort.levels(d$Chromosome)
    d <- d[order(d$Chromosome, d$Physical.Position),]
    tags <- .read.affy.header(file)
    do.call('structure', c(list(d), as.list(tags)))
}

.parse.affy.assays <- function(anno)
{
    alleles <- paste(anno$Allele.A, anno$Allele.B, sep='/')
    data.frame(assay.name=sub('-','_',rownames(anno)),
               alleles=alleles,
               probe.seq=sub('\\[.*\\]','_',anno$Flank),
               stringsAsFactors=FALSE)
}

.parse.affy.mapping <- function(anno)
{
    ploidy <- factor('A',levels=c('A','M','X','Y'))
    ploidy[is.na(anno$Chromosome)] <- NA
    ploidy[anno$Chromosome=='X'] <- 'X'
    ploidy[anno$Chromosome=='Y'] <- 'Y'
    ploidy[anno$Chromosome=='MT'] <- 'M'
    ploidy[anno$ChrX.pseudo.autosomal.region.1 |
           anno$ChrX.pseudo.autosomal.region.2] <- 'A'
    dbsnp.orient <-
        factor(anno$Strand.Versus.dbSNP,
               levels=c('reverse','same'), labels=c('-','+'))
    scaffold <- na.if(paste('chr', anno$Chromosome, sep=''), 'ChrNA')
    scaffold <- sub('chrMT','chrM',scaffold)
    data.frame(assay.name=sub('-','_',rownames(anno)),
               scaffold=scaffold,
               position=anno$Physical.Position,
               strand=anno$Strand,
               ploidy=ploidy,
               dbsnp.rsid=as.integer(sub('^rs','',anno$dbSNP.RS.ID)),
               dbsnp.orient=dbsnp.orient,
               stringsAsFactors=FALSE)
}

load.affy.platform <- function(anno, description, progress=TRUE)
{
    platform <- attr(anno,'chip_type')
    cat(sprintf("Creating platform '%s'...\n", platform))
    mk.platform(platform, description)
    mk.assay(platform, .parse.affy.assays(anno), progress=progress)
}

load.affy.mapping <- function(anno, progress=TRUE)
{
    platform <- attr(anno,'chip_type')
    version <- attr(anno, 'netaffx-annotation-netaffx-build')
    mapping <- sprintf('na%s', version)
    cat(sprintf("Creating mapping '%s' for '%s'...\n",
                mapping, platform))
    ncbi.build <- attr(anno, 'genome-version-ncbi')
    dbsnp.build <- attr(anno, 'dbsnp-version')
    desc <- sprintf("%s on NCBI b%s, dbSNP b%s",
                    ls.platform(platform)$description,
                    ncbi.build, dbsnp.build)
    mk.mapping(platform, mapping, desc,
               sub('\\(d+).*','\\1',ncbi.build))
    mk.assay.position(platform, mapping,
                      .parse.affy.mapping(anno),
                      progress=progress)
}

#---------------------------------------------------------------------

# read a single Affy CHP file (in text form) and convert individual
# genotypes to a GT.DB representation

.read.affy.chp <- function(file, names, row.names=TRUE)
{
    recode.gt <- function(x)
    {
        factor(x, levels=c('AA','AB','BB','NC'),
               labels=c('a','h','b','n'))
    }

    cc <- c('character','factor',rep('numeric',3),'factor')
    d <- read.table(file, sep='\t', header=TRUE,
                    row.names=1, colClasses=cc)
    if (!missing(names)) {
        if (!all(names %in% rownames(d)))
            stop('Probe set name(s) are not matched?')
        d <- d[names,]
    }
    if (!row.names) row.names(d) <- NULL

    colnames(d) <- tolower(colnames(d))
    d$forced.call <- recode.gt(d$forced.call)
    qscore <- round(-10*log10(d$confidence))
    qscore <- ifelse(is.finite(qscore), qscore, 0)
    raw.data <- rawToHex(matrix(.pack.chpdata(d), nrow=5))
    data.frame(genotype=recode.gt(d$call),
               qscore=as.integer(qscore),
               raw.data=raw.data,
               row.names=rownames(d),
               stringsAsFactors=FALSE)
}

# read a list of CHP files, and construct a set of temporary files
# made up of subsets of SNPs across all samples, to facilitate
# transposing the data to the form that will go into GT.DB

.reorg.chp.data <-
    function(files, names, slice=8192, progress=TRUE)
{
    nr <- nrow(.read.affy.chp(files[1], names))
    lo <- seq(1, nr, slice)
    hi <- c(lo[-1]-1, nr)
    nb <- length(lo)
    parts <- sprintf("part%03d.dat", 1:nb)
    if (length(parts) > 120)
        stop("too many slices")
    fd <- lapply(parts, gzfile, 'wb')
    sapply(fd, function(x) save(files, file=x))
    nf <- length(files)
    if (progress) progress.bar(0, nf)
    for (n in 1:nf) {
        d <- .read.affy.chp(files[n], names, row.names=(n==1))
        sapply(1:nb, function(x)
           { slice <- d[lo[x]:hi[x],]
             save(slice, file=fd[[x]]) })
        if (progress) progress.bar(n, nf)
    }
    sapply(fd, close)
    parts
}

.reshape.affy.block <- function(file)
{
    pack.by.assay <- function(data, col, ...)
    {
        # peel out selected column into a matrix, and transpose
        m <- t(sapply(data, function(x) x[,col]))
        # pack columns into strings
        if (is.character(m))
            apply(m, 2, paste, collapse='')
        else
            apply(m, 2, function(x) rawToHex(writeBin(x,raw(0),...)))
    }
    load.one <- function(fd)
    {
        tmpenv <- new.env(parent=emptyenv())
        name <- load(fd, tmpenv)
        get(name, envir=tmpenv)
    }

    fd <- gzfile(file, 'rb')
    files <- load.one(fd)
    d <- replicate(length(files), load.one(fd), simplify=FALSE)
    close(fd)
    data.frame(assay.name=sub('-','_',row.names(d[[1]])),
               genotype=pack.by.assay(d,'genotype'),
               qscore=pack.by.assay(d,'qscore',size=1),
               raw.data=pack.by.assay(d,'raw.data'),
               stringsAsFactors=FALSE)
}

load.affy.chp.data <- function(dataset, anno, files, progress=TRUE)
{
    cat("Reorganizing CHP data...\n")
    parts <- .reorg.chp.data(files, rownames(anno), progress=progress)
    cat("Loading data into database...\n")
    if (progress) progress.bar(0, length(parts))
    for (n in 1:length(parts)) {
        d <- .reshape.affy.block(parts[n])
        mk.assay.data(dataset, d)
        if (progress) progress.bar(n, length(parts))
    }
    unlink(parts)
}
