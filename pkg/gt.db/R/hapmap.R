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

# Functions for reading unphased HapMap genotype data files
#
# We only support the 'fwd' files, which do not specify the orientations
# of the rs clusters.  The 'rs' files contain orientation errors.

.fixup.hapmap.alleles <- function(alleles, geno)
{
    # fill in non-SNP allele lists based on the actual data
    bad <- grep('^[ACGT]/[ACGT]$', alleles, invert=TRUE)
    if (!length(bad)) return(alleles)
    gt <- c('A','C','G','T')
    tbl <- (ch.table(geno[bad], chars=gt) > 0)
    a12 <- apply(tbl, 1, function(x) c(gt[x],'?','?')[1:2])
    alleles <- as.character(alleles)
    alleles[bad] <- paste(a12[1,], a12[2,], sep='/')
    factor(alleles)
}

.read.hapmap.file <- function(file)
{
    # schlep the whole file into memory: 250 MB should be enough
    fd <- file(file, "r")
    filedata <- readChar(fd, nchars=250e6, useBytes=TRUE)
    close(fd)

    # parse just map information
    tc <- textConnection(gsub("(?m) ncbi_b\\d+ .*$", "",
                              filedata, perl=TRUE))
    map <- scan(tc, what=list('','','',1L,''),
                flush=TRUE, skip=1, quiet=TRUE)
    close(tc)
    names(map) <- c('assay.name','alleles','scaffold','position','strand')
    map$assay.name <- I(gsub('-','_',map$assay.name))
    map$scaffold <- sub('MT','M',map$scaffold)
    map <- as.data.frame(map)

    # now parse genotype data
    filedata <- gsub("(?m)^.* QC\\S+ ", "", filedata, perl=TRUE)
    geno <- strsplit(filedata, "\n", fixed=TRUE)[[1]]
    sample <- strsplit(geno[1], ' ')[[1]]
    geno <- geno[-1]
    map$alleles <- .fixup.hapmap.alleles(map$alleles, geno)
    map$genotype <- .Call('do_encode_gt', as.character(map$alleles),
                          geno, PACKAGE='gt.db')

    stopifnot(!is.na(map$genotype))
    stopifnot(nchar(map$genotype) == length(sample))
    structure(map, sample=sample)
}

.read.hapmap.r22phased <- function(file)
{
    basename <- sub('\\.phased?\\.gz', '', file)
    sample <- read.table(paste(basename,'_sample.txt.gz',sep=''),
                         header=FALSE, as.is=1)
    sample <- paste(rep(sample[,1],each=2), c('A','B'), sep='_')

    map <- read.table(paste(basename,'_legend.txt.gz',sep=''),
                      header=TRUE, as.is=1)
    names(map)[1] <- 'assay.name'
    chr <- sub('genotypes_(chr\\w\\d?).*', '\\1', file, perl=TRUE)
    map$strand <- factor('+')
    map$scaffold <- factor(chr)
    map$alleles <- factor(paste(map$X0, map$X1, sep='/'))
    map$X0 <- map$X1 <- NULL

    fd <- file(file,'r')
    gt <- readChar(fd, nchars=100e6, useBytes=TRUE)
    close(fd)
    gt <- chartr('01','ab',gsub(' ', '', gt, perl=TRUE))
    gt <- strsplit(gt, '\n', fixed=TRUE)[[1]]
    gt <- sapply(gt, charToRaw, USE.NAMES=FALSE)
    map$genotype <- apply(t(gt), 2, rawToChar)

    structure(map, sample = sample[1:ncol(gt)])
}

.resolve.hapmap.positions <- function(data)
{
    # if the same rsid appears at two positions, assign new names
    map <- unique(do.call("rbind", lapply(data, "[", -6)))
    bad <- map$assay.name[which(duplicated(map$assay.name))]
    for (rsid in bad) {
        pos <- unique(map[map$assay.name == rsid,]$position)
        if (length(pos) > 1) {
            warning(rsid, ": found at multiple map positions",
                    call.=FALSE, immediate.=TRUE)
            for (i in 1:length(data)) {
                rn <- which(data[[i]]$assay.name == rsid)
                data[[i]]$assay.name[rn] <-
                    paste(rsid, match(data[[i]]$position[rn], pos), sep='_')
            }
        }
    }
    data
}

.resolve.hapmap.alleles <- function(data)
{
    map <- unique(do.call("rbind", lapply(data, "[", -6)))
    bad <- map$assay.name[which(duplicated(map$assay.name))]
    for (rsid in bad) {
        m <- map[map$assay.name == rsid,]
        tbl <- ch.table(paste(m$alleles,collapse='/'),
                        chars=c('A','C','G','T'))[1,]
        a12 <- c(names(tbl[tbl>0]),'?')
        if (length(a12) > 3)
            stop(rsid,': more than two alleles')
        new <- paste(a12[1], a12[2], sep='/')
        for (i in 1:length(data)) {
            rn <- match(rsid, data[[i]]$assay.name)
            if (is.na(rn)) next
            # convert B/? to A/B genotypes?
            if (data[[i]]$alleles[rn] == sprintf('%s/?',a12[2]))
                data[[i]]$genotype[rn] <-
                    chartr('ab','ba',data[[i]]$genotype[rn])
            data[[i]]$alleles[rn] <- new
        }
    }
    data
}

.merge.hapmap.files <- function(files, verbose=TRUE)
{
    read.file <- if (all(grepl('_r22_.*phase', files)))
        .read.hapmap.r22phased
    else
        .read.hapmap.file
    wrap.read.file <- function(file)
    {
        part <- sub('.*(chr..?_[^_]+)_.*', '\\1', file)
        if (verbose) message('Reading ', part, '...')
        read.file(file)
    }
    all <- lapply(files, wrap.read.file)

    if (verbose) message('Merging...')

    # resolve inconsistencies in positions, alleles
    all <- .resolve.hapmap.positions(all)
    all <- .resolve.hapmap.alleles(all)

    # construct unified map
    map <- unique(do.call('rbind', lapply(all, '[', -6)))
    map <- map[order(map$position),]

    # verify that all files are consistent
    stopifnot(!anyDuplicated(map$assay.name))

    # concatenate genotype strings
    map$genotype <- ''
    sample <- character(0)
    for (d in all) {
        xx <- paste(rep('x', nchar(d$genotype[1])), collapse='')
        gt <- d$genotype[match(map$assay.name,d$assay.name)]
        map$genotype <- paste(map$genotype, if.na(gt,xx), sep='')
        sample <- c(sample, attr(d,'sample'))
    }
    structure(map, sample=sample)
}

load.hapmap.data <-
    function(files, project.name='HapMap', map=TRUE, verbose=TRUE)
{
    # split filenames at '_', and at '.' unless after a digit
    f <- strsplit(basename(files), '_|(?<!\\d)\\.', perl=TRUE)
    f <- as.data.frame(do.call('rbind', f)[,1:7])
    names(f) <- c('gt','chr','pop','rel','ver','build','strand')
    if (grepl('\\.phase', files[1]))
        f$ver <- paste(f$ver, '_phased', sep='')
    stopifnot(f$gt == 'genotypes', f$strand == 'fwd',
              f$rel == f$rel[1], f$ver == f$ver[1],
              f$build == f$build[1])
    nchr <- length(unique(f$chr))
    npop <- length(unique(f$pop))
    stopifnot(nrow(f) == nchr * npop)

    phase <- (if (grepl('phase3',f$rel[1])) 3 else 2)
    platform <- paste('HapMap', phase, sep='')
    mapping <- sub('phase3\\.','r',f$rel[1])
    dataset <- paste(platform, mapping, f$ver[1], sep='_')
    if (npop == 1) dataset <- paste(dataset, f$pop[1], sep='_')
    desc <- sprintf('HapMap Phase %d', phase)
    if (!nrow(ls.platform(platform))) {
        message("Creating platform '", platform, "'...")
        mk.platform(platform, desc)
    }
    desc <- sprintf('%s Release %s', desc, f$rel[1])
    if (!nrow(ls.mapping(platform, mapping))) {
        message("Creating mapping '", mapping, "'...")
        mk.mapping(platform, mapping, desc, sprintf('ncbi_%s', f$build[1]))
    }
    desc <- sprintf('%s (%s)', desc, f$ver[1])
    if (!nrow(ls.dataset(project.name, dataset))) {
        message("Creating dataset '", dataset, "'...")
        mk.dataset(dataset, project.name, platform, desc)
    }

    for (chr in levels(.sort.levels(f$chr))) {
        m <- .merge.hapmap.files(sort(files[f$chr == chr]), verbose)
        s <- ls.sample(dataset)
        if (!nrow(s)) {
            data(hapmap)
            sn <- attr(m,'sample')
            su <- sub('(NA\\d{5}).*', '\\1', sn)
            s <- data.frame(sample.name=sn, subject.name=su,
                            gender=hapmap.subjects[su,'gender.ref'],
                            position=1:length(sn))
            mk.sample(dataset, s)
        } else {
            stopifnot(s$sample.name == attr(m,'sample'))
        }
        rs <- ifelse(!grepl('rs\\d+',m$assay.name), NA,
                     sub('.*rs(\\d+).*','\\1',m$assay.name))
        m$dbsnp.rsid <- as.integer(rs)
        # FIXME: we'll need to fix up pseudo-autosomal regions
        m$ploidy <- switch(chr, 'chrX'='X', 'chrY'='Y', 'chrMT'='M', 'A')
        if (verbose) message("Loading ", nrow(m), " records...")
        if (map) {
            mk.assay(platform, m, progress=verbose)
            mk.assay.position(platform, mapping, m, progress=verbose)
        }
        mk.assay.data(dataset, m, progress=verbose)
    }
}
