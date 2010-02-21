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

.read.hapmap.file <- function(file)
{
    # read just map information
    map <- scan(file, what=list('','','',1L,''),
                flush=TRUE, skip=1, quiet=TRUE)
    names(map) <- c('assay.name','alleles','scaffold','position','strand')
    map$assay.name <- I(gsub('-','_',map$assay.name))
    map$scaffold <- sub('MT','M',map$scaffold)
    map <- as.data.frame(map)

    # now parse genotype data
    geno <- readLines(file)
    # note that we keep the leading space
    geno <- sub('.* QC[^ ]+', '', geno, perl=TRUE)
    sample <- strsplit(geno[1], ' ')[[1]][-1]
    geno <- geno[-1]
    for (a in levels(map$alleles)) {
        w <- (map$alleles == a)
        ab <- strsplit(a,'/')[[1]]
        geno[w] <- chartr(sprintf('%s%sN',ab[1],ab[2]),'abn',geno[w])
    }
    geno <- gsub(' ab| ba', 'h', geno, perl=TRUE)
    map$genotype <- gsub(' .', '', geno, perl=TRUE)
    stopifnot(all(nchar(map$genotype) == length(sample)))
    structure(map, sample=sample)
}

.merge.hapmap.files <- function(files, verbose=TRUE)
{
    wrap.read.file <- function(file)
    {
        part <- sub('.*(chr..?_[^_]+)_.*', '\\1', file)
        if (verbose) message('Reading ', part, '...')
        .read.hapmap.file(file)
    }
    all <- lapply(files, wrap.read.file)

    # construct unified map
    map <- do.call('rbind', lapply(all, '[', -6))
    map <- unique(map)
    map <- map[order(map$position),]

    # verify that all files are consistent
    stopifnot(!anyDuplicated(map$assay.name))

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
    function(files, project.name='HapMap', verbose=TRUE)
{
    # split filenames at '_', and at '.' unless after a digit
    f <- strsplit(basename(files), '_|(?<!\\d)\\.', perl=TRUE)
    f <- as.data.frame(do.call('rbind', f)[,1:7])
    names(f) <- c('gt','chr','pop','rel','ver','build','strand')
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
    if (!nrow(ls.platform(platform)))
        mk.platform(platform, desc)
    desc <- sprintf('%s Release %s', desc, f$rel[1])
    if (!nrow(ls.mapping(platform, mapping)))
        mk.mapping(platform, mapping, desc, sprintf('ncbi_%s', f$build[1]))
    desc <- sprintf('%s (%s)', desc, f$ver[1])
    if (!nrow(ls.dataset(project.name, dataset)))
        mk.dataset(dataset, project.name, platform, desc)

    for (chr in levels(.sort.levels(f$chr))) {
        m <- .merge.hapmap.files(sort(files[f$chr == chr]), verbose)
        s <- ls.sample(dataset)
        if (!nrow(s)) {
            data(hapmap)
            sn <- attr(m,'sample')
            s <- data.frame(sample.name=sn, subject.name=sn,
                            gender=hapmap.subjects[sn,'gender.ref'],
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
        mk.assay(platform, m)
        mk.assay.position(platform, mapping, m)
        mk.assay.data(dataset, m)
    }
}
