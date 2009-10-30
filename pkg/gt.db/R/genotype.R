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

.sort.levels <- function(x, decreasing=FALSE)
{
    # sort levels with numeric suffixes into numerical order
    n <- unique(as.character(x))
    a <- sub('([0-9]+)$', '000000000\\1', n)
    a <- sub('0+([0-9]{9})$', '\\1', a)
    factor(x, levels=n[order(a,decreasing=decreasing)])
}

fetch.gt.data <-
function(dataset.name, mapping.name, assay.name, dbsnp.rsid,
         part, parts, by=c('assay','position'), binsz=1, where,
         show.ids=FALSE, genotype=TRUE, qscore=FALSE, raw.data=FALSE)
{
    sql <-
     'select assay_data_id, a.assay_id, a.name assay_name, alleles,
        scaffold, position, strand, ploidy,
        dbsnp_rsid, dbsnp_orient, d.flags d_flags'

    db.mode <- .gt.db.options('db.mode')
    tx.mode <- .gt.db.options('tx.mode')
    if (db.mode == tx.mode) {
        cvt.fn <- ':blob:'
    } else if (db.mode == 'raw' && tx.mode == 'hex') {
        cvt.fn <- ':hex.blob:'
    } else {
        stop('unknown conversion!')
    }

    if (genotype)
        sql <- paste(sql, ', :clob:(genotype)')
    if (qscore)
        sql <- paste(sql, ', ', cvt.fn, '(qscore)', sep='')
    if (raw.data)
        sql <- paste(sql, ', ', cvt.fn, '(raw_data)', sep='')

    sql <- paste(sql,
     'from assay_data d, assay a, assay_position p
      where d.assay_id=a.assay_id
        and a.assay_id=p.assay_id
        and dataset_id=:1
        and mapping_id=:2')

    dset.id <- lookup.id('dataset', dataset.name)
    ds <- ls.dataset(,dataset.name)
    map.id <- lookup.mapping.id(ds$platform.name, mapping.name)
    dat <- .sql.prep.data(dset.id, map.id)
    if (!missing(assay.name)) {
        sql <- paste(sql, 'and a.name=:3')
        dat <- .sql.prep.data(dset.id, map.id, assay.name)
    } else if (!missing(dbsnp.rsid)) {
        sql <- paste(sql, 'and dbsnp.rsid=:3')
        dat <- .sql.prep.data(dset.id, map.id, dbsnp.rsid)
    } else if (!missing(part)) {
        by <- match.arg(by)
        if (by=='position') {
            sql <- paste(sql, 'and :mod:(:floor:(position/:3)%%:4)=:5-1')
        } else {
            sql <- paste(sql, 'and :mod:(:floor:(assay_data_id/:3)%%:4)=:5-1')
        }
        dat <- .sql.prep.data(dset.id, map.id, binsz, parts, part)
    }
    if (!missing(where)) sql <- paste(sql, 'and', where)

    d <- lapply(1:nrow(dat), function(x)
                sql.query(gt.db::.gt.db, sql, dat[x,]))
    d <- do.call('rbind', d)

    for (n in c('assay.data.id','assay.id','position','dbsnp.rsid'))
        d[,n] <- as.integer(d[,n])
    for (n in c('strand','dbsnp.orient'))
        d[,n] <- factor(d[,n], levels=c('+','-'))
    d$ploidy <- .fixup.ploidy(d$ploidy)
    d$scaffold <- .sort.levels(d$scaffold)

    nf <- match('d.flags', names(d))
    d <- cbind(d[1:(nf-1)],
               unpack.flags('assay_data', dataset.name, d$d.flags),
               d[-(1:nf)])

    s <- subset(ls.sample(dataset.name), !is.na(position))
    s <- s[order(s$position),]
    if (!identical(s$position,1:nrow(s)) ||
        any(nchar(d$genotype) != nrow(s))) {
        stop('inconsistent sample layout information!');
    }
    g <- c(F=as.mask(s$gender=='F'), M=as.mask(s$gender=='M'))

    d <- .filter.ids(data.frame(d,row.names=d$assay.name),
                     show.ids, 'assay.data.id')
    class(d) <- c('gt.data','data.frame')
    keep.attr(d, dataset.name=dataset.name,
              raw.layout=ds$raw.layout,
              platform.name=ds$platform.name,
              mapping.name=names(map.id),
              gender=g, sample.name=s$sample.name)
}

#---------------------------------------------------------------------

# convert between packed genotype strings and vectors

gt.split <-
function(s, convert=c('score.b','score.a','char','none'),
         alleles, na.codes=c(), strand=c('+','-'), sep='/')
{
    convert <- match.arg(convert)
    strand <- match.arg(strand)
    join <- function(x) paste(x, collapse=sep)
    lv <- c('a','h','b',na.codes)
    g <- factor(strsplit(s, '')[[1]], levels=lv)
    suppressWarnings(g[is.na(g)] <- 'n')
    if (convert=='char') {
        if (strand=='-') alleles <- revcomp(alleles)
        na.labels <- toupper(paste(na.codes, na.codes, sep=sep))
        la <- c(join(alleles[c(2,2)]), join(sort(alleles)),
                join(alleles[c(1,1)]), na.labels)
        factor(g, levels=lv, labels=la)
    } else if (convert=='score.a') {
        g <- as.integer(g) - 1
        ifelse(g<3, 2-g, g)
    } else if (convert=='score.b') {
        as.integer(g) - 1
    } else {
        g
    }
}

gt.paste <-
function(v, convert=c('score.b','score.a','char','none'),
         alleles, na.codes=c(), strand=c('+','-'), sep='/')
{
    convert <- match.arg(convert)
    strand <- match.arg(strand)
    join <- function(x) paste(x, collapse=sep)
    if (convert=='char') {
        if (strand=='-') alleles <- revcomp(alleles)
        na.labels <- toupper(paste(na.codes, na.codes, sep=sep))
        la <- c(join(alleles[c(2,2)]), join(alleles[c(1,2)]),
                join(alleles[c(2,1)]), join(alleles[c(1,1)]), na.labels)
        v <- factor(v, levels=la, labels=c('a','h','h','b', na.codes))
    } else if (convert=='score.a') {
        v <- factor(v, levels=c(2:0,3:6)[1:(length(na.codes)+2)],
                    labels=c('a','h','b',na.codes))
    } else if (convert=='score.b') {
        v <- factor(v, levels=0:(length(na.codes)+2),
                    labels=c('a','h','b',na.codes))
    }
    paste(if.na(as.character(v),'n'), collapse='')
}

#---------------------------------------------------------------------

unpack.gt.matrix <-
function(gt.data, names=gt.data$assay.name, ...)
{
    if ((nrow(gt.data)>1) && (length(names)==1)) {
        names <- paste(names,1:nrow(gt.data),sep='.')
    }
    a <- strsplit(gt.data$alleles, split='/', fixed=TRUE)
    fn <- function(x) gt.split(gt.data$genotype[x], ..., alleles=a[[x]])
    g <- lapply(1:nrow(gt.data), fn)
    names(g) <- names
    if (is.numeric(g[[1]])) g <- do.call('cbind', g)
    g <- as.data.frame(g, row.names=attr(gt.data,'sample.name'))
    keep.attr(g, ploidy=structure(gt.data$ploidy,names=names))
}

reshape.gt.data <- function(gt.data, ...)
{
    if (nrow(gt.data) > 1) {
        r <- lapply(1:nrow(gt.data), function(x)
                    reshape.gt.data(gt.data[x,], ...))
        return(do.call('rbind', r))
    }

    d <- data.frame(assay.name=gt.data$assay.name,
                    sample.name=attr(gt.data,'sample.name'))
    row.names(d) <- paste(d$assay.name,d$sample.name,sep='.')
    if (!is.null(gt.data$genotype)) {
        a <- strsplit(gt.data$alleles, split='/', fixed=TRUE)[[1]]
        d$genotype <- gt.split(gt.data$genotype, ..., alleles=a)
    }

    tx.mode <- .gt.db.options('tx.mode')
    cvt.fn <- function(x) stop('unknown conversion!')
    if (tx.mode == 'hex')
        cvt.fn <- hexToRaw

    if (!is.null(gt.data$qscore)) {
        d$qscore <- as.integer(cvt.fn(gt.data$qscore))
    }
    if (is.null(gt.data$raw.data))
        return(d)
    f <- attr(gt.data,'raw.layout')
    if (is.null(f))
        stop('raw data layout unavailable')
    nr <- nrow(d)
    if (f == 'signal') {
        i <- readBin(cvt.fn(gt.data$raw.data), what='int',
                     n=2*nr, size=2, signed=FALSE, endian='little')
        i <- na.if(i, 65535)
        dim(i) <- c(2,nr)
        cbind(d, signal.a=i[1,], signal.b=i[2,])
    } else if (f == 'seqread') {
        i <- readBin(cvt.fn(gt.data$raw.data), what='int',
                     n=4*nr, size=1, signed=FALSE, endian='little')
        i <- na.if(i, 255)
        dim(i) <- c(4,nr)
        cbind(d, fwd.a=i[1,], rev.a=i[2,], fwd.b=i[3,], rev.b=i[4,])
    } else {
        warning('unknown raw data layout')
        d
    }
}

.mask.dat <- function(str, mask, squeeze=FALSE)
{
    tx.mode <- .gt.db.options('tx.mode')
    if (tx.mode != 'hex')
        stop('unknown conversion!')
    mask <- un.mask(mask)
    mask <- rep(mask, each=nchar(str)/length(mask))
    mask.str(str, mask, if (squeeze) '' else 'F')
}

mask.gt.data <- function(gt.data, sample.mask, repack=FALSE)
{
    if (!is.null(gt.data$qscore))
        gt.data$qscore   <- .mask.dat(gt.data$qscore, sample.mask, repack)
    if (!is.null(gt.data$raw.data))
        gt.data$raw.data <- .mask.dat(gt.data$raw.data, sample.mask, repack)
    if (repack) {
        gt.data$genotype <-
            mask.str(gt.data$genotype, sample.mask, '')
        g <- mask.str(attr(gt.data,'gender'), sample.mask, '')
        attr(gt.data,'gender') <- c(F=g[1], M=g[2])
        if (!is.logical(sample.mask))
            sample.mask <- un.mask(sample.mask)
        attr(gt.data,'sample.name') <-
            attr(gt.data,'sample.name')[sample.mask]
    } else {
        gt.data$genotype <- mask.str(gt.data$genotype, sample.mask)
    }
    gt.data
}

#---------------------------------------------------------------------

gt.dataset <-
function(dataset.name, gt.filter=TRUE, parts=10,
         by='assay', binsz=1, progress=interactive())
{
    x <- list(dataset.name=dataset.name,
              gt.filter=substitute(gt.filter),
              env=parent.frame(),
              parts=parts, by=by, binsz=binsz,
              progress=progress)
    structure(x, class='gt.dataset')
}

apply.gt.dataset <-
function(gt.dataset, part.fn, aggr.fn, ..., GT.DATA='gt.data')
{
    part.fn <- match.fun(part.fn)
    aggr.fn <- match.fun(aggr.fn)
    part.list <- 1:gt.dataset$parts
    x <- NULL
    if (gt.dataset$progress) progress.bar(0, gt.dataset$parts)
    for (part in part.list) {
        g <- fetch.gt.data(gt.dataset$dataset.name, part=part,
                           parts=gt.dataset$parts, binsz=gt.dataset$binsz)
        attr(g,'part') <- part
        gf <- if.na(eval(gt.dataset$gt.filter, g, gt.dataset$env), FALSE)
        args <- list(gt=g[gf,], ...)
        names(args)[1] <- GT.DATA
        x <- aggr.fn(x, do.call(part.fn, args))
        if (gt.dataset$progress) progress.bar(part, gt.dataset$parts)
    }
    x
}

#---------------------------------------------------------------------

.empty.gt <- function(gt.data)
{
    nr <- nchar(attr(gt.data,'gender')['M'], 'bytes')
    if (length(nr)!=1) nr <- max(c(1,nchar(gt.data$genotype)))
    paste(rep('x',nr), collapse='')
}

summary.gt.data <- function(object, sample.mask, by.sample=FALSE, ...)
{
    count.by.col <- function(str, ch)
    {
        str <- str[!is.na(str) & nchar(str)]
        tot <- integer(max(nchar(str)))
        for (s in str) tot <- tot + (charToRaw(s) == charToRaw(ch))
        tot
    }
    if (!missing(sample.mask) && !identical(sample.mask, TRUE))
        object <- mask.gt.data(object, sample.mask, repack=by.sample)
    g <- object$genotype
    s <- object$ploidy
    if (any(is.na(s))) {
        warning('assays with unknown ploidy will be assumed autosomal')
        s[is.na(s)] <- 'A'
    }
    x <- .empty.gt(object)

    # construct diploid, haploid genotype subsets
    gm <- attr(object,'gender')
    dg <- ifelse(s=='A', g, x)
    hg <- ifelse(s=='M', g, x)
    if (any(s=='X'))
        dg[s=='X'] <- mask.str(dg[s=='X'], gm['F'])
    if (any(s %in% c('X','Y'))) {
        if (any(!(un.mask(gm['M']) | un.mask(gm['F']))))
            warning('some samples have unknown gender')
        hg[s %in% c('X','Y')] <- mask.str(g[s %in% c('X','Y')], gm['M'])
    }

    count <- if (by.sample) count.by.col else nsubstr
    NN <- count(dg, 'n') + count(hg, 'n')
    AA <- count(dg, 'a')
    AB <- count(dg, 'h')
    BB <- count(dg, 'b')
    A_ <- count(hg, 'a')
    B_ <- count(hg, 'b')
    gt.rate <- 1-NN/(NN+AA+AB+BB+A_+B_)

    if (by.sample) {
        d <- data.frame(NN, AA, AB, BB, A_, B_, gt.rate,
                        hz=AB/(AA+AB+BB))
        if (!is.null(attr(object,'sample.name')))
            rownames(d) <- attr(object,'sample.name')
        d
    } else {
        af <- (2*AA+AB+A_)/(2*(AA+AB+BB)+A_+B_)
        data.frame(NN, AA, AB, BB, A_, B_, gt.rate,
                   freq.a=af, freq.b=1-af,
                   hw.p.value=hwe.test(AA,AB,BB),
                   row.names=row.names(object))
    }
}

summary.gt.dataset <-
function(object, sample.mask=TRUE, by.sample=FALSE, ...)
{
    apply.gt.dataset(object, summary, rbind, sample.mask,
                     by.sample, GT.DATA='object')
}

#---------------------------------------------------------------------

.by.allele <- function(x, fn)
{
  sapply(x, function(s) paste(fn(strsplit(s,'/')[[1]]),collapse='/'))
}

match.gt.data <-
function(gt.data.1, gt.data.2, by=c('position','dbsnp.rsid'))
{
    by <- match.arg(by)
    if (by=='position') {
        cols <- c('assay.name','scaffold','position','strand','alleles')
        d1 <- subset(gt.data.1[cols], !is.na(position))
        d2 <- subset(gt.data.2[cols], !is.na(position))
        m <- merge(d1, d2, by=c('scaffold','position'))
        f <- (m$strand.x != m$strand.y)
    } else {
        cols <- c('assay.name','dbsnp.rsid','dbsnp.orient','alleles')
        d1 <- subset(gt.data.1[cols], !is.na(dbsnp.rsid))
        d2 <- subset(gt.data.2[cols], !is.na(dbsnp.rsid))
        m <- merge(d1, d2, by='dbsnp.rsid')
        f <- (m$dbsnp.orient.x != m$dbsnp.orient.y)
    }
    if (any(is.na(f)))
        warning(sum(if.na(f)), ' assays could not be oriented')

    if (any(if.na(f,FALSE))) {
        w <- which(f)
        m$alleles.y[w] <- sapply(m$alleles.y[w], .by.allele, revcomp)
    }
    s <- rep(NA, nrow(m))
    s[m$alleles.x == .by.allele(m$alleles.y, rev)] <- TRUE
    s[m$alleles.x == m$alleles.y] <- FALSE
    if (any(is.na(s)))
        warning(sum(if.na(s)), ' assays had inconsistent alleles')

    d <- data.frame(m$assay.name.x, m$assay.name.y,
                    is.flipped=f, is.swapped=s)
    names(d)[1] <- c(attr(gt.data.1,'dataset.name'),'gt.data.1')[1]
    names(d)[2] <- c(attr(gt.data.2,'dataset.name'),'gt.data.2')[1]
    subset(d, !is.na(f) & !is.na(s))
}

orient.gt.data <- function(gt.data, flip=FALSE, swap=FALSE)
{
    flip <- if.na(flip,FALSE)
    swap <- if.na(swap,FALSE)
    if (any(swap)) {
        gt.data$alleles[swap] <- .by.allele(gt.data$alleles[swap], rev)
        gt.data$genotype[swap] <- chartr('ab', 'ba', gt.data$genotype[swap])
    }
    if (any(flip)) {
        gt.data$strand[flip] <- chartr('-+', '+-', gt.data$strand[flip])
        gt.data$dbsnp.orient[flip] <-
            chartr('-+', '+-', gt.data$dbsnp.orient[flip])
        gt.data$alleles[flip] <-
            sapply(gt.data$alleles[flip], .by.allele, revcomp)
    }
    gt.data
}
