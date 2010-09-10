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

.scale.gt.matrix <- function(m, method='Price')
{
    m <- data.matrix(m)
    if (method == 'Price') {
        # As in Price et al (2006): center data, scale to uniform
        # variance, and replace missing genotypes with zeros
        p <- apply(m, 2, function(x) 0.5*mean(c(1,x),na.rm=TRUE))
        s <- sqrt(p*(1-p))
    } else if (method == 'sd') {
        s <- apply(m, 2, sd, na.rm=TRUE)
        s <- ifelse(if.na(s, 0), 1, s)
    } else {
        stop("Unknown scaling method")
    }
    t(if.na(scale(m, center=TRUE, scale=s), 0))
}

prcomp.gt.data <-
function(x, sample.mask=TRUE, nc=20, ...)
{
    m <- unpack.gt.matrix(x)
    rn <- row.names(m)
    m <- .scale.gt.matrix(m[sample.mask,])
    pc <- prcomp(m, center=FALSE, retx=FALSE)
    pv <- as.data.frame(pc$rotation[,1:nc])[rn,]
    sd <- structure(pc$sdev[1:nc], names=colnames(pv))
    list(loadings=structure(pv, row.names=rn), sdev=sd,
         assays=nrow(m), dataset.name=attr(x,'dataset.name'),
         call=match.call())
}

.sum.cov.wt <- function(c1, c2)
{
    if (is.null(c1)) return(c2)
    n1 <- c1$n.obs ; n2 <- c2$n.obs
    n.obs <- n1 + n2
    center <- (c1$center*n1 + c2$center*n2) / n.obs
    s <- n1 * crossprod(t(center-c1$center)) +
         n2 * crossprod(t(center-c2$center))
    if (length(s)==1) s <- as.vector(s)
    cov <- (c1$cov*(n1-1) + c2$cov*(n2-1) + s) / (n.obs-1)
    list(cov=cov, center=center, n.obs=n.obs)
}

prcomp.gt.dataset <-
function(x, sample.mask=TRUE, nc=20, ...)
{
    cov.gt.data <- function(gt.data)
    {
        m <- .scale.gt.matrix(unpack.gt.matrix(gt.data)[sample.mask,])
        cov.wt(m, center=FALSE)
    }
    cs <- apply.gt.dataset(x, cov.gt.data, .sum.cov.wt)
    pc <- princomp(covmat=cs)
    rn <- fetch.pt.data(x$dataset.name)$sample.name
    pv <- as.data.frame(pc$loadings[,1:nc])[rn,]
    colnames(pv) <- paste('PC', 1:nc, sep='')
    sd <- structure(pc$sdev[1:nc], names=colnames(pv))
    list(loadings=structure(pv, row.names=rn),
         sdev=sd, assays=cs$n.obs, dataset.name=x$dataset.name,
         gt.dataset=x, call=match.call())
}

#---------------------------------------------------------------------

.rescale.loadings <- function(x, n)
{
    sdev <- sqrt(colSums(x^2))
    s <- scale(x, center=FALSE, scale=sdev)
    list(loadings=as.data.frame(s), sdev=sdev/sqrt(n-1))
}

.snp.loadings.gt.data <- function(x, gt.data)
{
    m <- unpack.gt.matrix(gt.data, names=rownames(gt.data))
    dr <- !is.na(x$loadings[,1])
    p <- data.matrix(x$loadings)[dr,]
    s <- gt.data[, c('assay.name','scaffold','position',
                     'strand','alleles')]
    l <- .scale.gt.matrix(m[dr,]) %*% p
    r <- if (!is.null(attr(gt.data,'part')))
        list(loadings=as.data.frame(l))
    else
        .rescale.loadings(l, sum(dr))
    c(r, list(assays=s, samples=sum(dr),
              dataset.name=attr(gt.data,'dataset.name')))
}

.snp.loadings.gt.dataset <- function(x, gt.dataset)
{
    aggr.fn <- function(a, b)
    {
        if (is.null(a)) return(b)
        list(loadings=rbind(a$loadings,b$loadings),
             assays=rbind(a$assays,b$assays),
             samples=a$samples)
    }
    r <- apply.gt.dataset(gt.dataset, .snp.loadings.gt.data,
                          aggr.fn, x=x)
    c(.rescale.loadings(r$loadings, r$samples),
      list(assays=r$assays, samples=r$samples,
           dataset.name=gt.dataset$dataset.name))
}

snp.loadings <- function(x, data)
{
    if (missing(data) & !is.null(x$gt.dataset))
        data <- x$gt.dataset
    # ugly: dispatch based on second argument
    w <- inherits(data, c('gt.data','gt.dataset'), TRUE)
    if (!any(w)) stop('unknown data class')
    if (w[1])
        .snp.loadings.gt.data(x, data)
    else
        .snp.loadings.gt.dataset(x, data)
}

.apply.loadings.gt.data <- function(x, gt.data)
{
    # FIXME: check platform names
    i <- intersect(rownames(x$loadings), rownames(gt.data))
    if (length(i) > 0) {
        xl <- x$loadings[i,]
        g <- gt.data[i,]
    } else {
        m <- match.gt.data(x$assays, gt.data, by='position')
        xl <- x$loadings[m[,1],]
        g <- orient.gt.data(gt.data[m[,2],], m$flip, m$swap)
    }
    m <- unpack.gt.matrix(g)
    p <- as.data.frame(t(.scale.gt.matrix(m)) %*% data.matrix(xl))
    if (is.null(attr(gt.data,'part')))
        r <- .rescale.loadings(p, ncol(m))
    else
        r <- list(loadings=as.data.frame(p))
    c(r, list(assays=ncol(m),
              dataset.name=attr(gt.data,'dataset.name')))
}

.apply.loadings.gt.dataset <- function(x, gt.dataset)
{
    add.fn <- function(a, b)
    {
        if (is.null(a)) return(b)
        list(loadings=a$loadings+b$loadings,
             assays=a$assays+b$assays)
    }
    p <- apply.gt.dataset(gt.dataset, .apply.loadings.gt.data,
                          add.fn, x=x)
    c(.rescale.loadings(p$loadings, p$assays),
      list(assays=p$assays, dataset.name=gt.dataset$dataset.name))
}

apply.loadings <- function(x, data)
{
    w <- inherits(data, c('gt.data','gt.dataset'), TRUE)
    if (!any(w)) stop('unknown data class')
    if (w[1])
        .apply.loadings.gt.data(x, data)
    else
        .apply.loadings.gt.dataset(x, data)
}

#---------------------------------------------------------------------

big.loadings <- function(x, sigma=6)
{
    l <- x$loadings
    if (is.data.frame(x$assays)) {
        id <- x$assays
    } else {
        id <- data.frame(sample.name=rownames(l))
    }
    cc <- complete.cases(l)
    l <- l[cc,] ; id <- id[cc,]
    r <- do.call('rbind', lapply(1:ncol(l), function(n) {
        r <- data.frame(id, name=colnames(l)[n], value=l[,n])
        r$zscore <- (r$value-mean(r$value))/sd(r$value)
        r$var <- r$zscore^2/sum(r$zscore^2)
        subset(r, abs(zscore)>sigma)
    }))
    r$name <- factor(r$name, levels=names(l))
    structure(r, row.names=seq_along(r[,1]))
}

#---------------------------------------------------------------------

ls.prcomp <-
function(dataset.name, prcomp.name='%',
         show.all=FALSE, show.ids=FALSE)
{
    sql <-
     'select prcomp_id, name prcomp_name, description,
             fn_call, components, samples, assays,
             is_hidden, created_by, created_dt
      from prcomp
      where dataset_id=:1
        and name like :2
        and is_hidden<=:3'
    dset.id <- lookup.id('dataset', dataset.name)
    r <- sql.query(gt.db::.gt.db, sql, dset.id, prcomp.name, show.all)
    .filter.ids(r, show.ids)
}

lookup.prcomp.id <- function(dataset.name, prcomp.name)
{
    if (missing(prcomp.name)) {
        p <- ls.prcomp(dataset.name, show.ids=TRUE)
        if (nrow(p) != 1)
            stop('could not choose default prcomp', call.=FALSE)
        structure(p$prcomp.id, names=p$prcomp.name)
    } else {
        dset.id <- lookup.id('dataset', dataset.name)
        lookup.id('prcomp', prcomp.name, dataset.id=dset.id)
    }
}

rm.prcomp <- function(dataset.name, prcomp.name)
{
    id <- lookup.id('prcomp', prcomp.name,
                    dataset.id=lookup.id('dataset', dataset.name))
    sql <- 'delete from prcomp where prcomp_id=:1'
    sql.exec(gt.db::.gt.db, sql, id)
}

store.prcomp <-
function(x, prcomp.name, description,
         nc=ncol(x$loadings), is.hidden=FALSE)
{
    .check.name(prcomp.name)
    sql <-
     "insert into prcomp
      values (null,:1,:2,:3,:4,:5,:6,:7,:8,:user:,:sysdate:)"
    dset.id <- lookup.id('dataset', x$dataset.name)
    v <- subset(x$loadings, !is.na(x$loadings[,1]))
    sql.exec(gt.db::.gt.db, sql, dset.id, prcomp.name, description,
             .as.text(x$call), nc, nrow(v), x$assays, is.hidden)
    id <- lookup.id('prcomp', prcomp.name, dataset.id=dset.id)
    sql.exec(gt.db::.gt.db, 'insert into prcomp_component values (:1,:2,:3)',
             id, 1:nc, x$sdev)
    samp.id <- lookup.id('sample', rownames(v), dataset.id=dset.id)
    for (n in 1:nc) {
        sql.exec(gt.db::.gt.db,
                 'insert into prcomp_loading values (:1,:2,:3,:4)',
                 id, n, samp.id, v[,n])
    }
    nc
}

fetch.prcomp <- function(dataset.name, prcomp.name, nc)
{
    pca.id <- lookup.prcomp.id(dataset.name, prcomp.name)
    info <- ls.prcomp(dataset.name, names(pca.id))
    if (missing(nc))
        nc <- info$components
    sql <- 'select * from prcomp_component
            where prcomp_id=:1 and component<=:2
            order by component'
    comp <- sql.query(gt.db::.gt.db, sql, pca.id, nc)

    s <- subset(ls.sample(dataset.name), !is.na(position))
    s <- s[order(s$position),]

    sql <- 'select position, loading
            from sample s, prcomp_loading p
            where s.sample_id=p.sample_id
              and prcomp_id=:1
              and component=:2'
    loadings <- NULL
    for (n in 1:nc) {
        pc <- sql.query(gt.db::.gt.db, sql, pca.id, n)
        v <- pc$loading[match(s$position,pc$position)]
        loadings <- cbind(loadings, v)
    }
    colnames(loadings) <- cn <- paste('PC', comp$component, sep='')
    rownames(loadings) <- s$sample.name
    list(loadings=as.data.frame(loadings),
         sdev=structure(comp$sdev,names=cn),
         assays=info$assays, dataset.name=dataset.name,
         call=parse(text=info$fn.call)[[1]])
}

#---------------------------------------------------------------------

gplot.prcomp <-
function(x, col=1:5, aggr.fn=function(x) sum(x^2), rescale=TRUE,
         xlab=NULL, ylab=NULL, ...)
{
    d <- x$loadings[col]
    cl <- names(d)
    nc <- ncol(d)
    for (n in 1:nc) {
        gr=list(inside=list(x=0.93,y=0.1,corner=c(1,1)/2,fun=textGrob(cl[n])))
        print(gplot(~d[cl[n]], x$assays, aggr.fn=aggr.fn,
                    rescale=rescale, legend=gr, colorkey=FALSE,
                    xlab=xlab, ylab=ylab, draw.at=FALSE, ...),
              newpage=(n==1), split=c(1,n,1,nc))
    }
}

qqprcomp <- function(x, col=1:6, layout, ...)
{
    x <- x$loadings
    cl <- names(x)[col]
    np <- length(cl)
    if (missing(layout)) {
        ddim <- par("din")
        ny <- max(1, round(sqrt(np*ddim[2]/ddim[1])))
        nx <- ceiling(np/ny)
        ny <- ceiling(np/nx)
    } else {
        nx <- layout[1]
        ny <- layout[2]
    }
    for (n in 1:np) {
        ix <- (n-1) %% nx + 1 ; iy <- (n-1) %/% nx + 1
        gr=list(inside=list(x=0.95, y=0.05, corner=c(1,0),
                            fun=textGrob(cl[n])))
        print(qqthin(~x[,cl[n]], xlab=NULL, ylab=NULL, aspect=1,
                     scales=list(draw=0), legend=gr, ...),
              newpage=(n==1), split=c(ix,iy,nx,ny))
    }
}

score.prcomp <- function(formula, pc.data, pt.data, ...)
{
    fn <- function(PC)
    {
        pt.data$PC <- PC
        s <- summary(lm(formula, pt.data, ...))
        p <- do.call('pf',c(unname(as.list(s$fstatistic)),lower.tail=FALSE))
        c(s$adj.r.squared, s$fstatistic[1], p)
    }
    d <- as.data.frame(t(sapply(pc.data$loadings, fn)))
    colnames(d) <- c('R^2', 'F value', 'Pr(>F)')
    class(d) <- c('anova','data.frame')
    d
}
