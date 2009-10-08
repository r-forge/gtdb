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

#-----------------------------------------------------------------------

.as.text <- function(expr)
{
    gsub('\\s+', ' ', paste(deparse(expr), collapse=''))
}

.null.model <- function(formula, term='genotype')
{
    txt <- gsub(paste('[*+:]',term), '', .as.text(formula))
    eval(parse(text=gsub(term, '1', txt)))
}

.flat.model <- function(formula, term='genotype')
{
    txt <- .as.text(.null.model(formula, term))
    eval(parse(text=paste(txt, '+', term)))
}

.lhs.formula <- function(formula) as.character(formula[2])
.rhs.formula <- function(formula) as.character(formula[3])

.recode.gt <- function(genotype, mode)
{
    if (mode=='additive')
        return(genotype)
    else if (mode=='recessive')
        return(as.numeric(genotype==2))
    else if (mode=='dominant')
        return(sign(genotype))
    else if (mode=='general')
        return(factor(genotype,levels=0:2))
    else
        stop('unknown mode-of-action')
}

#-----------------------------------------------------------------------

score.gt.data <-
function(formula, pt.data, gt.data, score.fn=NULL,
         pt.filter=TRUE, gt.filter=TRUE, progress=FALSE, ...)
{
    pt.filter <- eval(substitute(pt.filter), pt.data, parent.frame())
    gt.filter <- eval(substitute(gt.filter), gt.data, parent.frame())
    gt.data <- gt.data[if.na(gt.filter,FALSE),]

    # parse formula to get list of variables used
    ft <- attr(terms(formula),'factor')
    col <- dimnames(ft)[[1]]
    dsub <- get_all_vars(formula, cbind(pt.data, genotype=NA))
    dsub$genotype <- NULL

    # special case: chisq.2x2 needs gender information
    if (deparse(substitute(score.fn))=='score.chisq.2x2')
        dsub$gender <- pt.data$gender

    # guess scoring function based on model and outcome
    if (is.null(score.fn)) {
        binary <- (length(table(pt.data[,col[1]])) == 2)
        groups <- (length(grep(':genotype',dimnames(ft)[[2]])))
        fn <- if (!groups && binary && (length(col)==2))
            'score.trend'
        else paste('score', c('.lm','.glm')[binary+1],
                   c('','.groups')[groups+1], sep='')
        warning("Using score.fn = ", fn, call.=FALSE, immediate.=TRUE)
        score.fn <- eval(parse(text=fn))
    }

    nr <- nrow(gt.data)
    wrap.score.fn <- function(n) {
        gn <- gt.data[n,]
        my.error <- function(e) {
            e$message <- paste('assay', gn$assay.name, ':', e$message)
            .signalSimpleWarning(e$message,e$call)
        }
        my.warn <- function(e) {
            my.error(e)
            invokeRestart('muffleWarning')
        }
        d <- cbind(dsub, unpack.gt.matrix(gn, 'genotype'))
        keep <- if.na(pt.filter,FALSE) & complete.cases(d)
        r <- tryCatch(withCallingHandlers(
            score.fn(formula, d[keep,], gn$ploidy, ...),
            warning=my.warn), error=my.error)
        if (progress) progress.bar(n, nr)
        if (!is.null(r)) {
            id <- gn[c('assay.data.id','assay.name')]
            keep.attr(merge(id,r), .Attr=kept.attr(r))
        }
    }
    if (progress) progress.bar(0, nr)
    r <- lapply(1:nr, wrap.score.fn)

    names(r) <- rownames(gt.data)
    r <- do.call('rbind', r)
    if (progress) cat("Scored", (if (is.null(r)) 0 else nrow(r)),
                      "out of", nr, "assays.\n")

    if (is.null(r)) return(data.frame())
    keep.attr(r, model=.as.text(formula), call=.as.text(match.call()),
              dataset.name=attr(gt.data,'dataset.name'))
}

ls.test <-
function(dataset.name, test.name='%',
         show.all=FALSE, show.ids=FALSE)
{
    sql <-
     'select test_id, name test_name, description, fit, model,
             term, is_hidden, created_by, created_dt
      from test
      where dataset_id=:1
        and name like :2
        and is_hidden<=:3'
    dset.id <- lookup.id('dataset', dataset.name)
    r <- sql.query(gt.db::.gt.db, sql, dset.id, test.name, show.all)
    .filter.ids(r, show.ids)
}

rm.test <- function(dataset.name, test.name)
{
    id <- lookup.id('test', test.name,
                    dataset.id=lookup.id('dataset', dataset.name))
    sql <- 'delete from test where test_id=:1'
    sql.exec(gt.db::.gt.db, sql, id)
}

store.test.scores <-
function(x, test.name, description, is.hidden=FALSE)
{
    .check.name(test.name)
    ut <- if (is.null(x$term)) NA else unique(x$term)
    dn <- attr(x, 'dataset.name')
    dset.id <- lookup.id('dataset', dn)
    sql <-
     'insert into test
      values (null,:1,:2,:3,:4,:5,:6,:7,:user:,:sysdate:)'
    try(sql.exec(gt.db::.gt.db, sql, dset.id, test.name, description,
                 attr(x,'fit'), attr(x,'model'),
                 ut, is.hidden), silent=TRUE)
    tid <- ls.test(dn, test.name, show.ids=TRUE)
    row.names(tid) <- if.na(tid$term, 'NA')

    for (n in c('term','effect','stderr')) {
        if (is.null(x[[n]])) x[[n]] <- NA
    }
    fmt <- function(y,n) {
        y[!is.finite(y)] <- NA
        if.na(y, NA, format(y, trim=TRUE, scientific=-20, digits=7))
    }
    for (n in c('pvalue','effect','stderr')) {
        x[[n]] <- fmt(x[[n]])
    }

    sql <- 'insert into test_result values(:1,:2,:3,:4,:5)'
    sql.exec(gt.db::.gt.db, sql, tid[if.na(x$term,'NA'),'test.id'],
             x[c('assay.data.id','pvalue','effect','stderr')])
}

fetch.test.scores <-
function(dataset.name, test.name, term, max.pvalue)
{
    dset.id <- lookup.id('dataset', dataset.name)
    tid <- ls.test(dataset.name, test.name, show.ids=TRUE)
    if (!missing(term))
        tid <- tid[if.na(tid$term,'NA') %in% if.na(term,'NA'),]

    sql <-
     'select r.assay_data_id, a.name assay_name,
             term, pvalue, effect, stderr
      from test t, test_result r, assay_data d, assay a
      where t.test_id=:1
        and t.test_id=r.test_id
        and r.assay_data_id=d.assay_data_id
        and d.assay_id=a.assay_id'
    dat <- list(tid$test.id)
    if (!missing(max.pvalue)) {
        sql <- paste(sql, 'and pvalue<=:2')
        dat <- list(tid$test.id, max.pvalue)
    }
    d <- sql.query(gt.db::.gt.db, sql, dat)
    keep.attr(d, test.name=test.name)
}

score.and.store <-
function(dataset.name, test.name, description, formula,
         score.fn=NULL, pt.filter=TRUE, gt.filter=TRUE, pca=FALSE,
         part=1:parts, parts=200, dryrun=FALSE, ...)
{
    p <- fetch.pt.data(dataset.name, pca=pca)
    pf <- eval(substitute(pt.filter), p, parent.frame())

    step <- as.numeric(Sys.getenv('LSB_JOBINDEX_END'))
    if (!is.na(step) && (step > 1)) {
        start <- as.numeric(Sys.getenv('LSB_JOBINDEX'))
        part <- seq(start, parts, step)-1
    }
    nr <- 0
    for (n in part) {
        cat(sprintf("%s: %.0f/%.0f\n", date(), n, parts+0))
        flush.console()
        g <- fetch.gt.data(dataset.name, part=n, parts=parts)
        gf <- eval(substitute(gt.filter), g, parent.frame())
        r <- score.gt.data(formula, p, g, score.fn, pf, gf, ...)
        if (nrow(r) && !dryrun)
	    nr <- nr + store.test.scores(r, test.name, description)
        else
            nr <- nr + nrow(r)
    }
    cat(date(), ': done, ', nr, ' rows inserted\n', sep='')
}
