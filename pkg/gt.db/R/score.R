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

# helper functions for parsing out terms for factor interactions

.subgroups <- function(term, model, fn=as.formula(model), data)
{
    f <- attr(terms(fn), 'factors')
    vars <- setdiff(names(which(rowSums(f[,f[term,]>0])>0)), term)
    if (missing(model)) {
        fn <- function(x) levels(data[,x])
    } else if (!is.null(model$xlevels)) {
        fn <- function(x) {
            g <- model$xlevels[[x]]
            if (is.null(g)) c('FALSE','TRUE') else g
        }
    } else {
        fn <- function(x) dimnames(model$contr[[x]])[[1]]
    }
    sapply(vars, fn, simplify=FALSE)
}

.subgroup.effects <- function(model, term)
{
    grp.coef <- function(n)
    {
        contr <- lapply(names(grps), function(i)
                        contr.treatment(grps[[i]], unclass(gr[n,i])))
        names(contr) <- names(grps)
        u <- update(model, formula=formula(model),
                    data=model.frame(model), contrasts=contr)
        d <- merge(model.frame(model), gr[n,,FALSE])
        d <- subset(d, complete.cases(d))
        if (length(table(d[,term])) < 2)
            c(NA,NA,NA,NA)
        else coef(summary(u))[term,]
    }
    grps <- .subgroups(term, model=model)
    gr <- do.call('expand.grid', grps)
    d <- data.frame(gr, t(sapply(1:nrow(gr), grp.coef)))
    names(d) <- c(names(gr),'effect','stderr','t.value','pvalue')
    rn <- lapply(names(grps), function(x) paste(x, gr[,x], sep=''))
    structure(d, row.names=sub(' ',':',do.call('paste',rn)))
}

#-----------------------------------------------------------------------

score.kruskal <-
function(formula, data, ploidy,
         mode=c('general','recessive','dominant'))
{
    data$genotype <- .recode.gt(data$genotype, match.arg(mode))
    t <- kruskal.test(formula, data)
    keep.attr(data.frame(pvalue=t$p.value), fit='kruskal')
}

score.jt <- function(formula, data, ploidy, ...)
{

    jt <- jt.test(data[,.lhs.formula(formula)],
                  as.ordered(data$genotype), ...)
    keep.attr(with(jt, data.frame(
        pvalue=jt$p.value,
        effect=6*(jt$statistic-jt$EH)/jt$N^2,
        stderr=6*sqrt(jt$VH)/jt$N^2
    )), fit='jt')
}

score.chisq.2x2 <- function(formula, data, ploidy=NA)
{
    pt <- data[,.lhs.formula(formula)]
    gt <- factor(data$genotype,levels=0:2)
    if (ploidy %in% c('A',NA)) {
        t <- table(pt, gt)
        n <- 2*t[,c(1,3)] + t[,2]
    } else if (ploidy == 'X') {
        t <- table(pt, gt, data$gender)
        n <- t[,c(1,3),'M'] + 2*t[,c(1,3),'F'] + t[,2,'F']
    } else if (ploidy %in% c('Y','M')) {
        t <- table(pt, gt)
        n <- t[,c(1,3)]
    }
    stopifnot(nrow(n) == 2, ncol(n) == 2)
    t <- chisq.test(n, correct=FALSE)
    lo <- log(n[1,1]*n[2,2]/(n[2,1]*n[1,2]))
    se <- sqrt(sum(1/n))
    #phi2 <- t$statistic / sum(n)
    keep.attr(data.frame(
        pvalue=t$p.value, effect=lo, stderr=se
    ), fit='chisq.2x2')
}

score.chisq <-
function(formula, data, ploidy,
         mode=c('general','recessive','dominant'))
{
    data$genotype <- .recode.gt(data$genotype, match.arg(mode))
    l <- factor(data[,.lhs.formula(formula)])
    r <- factor(data[,.rhs.formula(formula)])
    n <- table(l,r)
    stopifnot(nrow(n) > 1, ncol(n) > 1)
    t <- chisq.test(n, correct=FALSE)
    # square of Cramer's V (reduces to phi^2 for 2x2 tables)
    #v2 <- unname(t$statistic) / (sum(n) * (min(nrow(n),ncol(n))-1))
    if (nrow(n) == 2 && ncol(n) == 2) {
        lo <- log(n[1,1]*n[2,2]/(n[2,1]*n[1,2]))
        se <- sqrt(sum(1/n))
    } else {
        lo <- NA
        se <- NA
    }
    keep.attr(data.frame(
        pvalue=t$p.value, effect=lo, stderr=se
    ), fit='chisq')
}

score.fisher <-
function(formula, data, ploidy,
         mode=c('general','recessive','dominant'), ...)
{
    data$genotype <- .recode.gt(data$genotype, match.arg(mode))
    l <- factor(data[,.lhs.formula(formula)])
    r <- factor(data[,.rhs.formula(formula)])
    t <- fisher.test(table(l,r), ...)
    or <- ifelse(is.null(t$estimate), NA, t$estimate)
    keep.attr(data.frame(
        pvalue=t$p.value, effect=log(or)
    ), fit='fisher')
}

score.trend <- function(formula, data, ploidy)
{
    x <- table(data$genotype, data[,.lhs.formula(formula)])
    stopifnot(ncol(x)==2)
    n <- x[,1] + x[,2] ; p <- x[,2]/n
    p0 <- sum(x[,2])/sum(n)
    score <- 1:length(p)
    m <- lm(p~score, weights=n/p0/(1-p0))
    a <- anova(m)
    x <- unname(coef(summary(m))['score',])
    cs <- a['score','Sum Sq']
    pv <- pchisq(cs, 1, lower.tail=FALSE)
    keep.attr(data.frame(
        pvalue=pv, effect=x[1], stderr=x[2]
    ), fit='trend')
}

#-----------------------------------------------------------------------

score.lm <-
function(formula, data, ploidy, drop='genotype',
         mode=c('additive','recessive','dominant','general'))
{
    data$genotype <- .recode.gt(data$genotype, match.arg(mode))
    m <- lm(formula, data)
    nf <- .null.model(formula, term=drop)
    pvalue <- anova(update(m,nf), m)[2,6]
    x <- coef(summary(m))
    x <- if (drop %in% row.names(x)) unname(x[drop,]) else c(NA,NA)
    keep.attr(data.frame(
        pvalue, effect=x[1], stderr=x[2]
    ), fit='lm')
}

score.lm.general <- function(formula, data, ploidy)
{
    if (length(unique(data$genotype)) < 3) {
        r <- score.lm(formula, data, ploidy)
        return(keep.attr(cbind(term=NA, r), fit='lm'))
    }
    data$genotype <- factor(data$genotype)
    contrasts <- list(genotype=matrix(c(0,1,2,0,1,0),3,2))
    m <- lm(formula, data, contrasts=contrasts)
    nf <- .null.model(formula)
    pvalue <- anova(update(m,nf,contrasts=), m)[2,6]
    x <- unname(coef(summary(m))[-1,])
    keep.attr(rbind(
        data.frame(term=NA, pvalue, effect=NA, stderr=NA),
        data.frame(term=c('additive','dominance'),
                   pvalue=x[,4], effect=x[,1], stderr=x[,2])
    ), fit='lm')
}

score.lm.groups <-
function(formula, data, ploidy, quick=FALSE,
         mode=c('additive','recessive','dominant'))
{
    nested <- function(term, new.form)
    {
        a <- anova(update(m,new.form), m)
        data.frame(term, pvalue=a[2,6], effect=NA, stderr=NA)
    }
    data$genotype <- .recode.gt(data$genotype, match.arg(mode))
    m <- lm(formula, data)
    if (quick)
        return(nested('Full vs Null Model', .null.model(formula)))
    e <- .subgroup.effects(m, 'genotype')
    if (!any(!is.na(e$pvalue))) return(NULL)
    term <- paste(row.names(e), 'genotype', sep=':')
    keep.attr(rbind(
        nested(NA, .null.model(formula)),
        nested('interaction', .flat.model(formula)),
        data.frame(
            term, pvalue=e$pvalue, effect=e$effect, stderr=e$stderr
        )
    ), fit='lm')
}

#-----------------------------------------------------------------------

score.glm <-
function(formula, data, ploidy, drop='genotype', test=c('LR','Wald','Rao'),
         mode=c('additive','recessive','dominant','general'))
{
    data$genotype <- .recode.gt(data$genotype, match.arg(mode))
    m <- glm(formula, binomial, data)
    x <- coef(summary(m))
    test <- match.arg(test)
    if (test == 'LR') {
        nf <- .null.model(formula, term=drop)
        pvalue <- anova(update(m,nf), m, test='Chisq')[2,5]
    } else if (test == 'Rao') {
        require(statmod)
        nf <- .null.model(formula, term=drop)
        pvalue <- 2*pnorm(-abs(glm.scoretest(update(m,nf), data[,drop])))
    } else {
        pvalue <- x[drop,4]
    }
    x <- if (drop %in% row.names(x)) unname(x[drop,]) else c(NA,NA)
    keep.attr(data.frame(
        pvalue, effect=x[1], stderr=x[2]
    ), fit='glm')
}

score.glm.general <- function(formula, data, ploidy)
{
    if (length(unique(data$genotype)) < 3) {
        r <- score.glm(formula, data, ploidy)
        return(keep.attr(cbind(term=NA, r), fit='glm'))
    }
    data$genotype <- factor(data$genotype, levels=0:2)
    contrasts <- list(genotype=matrix(c(0,1,2,0,1,0),3,2))
    m <- glm(formula, binomial, data, contrasts=contrasts)
    nf <- .null.model(formula)
    pvalue <- anova(update(m,nf,contrasts=), m, test='Chisq')[2,5]
    x <- unname(coef(summary(m))[-1,])
    keep.attr(rbind(
        data.frame(term=NA, pvalue, effect=NA, stderr=NA),
        data.frame(term=c('additive','dominance'),
                   pvalue=x[,4], effect=x[,1], stderr=x[,2])
    ), fit='glm')
}

score.glm.groups <-
function(formula, data, ploidy, quick=FALSE,
         mode=c('additive','recessive','dominant'))
{
    nested <- function(term, new.form)
    {
        a <- anova(update(m,new.form), m, test='Chisq')
        data.frame(term, pvalue=a[2,5], effect=NA, stderr=NA)
    }
    data$genotype <- .recode.gt(data$genotype, match.arg(mode))
    m <- glm(formula, binomial, data)
    if (quick)
        return(nested('Full vs Null Model', .null.model(formula)))
    e <- .subgroup.effects(m, 'genotype')
    if (!any(!is.na(e$pvalue))) return(NULL)
    term <- paste(row.names(e), 'genotype', sep=':')
    keep.attr(rbind(
        nested(NA, .null.model(formula)),
        nested('interaction', .flat.model(formula)),
        data.frame(
            term, pvalue=e$pvalue, effect=e$effect, stderr=e$stderr
        )
    ), fit='glm')
}

#-----------------------------------------------------------------------

.gls.effects <- function(model, term)
{
    grps <- .subgroups(term, model=model)
    gr <- do.call('expand.grid', grps)
    rn <- lapply(names(grps), function(x) paste(x, gr[,x], sep=''))
    cn <- names(coef(model))
    L <- matrix(0,nrow(gr),length(cn))
    L[,cn==term] <- 1
    for (nr in 1:nrow(gr)) {
        off <- which(regexpr(paste(term,'$',sep=''),cn) < 0)
        for (nc in 1:ncol(gr)) {
            a <- grep(paste('\\b',names(gr)[nc],sep=''),cn)
            b <- grep(paste('\\b',rn[[nc]][nr],'\\b',sep=''),cn)
            off <- union(off, setdiff(a,b))
            L[nr,b] <- 1
        }
        L[nr,off] <- 0
    }
    e <- do.call('rbind',lapply(1:nrow(L), function(x)
                                anova(model, L=L[x,])))
    v <- as.vector(L %*% coef(model))
    data.frame(
        pvalue=e[,3], effect=v, stderr=abs(v/e[,2]^0.5),
        row.names=sub(' ',':',do.call('paste',rn))
    )
}

score.gls.groups <- function(formula, data, ploidy, ...)
{
    fix <- (sapply(data,class)=='logical')
    data[,fix] <- lapply(data[,fix,drop=FALSE], as.factor)
    do.gls <- function(formula)
        eval(parse(text=paste('gls(',.as.text(formula),',data,...)')))
    nested <- function(term, new.form)
    {
        a <- anova(do.gls(new.form), m)
        data.frame(term, pvalue=a[2,9], effect=NA, stderr=NA)
    }
    m <- do.gls(formula)
    e <- .gls.effects(m, 'genotype')
    if (!any(!is.na(e$pvalue))) return(NULL)
    term <- paste(row.names(e), 'genotype', sep=':')
    keep.attr(rbind(
        nested(NA, .null.model(formula)),
        nested('interaction', .flat.model(formula)),
        data.frame(
            term, pvalue=e$pvalue,
            effect=e$estimate, stderr=e$std.error
        )
    ), fit='gls')
}

