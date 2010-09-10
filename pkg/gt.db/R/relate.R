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

gt.dist <-
function(gt1, gt2=gt1, operation='==', aggregator='sum', na.value=NA)
{
    .Call("do_gt_dist", gt1, gt2, operation, aggregator, na.value,
          PACKAGE='gt.db')
}

ibs.gt.data <- function(x, sample.mask=TRUE)
{
    g <- x$genotype
    if (!identical(sample.mask,TRUE))
        g <- mask.str(g, sample.mask)
    list("IBS=0"=gt.dist(g, operation='ibs==0'),
         "IBS=1"=gt.dist(g, operation='ibs==1'),
         "IBS=2"=gt.dist(g, operation='ibs==2'))
}

ibs.gt.dataset <- function(x, sample.mask=TRUE)
{
    sum.fn <- function(a, b)
    {
        if (is.null(a)) return(b)
        mapply('+', a, b, SIMPLIFY=FALSE)
    }
    apply.gt.dataset(x, ibs.gt.data, sum.fn, sample.mask, GT.DATA='x')
}

ibs <- function(x, sample.mask=TRUE) UseMethod('ibs')

grr <- function(x, sample.mask=TRUE)
{
    i <- ibs(x, sample.mask)
    n <- i[[1]]+i[[2]]+i[[3]]
    m <- (i[[2]]+2*i[[3]])/n
    v <- ((m-0)^2*i[[1]]+(m-1)^2*i[[2]]+(m-2)^2*i[[3]])/n
    list(mu=m, sd=sqrt(v))
}

ibd.gt.data <-
function(gt.data, binsz=1e6, min.snps=25, max.snps=50,
         min.gt=0.8, ibs.limit=2)
{
    nc <- max(nchar(gt.data$genotype))
    m0 <- m1 <- m2 <- matrix(0, nc, nc)
    if (any(gt.data$ploidy != 'A')) {
        gt.data <- subset(gt.data, ploidy=='A')
        warning('non-autosomal data will be excluded')
    }

    # remove outlier samples
    gt <- with(summary.gt.data(gt.data, by.sample=TRUE), AA+AB+BB)
    bs <- (gt < min.gt*max(gt))
    if (any(bs)) {
        cat("excluding samples with insufficient genotypes:\n")
        print(which(bs))
        m0[bs,] <- m1[bs,] <- m2[bs,] <- NA
        m0[,bs] <- m1[,bs] <- m2[,bs] <- NA
    }

    bin <- paste(gt.data$scaffold, trunc(gt.data$position/binsz))
    bin.ok <- names(which(table(bin) >= min.snps))
    cat("Distribution of SNP counts for selected bins:\n");
    print(summary(table(bin)[bin.ok]))
    stopifnot(length(bin.ok) > 0)

    n1 <- skip <- numeric(0)
    progress.bar(0, length(bin.ok))
    for (i in 1:length(bin.ok)) {
        g <- head(gt.data$genotype[bin==bin.ok[i]], n=max.snps)
        s <- gt.dist(g, operation='ibs', aggregator='min')
        s[is.na(s)] <- -1
        n1[i] <- sum(s==1, na.rm=TRUE)
        skip[i] <- (n1[i] > ibs.limit*median(n1))
        if (!skip[i]) {
            m0 <- m0 + (s==0)
            m1 <- m1 + (s==1)
            m2 <- m2 + (s==2)
        }
        progress.bar(i, length(bin.ok))
    }

    # fixup outlier list based on final median IBS=1
    last <- (n1 > ibs.limit*median(n1))
    for (i in which(xor(skip, last))) {
        s <- gt.dist(gt.data$genotype[bin==bin.ok[i]],
                     operation='ibs', aggregator='min')
        s[is.na(s)] <- -1
        x <- if (skip[i]) 1 else -1
        m0 <- m0 + x * (s==0)
        m1 <- m1 + x * (s==1)
        m2 <- m2 + x * (s==2)
    }
    cat("Excluded", sum(last), "of", length(bin.ok), "bins\n");

    n <- m0+m1+m2
    id <- sample.names(gt.data)
    dimnames(m1) <- dimnames(m2) <- list(id,id)
    list("IBD=1"=m1/n, "IBD=2"=m2/n)
}

ibd.dataset <-
function(dataset.name, binsz=1e6, part=1, parts=10,
         maf.min=0.1, hw.p.min=0.01, gt.rate.min=0.98, ...)
{
    g <- fetch.gt.data(dataset.name, part=part, parts=parts,
                       by='position', binsz=binsz)
    g <- subset(g, ploidy=='A')
    g <- with(summary.gt.data(g), subset(g, gt.rate>=gt.rate.min &
              pmin(freq.a,freq.b)>=maf.min & hw.p.value>=hw.p.min))
    ibd.gt.data(g, binsz=binsz, ...)
}

ibd.summary <- function(ibd, alpha=0.05)
{
    u <- ibd[[1]][lower.tri(ibd[[1]])]
    v <- ibd[[2]][lower.tri(ibd[[2]])]

    qv <- abs(qnorm(alpha/sum(!is.na(u))))
    grp <- (abs(scale(u)) > qv) + 2*(v>0.05)
    u2 <- ifelse(grp,NA,u)
    q2 <- mean(u2,na.rm=TRUE) + qv*sd(u2,na.rm=TRUE)
    grp <- grp + if.na(u2 > q2, 0)
    cat("Outlier cutoff for IBS = 1:", q2, "\n")
    table(factor(grp,levels=0:3, labels=c('None','IBD=1','IBD=2','Both')))
}

ibd.outliers <- function(ibd, alpha=0.05)
{
    u <- ibd[[1]][lower.tri(ibd[[1]])]
    v <- ibd[[2]][lower.tri(ibd[[2]])]
    qv <- abs(qnorm(alpha/sum(!is.na(u))))
    u[(abs(scale(u)) > qv) + 2*(v>0.05)] <- NA
    q2 <- mean(u,na.rm=TRUE) + qv*sd(u,na.rm=TRUE)

    w1 <- which(ibd[[1]] > q2 & ibd[[2]] < 0.05, TRUE)
    w2 <- which(ibd[[1]] < q2 & ibd[[2]] > 0.05, TRUE)
    w3 <- which(ibd[[1]] > q2 & ibd[[2]] > 0.05, TRUE)
    w <- rbind(cbind(w1,grp=rep(1,nrow(w1))),
               cbind(w2,grp=rep(2,nrow(w2))),
               cbind(w3,grp=rep(3,nrow(w3))))

    n <- rownames(ibd[[1]])
    w <- w[n[w[,1]] < n[w[,2]],]
    d <- data.frame(sample.name.1=n[w[,1]], sample.name.2=n[w[,2]],
                    ibd1=ibd[[1]][w[,-3]], ibd2=ibd[[2]][w[,-3]])
    d$grp <- factor(w[,3],levels=1:3,labels=c('IBD=1','IBD=2','Both'))
    d[,1] <- as.character(d[,1])
    d[,2] <- as.character(d[,2])
    d <- d[order(d$grp,d$sample.name.1,d$sample.name.2),]
    row.names(d) <- 1:nrow(d)
    d
}

ibd.plot <- function(ibd, jitter.amount, split=0.3,
                     panel.width=list(0.8,'npc'))
{
    u <- ibd[[1]][lower.tri(ibd[[1]])]
    v <- ibd[[2]][lower.tri(ibd[[2]])]

    qv <- abs(qnorm(0.05/sum(!is.na(u))))
    grp <- (abs(scale(u)) > qv) + 2*(v>0.05)
    u2 <- ifelse(grp,NA,u)
    q2 <- mean(u2,na.rm=TRUE) + qv*sd(u2,na.rm=TRUE)
    grp <- grp + if.na(u2 > q2, 0)

    s <- unique(c(sample(length(u), 10000), which(grp>0)))
    if (!missing(jitter.amount)) v <- jitter(v, amount=jitter.amount)
    par <- list(layout.heights=list(bottom.padding=0,xlab=0))
    p1 <- xyplot(v[s]~u[s], groups=grp[s], xlim=c(-0.05,1.05),
                 ylim=c(-0.05,1.05), xlab=NULL, ylab='IBD = 2', pch=1,
                 par.settings=par, scales=list(x=list(alternating=0)))
    sc <- list(y=list(at=seq(0,100,10), labels=format(seq(0,1,0.1))))
    par <- list(layout.heights=list(top.padding=0,main=0,axis.top=0))
    p2 <- histogram(~u, breaks=seq(0,1,0.01), type='density', scales=sc,
                    par.settings=par, xlim=c(-0.05,1.05), xlab='IBD = 1')
    print(p1, position=c(0.0,split,1.0,1.0),
          panel.width=panel.width, more=TRUE)
    print(p2, position=c(0.0,0.0,1.0,split), panel.width=panel.width)
}
