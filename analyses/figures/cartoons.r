# 
library(wesanderson)
library(diagram)
library(tikzDevice)

textwidth <- 5.91667  # 426 pts to inches

ssf <- function(x) sprintf('\\textsf{%s}', x)

stack <- TRUE

if (!interactive())  {
  if (!stack) {
    tikz('figure-1.tex', width=6.5, height=3.2, standAlone=TRUE)
  } else {
    tikz('figure-1.tex', width=6.5/2, height=3.2*1.5, standAlone=TRUE)
  }
}
#pdf('figure-1.pdf', width=6.5/2, height=2*3.2)
#pdf('figure-1.pdf', width=6.5, height=3.2)

# loadfonts()

cols <- wes_palette('Darjeeling1')
axs_col <- 'gray38'
axs_lwd <- 2.
dashed_col <- 'gray12'
vertlty <- 1

### FIGURE A
opar <- par(no.readonly=TRUE)
par(mar=c(3.1, 3.1, 2.1, 1.1))

logistic <- function(p0, t, s=0.1) p0/(p0 + (1-p0)*exp(-s*t))

if (!stack) {
  layout(matrix(c(1, 2), nrow=1), widths=c(0.5, 0.3))
} else {
  layout(matrix(c(1, 2), nrow=2), heights=c(0.7, 0.4))
}

gens <- seq(0, 50, length.out=5)
covs <- logistic(0.5, gens, 0.08)
plot(gens, covs, type='b',  lwd=3, col=cols[5],
     axes=FALSE, ylab='', xlab='', ylim=c(0.4, 1))

axis(1, at=gens, labels=ssf(seq_along(gens)-1), lwd=axs_lwd, 
     cex.axis=0.8,
     col=axs_col, tck=0.01, lend=2, padj=-2)
axis(2, at=seq(0.4, 1, length.out=4), labels=ssf(seq(0.4, 1, length.out=4)-0.1), 
     lwd=axs_lwd, hadj=0.3, cex.axis=0.8,
     col=axs_col, tck=0.01, las=1, lend=2)
mtext(ssf("generation"), 1, col=dashed_col, line=1)
mtext(ssf("frequency"), 2, col=dashed_col, line=1.4)


t <- 1
segments(gens[t], covs[t], x1=gens[t+1], y1=covs[t], lty=3, col=dashed_col)
segments(gens[t+1], covs[t+1], x1=gens[t+1], y1=covs[t], lty=vertlty, 
         col=dashed_col)

t <- 3
segments(gens[t], covs[t], x1=gens[t+1], y1=covs[t], lty=3, col=dashed_col)
segments(gens[t+1], covs[t+1], x1=gens[t+1], y1=covs[t], lty=vertlty, 
         col=dashed_col)

points(gens, covs, pch=19, lwd=3, col=cols[5])

if (TRUE) {
  t <- 1
  segments(gens[t], 1.5-covs[t], x1=gens[t+1], y1=1.5-covs[t], 
           lty=3, col=dashed_col)
  segments(gens[t+1], 1.5-covs[t+1], x1=gens[t+1], y1=1.5-covs[t], 
           lty=vertlty, col=dashed_col)

  t <- 3
  segments(gens[t], 1.5-covs[t], x1=gens[t+1], y1=1.5-covs[t], 
           lty=3, col=dashed_col)
  segments(gens[t+1], 1.5-covs[t+1], x1=gens[t+1], y1=1.5-covs[t], 
           lty=vertlty, col=dashed_col)
  points(gens, 1.5-covs, pch=19, lwd=3, col=cols[4], type='b')
  text(gens[2]+4, 1.5-covs[1] - 0.09, '$\\Delta p_0$', col=cols[4])
  text(gens[4]+4, 1.5-covs[3] - 0.03, '$\\Delta p_2$', col=cols[4])
}



# text(3.2, 0.75, latex2exp::TeX('$\\Delta p_2$'))
# text(1.2, 0.5, latex2exp::TeX('$\\Delta p_1$'))
text(gens[2]+4, covs[1] + 0.1, '$\\Delta p_0$', col=cols[5])
text(gens[4]+4, covs[3] + 0.03, '$\\Delta p_2$', col=cols[5])
# text(gens[3]+4, covs[1] + 0.3, 'cov$(\\Delta p_0, \\Delta p_2) > 0$')
legend(10, 0.5, ssf(c('high fitness background', 'low fitness background')), 
       border=0, x.intersp=0.4,
       fill=c(cols[5], cols[4]), lty=0, box.lty=0, cex=0.9, text.col=dashed_col)

text(-5, 1.06, ssf("A"), las=1, cex=1.4, font=2, col=dashed_col, xpd=TRUE)

### FIGURE B
par(mar=c(1, 1, 0, 1))
plot.new()
plot.window(c(0, 10), c(0, 3))
text(0, 2.5, ssf("B"), las=1, cex=1.4, font=2, col=dashed_col)

arrow_ld <- 1.4
arrow_lty <- 2
lwd <- 0.13/2

selsites <- c(0.5, 3, 9, 7, 2)

rect(0, 1-0.1, 10, 1.1-0.1, col='gray28', lty=0)


crv <- 0.12
for (s in selsites) {
  curve <- -crv
  if (s > 10/2)
    curve <- crv
  curvedarrow(c(s, 1), c(10/2, 1), curve=curve, lwd=arrow_ld, 
              arr.width=0.15,
              lcol=axs_col, arr.pos=0.5, lty=arrow_lty)
}


# selected sites
for (s in selsites) {
  rect(s-lwd, 1-0.1, s+lwd, 1+0.08, col=cols[3], lty=0)
}

rect(10/2 - lwd, 1-0.1, 10/2 + lwd, 1+0.08, col=cols[5], lty=0)

text(10/2, 1-0.6, '$R$', cex=2, 
     font=3,
     col=rgb(47,47,47, maxColorValue=255))

segments(0, 1 - 0.6, 10/2 - 0.6, 1 - 0.6)
segments(10/2 + 0.6, 1 - 0.6, 10, 1 - 0.6)
segments(0, 1-0.6-0.05, 0, 1-0.6+0.05)
segments(10, 1-0.6-0.05, 10, 1-0.6+0.05)

legend(1.5, 2.4, ssf(c('focal neutral site', 'selected sites')), 
       border=0, x.intersp=0.4,
       ncol=2,
       fill=c(cols[5], cols[3]), lty=0, box.lty=0, cex=0.9, text.col=dashed_col)


par(opar)

if (!interactive()) dev.off()
