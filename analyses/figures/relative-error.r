## plot of relative error 
data(mom_rel_error)

library(RColorBrewer)
textwidth <- 7 # 426 pts to inches
mtext_cex <- 1.2
opar <- par(no.readonly=TRUE)
phi <- (1 + sqrt(5))/2
title_col='grey12'
dashed_gray <- 'gray58'

pdf('supp-relative-error.pdf', width=textwidth, height=0.8 * textwidth/phi)
brewercols <- brewer.pal(5,"Dark2")
layout(matrix(c(1,2), ncol=2, byrow=TRUE))
par(oma=c(0, 0, 0, 0), mar=c(3, 3, 1, 2))
ymax <- 0.8

plot.default(as.factor(mom_rel_error$sample_size), mom_rel_error$N_error, 
               pch=19, type='n', bty='n', axes=FALSE,
               ylab='median relative error', xlab='sample size',
               line=1.8,
               ylim=c(0, ymax),
               col=brewercols[mom_rel_error$Va])

for (va in unique(mom_rel_error$Va)) {
  rel <- mom_rel_error[mom_rel_error$Va == va,]
  points.default(as.factor(rel$sample_size), rel$N_error, 
               pch=19, type='b', bty='n', axes=FALSE,
               col=brewercols[rel$Va])
}



axes_col <- 'gray12'
axis_cex <- 1*0.7
ss <- sort(unique(mom_rel_error$sample_size))
axis(1, at=1:4, labels=ss, 
     padj=-1.5,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(2, at=seq(0, ymax, length.out=5),
     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)

legend(2.5, 0.72, inset=0, legend=unique(mom_rel_error$Va), fill=brewercols,
       title=latex2exp:::TeX('additive genetic variation'),
       cex=0.6,
       ncol=2,
       xpd=TRUE,
       bty='n', text.col=title_col,
       border=0)

mtext(latex2exp:::TeX('$\\widehat{N}$'), 3, line=-1, cex=mtext_cex)



plot.default(as.factor(mom_rel_error$sample_size), mom_rel_error$va_error, 
               pch=19, type='n', bty='n', axes=FALSE,
               line=1.8,
               ylim=c(0, 10),
               ylab='median relative error', xlab='sample size',
               col=brewercols[mom_rel_error$Va])

for (va in unique(mom_rel_error$Va)) {
  rel <- mom_rel_error[mom_rel_error$Va == va,]
  points.default(as.factor(rel$sample_size), rel$va_error, 
               pch=19, type='b', bty='n', axes=FALSE,
               col=brewercols[rel$Va])
}

axis(1, at=1:4, labels=ss, 
     padj=-1.5,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(2, at=seq(0, 10, 2), 
     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)

legend(2.5, 9, inset=0, legend=unique(mom_rel_error$Va), fill=brewercols,
       title=latex2exp:::TeX('additive genetic variation'),
       cex=0.6,
       ncol=2,
       xpd=TRUE,
       bty='n', text.col=title_col,
       border=0)

mtext(latex2exp:::TeX('$\\widehat{V_A(1)}$'), 3, line=-1, cex=mtext_cex)


dev.off()

par(opar)
