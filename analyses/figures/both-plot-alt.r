library(RColorBrewer)

textwidth <- 7 # 426 pts to inches
mtext_cex  <- 0.8
opar <- par(no.readonly=TRUE)
phi <- (1 + sqrt(5))/2
title_col='grey12'
dashed_gray <- 'gray58'


pretty_log <- function(x, use_one=TRUE) {
  if (use_one && x==0) 
    return('1')
  else 
    return(latex2exp::TeX(sprintf('$10^{%d}$', x)))
}


## Estimation Figure -- both Va and N estimation, with N boxplot grouped by Va
pdf('mom-fits-both-alt.pdf', width=textwidth, height=textwidth/2 * 0.9)

#data(mom_fitsd_finite)
mom_fitsd <- mom_fitsd_finite
mom_fitsd_subset <- mom_fitsd %>% 
  filter(Va %in% c(0.001, 0.005, 0.01, 0.05, 0.1)) #%>%
  #filter(genlen %in% c(0.1, 0.5, 1.5))

# Layout
layout(matrix(c(1, 2, 3, 4), ncol=2, byrow=TRUE), heights=c(0.1, 0.90))
opar <- par(no.readonly=TRUE)
par(mar=c(0,0,0,0), xpd=TRUE)
plot.new()
text(0.1, 0.3, "A", cex=1.5)

plot.new()
text(0.1, 0.3, "B", cex=1.5)

wescols <- wesanderson::wes_palette('Darjeeling1')[c(4, 1, 2, 3, 5)]
brewercols <- brewer.pal(5,"Dark2")

mom_fitsd$genlen <- droplevels(mom_fitsd$genlen)

par(mar=c(4.1, 3.1, 0.1, 2.1))
with(mom_fitsd[mom_fitsd$va0_est > 1e-4 & mom_fitsd$rep %in% 1:10, ],
     plot(emp_va, va0_est, col=wescols[as.factor(genlen)],
          xlab='', ylab='',
          pch=19, cex=0.5, axes=FALSE,
          ylim=c(1e-4, 10),
          xlim=c(1e-4, 1),
          lwd=0, log='xy'))
#abline(a=0, b=1, col=dashed_gray, lty=2)
segments(1e-4, 1e-4, 1, 1, col=dashed_gray, lty=2)
axes_col <- 'gray12'
axis_cex <- 1*0.7
axis(1, at=10^(-(4:0)), labels=pretty_log(-(4:0)), 
     padj=-1.5,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(2, at=10^(c(-(4:0), 1)), labels=pretty_log(c(-(4:0), 1)), 
     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)
mtext(latex2exp:::TeX('empirical $V_A(1)$'), 1, line=1.3, cex=mtext_cex)
mtext(latex2exp:::TeX('estimated $V_A(1)$'), 2, line=1.5, cex=mtext_cex)
legend(1e-4, 5, inset=0, legend=unique(mom_fitsd$genlen), fill=wescols,
           title=latex2exp:::TeX('level of\nrecombination (Morgans)'),
           cex=0.6,
           ncol=2,
           bty='n', text.col=title_col,
           border=0)
#text(1e-4, 11, 'A', line=2, cex=2, xpd=TRUE)

ymax <- 3e3
boxplot(Nest ~ Va + genlen, mom_fitsd_subset,
        axes=FALSE, ylim=c(0,ymax), outline=FALSE,
        col=brewercols,
        border=brewercols,
        pars = list(boxcol = "transparent", medlty = "blank", 
                    medpch=16, whisklty = c(1, 1),
                    medcex = 0.5, outcex = 0, staplelty = "blank"))

#abline(h=1e3, col=dashed_gray, lty=2)
y <- mom_fitsd_subset %>% group_by(Va, genlen) %>% 
  summarize(Ne_est=mean(Ne_est, na.rm=TRUE)) %>%
  arrange(genlen, Va) %>% pull(Ne_est)
x <- 1:length(y)
segments(0, 1e3, length(x), 1e3, col=dashed_gray, lty=2)
points(x, y, pch='-', col=brewercols)
hadj  <- 2.3
axis(1, at=c(1, 5), labels=c("", "0.01"),
     hadj=hadj-0.6,
     padj=-2.4,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(1, at=c(6, 10), 
     labels=c("", "0.1"), hadj=hadj,
     padj=-2.4,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(1, at=c(11, 15), 
     labels=c("", "0.5"), hadj=hadj,
     padj=-2.4,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(1, at=c(16, 20), 
     labels=c("", "1.5"), hadj=hadj,
     padj=-2.4,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)

axis(2, at=seq(0, ymax, 500), 
     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)
mtext(latex2exp:::TeX('level of recombination (Morgans)'), 
      1, line=1.2, cex=mtext_cex)
mtext(latex2exp:::TeX('estimated N'), 2, line=2, cex=mtext_cex)

legend(11, 1500, inset=0, legend=unique(mom_fitsd_subset$Va), fill=brewercols,
       title=latex2exp:::TeX('additive genetic variation'),
       cex=0.6,
       ncol=3,
       xpd=TRUE,
       bty='n', text.col=title_col,
       border=0)


#boxplot(Nest ~ genlen, mom_fitsd, axes=FALSE, ylim=c(0, 2e3), outline=FALSE,
#        pars = list(boxcol = "transparent", medlty = "blank", 
#                    medpch=16, whisklty = c(1, 1),
#                    medcex = 0.7,  outcex = 0, staplelty = "blank"))
##abline(h=1e3, col=dashed_gray, lty=2)
#segments(0.3, 1e3, 4.4, 1e3, col=dashed_gray, lty=2)
#axis(1, at=unique(mom_fitsd$genlen), labels=levels(mom_fitsd$genlen), 
#     padj=-2.4,
#     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#axis(2, at=c(0, 500, 1e3, 1.5e3, 2e3), 
#     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#mtext(latex2exp:::TeX('level of recombination (Morgans)'), 
#      1, line=1.2, cex=mtext_cex)
#mtext(latex2exp:::TeX('estimated N'), 2, line=2, cex=mtext_cex)
#text(0, 2010, 'B', line=2, cex=2, xpd=TRUE)
dev.off()


