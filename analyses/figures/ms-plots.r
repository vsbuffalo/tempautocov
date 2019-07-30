# plots.r -- create all plots for the manuscript (other than cartoons)
devtools::load_all()
library(viridis)

pretty_log <- function(x, use_one=TRUE) {
  if (use_one && x==0) 
    return('1')
  else 
    return(latex2exp::TeX(sprintf('$10^{%d}$', x)))
}
textwidth <- 7 # 426 pts to inches
mtext_cex  <- 0.8
opar <- par(no.readonly=TRUE)
phi <- (1 + sqrt(5))/2
title_col='grey12'
legend_cex <- 1.1
axs_col <- 'gray12'

## OLD VERSION
# ## Figure 2 
# data(predf)
# data(pd5f)
# pd5f$L <- droplevels(pd5f$L)

# pdf('sim-pred-covs-varyl.pdf', width=textwidth, height=textwidth/phi)

# panel_plot(pd5f, predf, gen, cov, Va, genlen, L, 
#            row_lab='$V_a = $', col_lab = '$R = $', 
#            xlab='generation', 
#            y_axis_line=3,
#            x_axis_line=3,
#            pch_cex=0.8,
#            point_alpha=0.6,
#            legend_prop=0.7,
#            axis_lwd=1.2,
#            ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 
#            legend_title='number\nof loci')

# dev.off()

## Figure 2 from different starting generation
data(pred13f)
data(pd13f)
pd13f$L <- droplevels(pd13f$L)

pdf('sim-pred13-covs-varyl.pdf', width=textwidth, height=textwidth/phi)

panel_plot(pd13f, pred13f, gen, cov, Va, genlen, L, 
           row_lab='$V_A = $', col_lab = '$R = $', 
           xlab='generation', 
           y_axis_line=3,
           x_axis_line=3,
           pch_cex=0.8,
           point_alpha=0.6,
           grp_col='gray40',
           legend_prop=0.1,
           axis_lwd=1.2,
           ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 
           legend_title='number\nof loci')

dev.off()


## Figure 2 with varying N (appendix)
data(predf_varyn)
data(pd5f_varyn)

predf_varyn <- predf_varyn %>% filter(L == 500, gen > 5, type=='pred_zvar') %>% mutate(N=as.factor(N))
pd5f_varyn <- pd5f_varyn %>% filter(L == 500) %>% mutate(N=as.factor(N))

pdf('sim-pred-covs-varyn.pdf', width=textwidth, height=textwidth/phi)

panel_plot(pd5f_varyn, predf_varyn, gen, cov, Va, genlen, N, 
           row_lab='$V_A = $', col_lab = '$R = $', 
           xlab='generation', 
           y_axis_line=3,
           x_axis_line=3,
           pch_cex=0.8,
           point_alpha=0.6,
           legend_prop=0.7,
           axis_lwd=1.2,
           ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 
           legend_title='pop size')
dev.off()


## Figure 3 alt
data(cumcov_plot_all_df)
pdf('cummulative-cov-var-all.pdf', width=textwidth, height=textwidth/1.9, 
     pointsize=12)
cumcov_panels(cumcov_plot_all_df, TRUE,  ylab="var$(p_{15} - p_{5})$", 
	     lx=0.5/2, ly=0.5, lwd=1.2, ymax=0.5, basic=FALSE)
dev.off()
 
## Figure 3
data(cumcov_plot_df)
pdf('cummulative-cov-var.pdf', width=textwidth, height=textwidth/1.9, pointsize=12)
cumcov_panels(cumcov_plot_df, TRUE,  ylab="var$(p_{15} - p_{5})$", 
             lx=1.8, ly=0.57, lwd=1.2, ymax=0.5)
dev.off()

## Figure 2 ALT -- This is currently the one in the main text
data(predf_ssh)
# NOTE: we use pd5falt for the data, as this has all L
#data(pd5f_ssh)

pdf('sim-pred-covs-varyl-alt.pdf', width=textwidth, 
    height=textwidth/1.9, pointsize=12)
# we only plot points not too far below zero (this happens
# due to random fixations)

pd5falt <- pd5f %>% filter(cov > -0.001) %>% mutate(L=droplevels(L))
predf_ssh_alt <- predf_ssh
#predf_ssh_alt$var_type <- predf_ssh_alt$var_type %>% 
  #fct_relevel('neutral SSH', 'additive genetic', 'additive genic')

#levels(predf_ssh_alt$var_type) <- c('$V_a$', '$V_{a,ssh}$', '$V_A$')
levels(predf_ssh_alt$var_type) <- c('additive genic', 'neutral SSH', 'additive genetic')

#<- recode_factor(predf_ssh_alt$var_type, neut_ssh_proxy='neut. SSH', zvar='$V_A$', genic_va='$V_a$')

panel_plot2(pd5falt, predf_ssh_alt, gen, cov, Va, genlen, L, var_type,
           row_lab='$V_A = $', col_lab = '$R = $', 
           xlab='generation', 
           y_axis_line=3,
           x_axis_line=3,
           pch_cex=0.8,
           lwd=1.7,
           point_alpha=0.6,
           line_alpha=0.8,
           legend_prop=1.3,
           lg1y=0.5,
           axis_lwd=1.2,
           ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 
           legend_title1='number of loci', 
           legend_title2='variance')

dev.off()


## Figure 2, VAs over orders of magnitude (currently. supp.)
data(predf_ssh2)
data(pd5f_ssh2)
# NOTE: we use pd5falt for the data, as this has all L
pdf('sim-pred-covs-varyl-va-oom.pdf', width=textwidth, 
    height=textwidth/1.9, pointsize=12)
# we only plot points not too far below zero (this happens
# due to random fixations)

predf_ssh_alt2 <- predf_ssh2
#predf_ssh_alt$var_type <- predf_ssh_alt$var_type %>% 
  #fct_relevel('neutral SSH', 'additive genetic', 'additive genic')

#levels(predf_ssh_alt$var_type) <- c('$V_a$', '$V_{a,ssh}$', '$V_A$')
levels(predf_ssh_alt2$var_type) <- c('additive genic', 'neutral SSH', 'additive genetic')

#<- recode_factor(predf_ssh_alt2$var_type, neut_ssh_proxy='neut. SSH', zvar='$V_A$', genic_va='$V_a$')

panel_plot2(pd5f_ssh2, predf_ssh_alt2, gen, cov, Va, genlen, L, var_type,
           row_lab='$V_A = $', col_lab = '$R = $', 
           xlab='generation', 
           y_axis_line=3,
           x_axis_line=3,
           pch_cex=0.8,
           lwd=1.7,
           point_alpha=0.6,
           line_alpha=0.8,
           legend_prop=1.2,
           lg1y=0.5,
           axis_lwd=1.2,
           ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 
           legend_title1='number of loci', 
           legend_title2='variance')

dev.off()



### Figure 2 ALT -- no ssh
#data(predf_ssh)
#data(pd5f_ssh)

#pdf('sim-pred-covs-varyl-alt2.pdf', width=textwidth, 
#    height=textwidth/1.9, pointsize=12)
## we only plot points not too far below zero (this happens
## due to random fixations)

#pd5falt <- pd5f %>% filter(cov > -0.001)
#predf_ssh_alt <- predf_ssh %>% filter(var_type != 'neutral SSH')
##predf_ssh_alt$var_type <- predf_ssh_alt$var_type %>% 
#  #fct_relevel('neutral SSH', 'additive genetic', 'additive genic')

#predf_ssh_alt$var_type <- recode_factor(predf_ssh_alt$var_type, neut_ssh_proxy='neut. SSH', zvar='$V_A$', genic_va='$V_a$')

#levels(predf_ssh_alt$var_type) <- droplevels(predf_ssh_alt$var_type)
#levels(predf_ssh_alt$var_type) <- c('$V_A$', '$V_a$')

#panel_plot2(pd5falt, predf_ssh_alt, gen, cov, Va, genlen, L, var_type,
#           row_lab='$V_A = $', col_lab = '$R = $', 
#           xlab='generation', 
#           y_axis_line=3,
#           x_axis_line=3,
#           pch_cex=0.8,
#           lwd=1.7,
#           point_alpha=0.6,
#           line_alpha=0.8,
#           legend_prop=0.99,
#           lg1y=0.5,
#           axis_lwd=1.2,
#           ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 
#           legend_title1='number of loci', 
#           legend_title2='variance type')

#dev.off()


## Figure 4 -- Va/R plot
data(pd5_vard)

jitter_log10 <- function(x, scaler=0.05) x + rnorm(length(x), 0, x*scaler)


pdf('va-r-cov.pdf', width=textwidth/2, height=textwidth/2)

par(opar)
par(oma=c(0, 0, 0, 0), mar=c(2, 4.5, 1, 1))
pd5_dfs <- pd5_vard %>% filter(t1 < 30, va_r<=10, va_r>1e-3) 
n <- length(unique(pd5_dfs$gen))
cols <- viridis(n)
#cols <- wesanderson::wes_palette("Zissou1", n, type = "continuous")

## new version of the plot uses a subset of gens for clarity
KEEP_GENS <- seq(6, 24, 4)
pd5_dfs <- pd5_dfs %>% filter(gen %in% KEEP_GENS)

pd5_vard_lines <- pd5_vard %>% filter(gen %in% KEEP_GENS, is.finite(va_r), is.finite(cov)) %>% 
  filter(t1 < 30) %>% 
  group_by(gen) %>% 
  #summarize(cov=mean(cov)) %>% 
  group_by(gen) %>% nest()


with(pd5_dfs, 
     plot(jitter_log10(va_r), cov, col=adjustcolor(cols[gen], alpha.f=0.6), 
          log='x', 
          xlab='', ylab='',
          pch=19, cex=0.9, axes=FALSE,
          lwd=0,
          xlim=c(1e-3, 1e1),
          ylim=c(-0.001, 0.008)))
segments(1e-3, 0, 1e1, 0, col=dashed_gray, lty=2)
axes_col <- 'gray12'
axis_cex <- 1*0.7
axis(1, at=10^c(-3:1), 
     labels=sapply(-3:1, pretty_log, use_one=FALSE),
     padj=-1.2,
     tck=0.01, line=-1.4, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(2, tck=0.018, hadj=0.7, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)


for (i in seq_along(pd5_vard_lines$data)) {
  dat <- pd5_vard_lines$data[[i]]
  fit <- with(dat, loess(cov ~ va_r, span=0.9))
  va_r <- seq(1e-3, 1e1, length.out=100)
  lines(va_r, predict(fit, data.frame(va_r=va_r)), col=cols[pd5_vard_lines$gen[i]], lwd=2)
}

#mtext(latex2exp:::TeX('$\\frac{V_A}{R}$'), 1, line=3)
mtext(latex2exp:::TeX('$V_A / R$'), 1, line=0.5, cex=mtext_cex)
mtext(latex2exp:::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 2, line=3, 
      cex=mtext_cex)

legend(10^{-2.8}, 0.0078, KEEP_GENS, fill=cols[KEEP_GENS-5], title='generation (s)', bty='n', text.col=title_col, border=0, cex=0.6)


# keep <- seq(1, n, length.out=6)
# legend_image <- as.raster(matrix(viridis(n), ncol=1))
# rasterImage(legend_image, 4e-3, 0.002, 4e-3 + 0.002, 0.007)
# text(x=4e-3 + 0.0025, y = seq(0.002+0.0002, 0.007-0.0002, l=5), 
#      adj=c(0, 1/2),  cex=0.9,
#      labels = floor(seq(24,  1, l=5)))
#text(0.6e-2, 0.0075, 'generation', adj=c(0.5, 0.5))


# legend(1e-3 + 0.0001, 0.0075, unique(pd5_dfs$gen)[keep], 
#        fill=viridis(n)[keep], ncol=2, bty='n', 
#        title=latex2exp:::TeX('generation $s$'), border=0, cex=0.8)

dev.off()
par(opar)


## Supplementary figure: single locus stuff
data(sl_pd5f)
data(sl_predf)
pdf('sim-pred-covs-sl.pdf', width=textwidth, height=textwidth/phi)
panel_plot(sl_pd5f, sl_predf, gen, cov, alpha, rec_frac, N, 
           row_lab='$\\alpha = $', col_lab = '$r = $', 
           xlab='generation', 
           y_axis_line=3,
           x_axis_line=3,
           pch_cex=0.8,
           point_alpha=0.6,
           axis_lwd=1.2,
           grp_cols='gray28',#wesanderson::wes_palette('Darjeeling1')[4],
           ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'))
dev.off()



data(vark_types)
pdf('expfit-vark-types.pdf', width=textwidth, height=textwidth/phi)

single_panel_plot(vark_types, row=Va, col=genlen, x=gen, 
                  y=variance, groups=var_type, type=c('p', 'l', 'l', 'l', 'l'), 
                  row_lab="$V_A = $",
                  col_lab='$R = $',
                  axis_lwd=1.2,
                  y_axis_line=3,
                  x_axis_line=3,
                  row_lab_cex=1.1,
                  col_lab_cex=1.1,
                  N=1e3,
                  grp_cols=wesanderson::wes_palette('Darjeeling1')[c(4, 1, 2, 3, 5)],
                  ylab='variance', xlab='generation',
                  legend_prop=1.2, legend_title='variance type')

dev.off()



## Figure 2 -- alt with SSH
data(predf_ssh)
data(pd5f_ssh)

pdf('sim-pred-covs-L500-ssh.pdf', width=textwidth, height=textwidth/phi)
panel_plot(pd5f_ssh, predf_ssh, gen, cov, Va, genlen, var_type, 
           row_lab='$V_A = $', col_lab = '$r = $', 
           xlab='generation', 
           y_axis_line=3,
           x_axis_line=3,
           pch_cex=0.8,
           point_alpha=0.6,
           legend_prop=1.4,
           axis_lwd=1.2,
           ylab=latex2exp::TeX('$cov(\\Delta p_5, \\; \\Delta p_s)$'), 
           legend_title='variance type', data_no_group=TRUE)
dev.off()


## Supp Figure -- single locus var(w) trajectory plot

data(sl_wvar_trajs)
pdf('wvar-sl-trajectories.pdf', width=textwidth/2, height=textwidth/2)

p0_bins <- unique(sl_wvar_trajs$p0_bin)

par(mar=c(4.1, 4.1, 4.1, 2.1))
plot.new()
plot.window(xlim=c(0, 100), ylim=c(0, 0.0052))

cols <- wesanderson::wes_palette('Darjeeling1')

for (i in seq_along(p0_bins)) {
  p0b  <- p0_bins[i]
  p0d <- filter(sl_wvar_trajs, p0_bin==p0b)
  nboots <- unique(p0d$.id)
  for (boot in nboots) {
    p0d_boot <- filter(p0d, .id==boot)
    #lines(lowess(p0d_boot$gen, p0d_boot$wvar), col=cols[i]) 
    lines(predict(smooth.spline(p0d_boot$gen, p0d_boot$wvar), 
                  data=data.frame(gen=0:100)), 
          col=adjustcolor(cols[i], alpha.f=0.4), lwd=1.6)
  }
}
cex <- 1

axis(1, las=1, 
     at=seq(0, 100, length.out=5),
     line=0,
     tck=0.01,
     lwd=1.2,
     padj=-2.2,
     cex.axis=0.7*cex, 
     col.axis=axs_col)

axis(2, at=c(0, 0.001, 0.002, 0.003, 0.004, 0.005), col.axis=axs_col, 
     cex.axis=cex*0.7, col=axs_col, 
     line=0,
     lwd=1.2, tck=0.018, las=1, hadj=0.6)

mtext('generation', 1, line=1.8, col='gray10')
mtext('genetic variance', 2, line=3, col='gray10')

legend(66, 0.0049, legend=p0_bins, fill=cols[1:length(nboots)],
       title='selected allele\ninitial frequency',
       # adj=c(0, 0.45),
       cex=0.6,
       bty='n', text.col='gray10',
       border=0)
dev.off()

## Figure -- MoM fits 
data(mom_fitsd)
wescols <- wesanderson::wes_palette('Darjeeling1')[c(4, 1, 2, 3, 5)]
pdf('mom-fits-va1.pdf', width=textwidth/2, height=textwidth/2)
#mom_fitsd$genlen <- droplevels(mom_fitsd$genlen)

with(mom_fitsd,
     plot(emp_va, va0_est, col=wescols[as.factor(genlen)],
          xlab='', ylab='',
          pch=19, cex=0.5, axes=FALSE,
          ylim=c(1e-4, 10),
          xlim=c(1e-4, 1),
          lwd=0, log='xy'))
abline(a=0, b=1, col=dashed_gray, lty=2)
axes_col <- 'gray12'
axis_cex <- 1*0.7
axis(1, at=10^(-(4:0)), labels=pretty_log(-(4:0)), 
     padj=-1.5,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(2, at=10^(c(-(4:0), 1)), labels=pretty_log(c(-(4:0), 1)), 
     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)
mtext(latex2exp:::TeX('empirical $V_A(1)$'), 1, line=1.3, cex=mtext_cex)
mtext(latex2exp:::TeX('estimated $V_A(1)$'), 2, line=1.5, cex=mtext_cex)
title_col='grey12'
legend(1e-4, 5, inset=0, legend=unique(mom_fitsd$genlen), fill=wescols,
           title=latex2exp:::TeX('level of\nrecombination (Morgans)'),
           cex=0.6,
           ncol=2,
           bty='n', text.col=title_col,
           border=0)
#text(1e-4, 11, 'A', line=2, cex=2, xpd=TRUE)
dev.off()



#pdf('mom-fits-N-Va.pdf', width=textwidth/2, height=textwidth/2)
#wescols <- wesanderson::wes_palette('Zissou1')
#wescols <- brewer.pal(5,"Dark2")
#mom_fitsd_subset <- mom_fitsd %>% 
#  filter(Va %in% c(0.001, 0.005, 0.01, 0.05, 0.1)) #%>%
#  #filter(genlen %in% c(0.1, 0.5, 1.5))
#boxplot(Nest ~ Va + genlen, mom_fitsd_subset,
#        axes=FALSE, ylim=c(0, 1.5e3), outline=FALSE,
#        #col=wescols,
#        border=wescols,
#        pars = list(boxcol = "transparent", medlty = "blank", 
#                    medpch=16, whisklty = c(1, 1),
#                    medcex = 0.5, outcex = 0, staplelty = "blank"))
#abline(h=1e3, col=dashed_gray, lty=2)
#x <- 1:20
#y <- mom_fitsd_subset %>% group_by(Va, genlen) %>% 
#  summarize(Ne_est=mean(Ne_est, na.rm=TRUE)) %>%
#  arrange(genlen, Va) %>% pull(Ne_est)
#points(x, y, pch='-', col=wescols)
#axis(1, at=c(1, 5), labels=c("", "0.01"),
#     hadj=3,
#     padj=-2.4,
#     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#axis(1, at=c(6, 10), 
#     labels=c("", "0.1"), hadj=4,
#     padj=-2.4,
#     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)

#axis(1, at=c(11, 15), 
#     labels=c("", "0.5"), hadj=4,
#     padj=-2.4,
#     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#axis(1, at=c(16, 20), 
#     labels=c("", "1.5"), hadj=4,
#     padj=-2.4,
#     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#axis(2, at=c(0, 500, 1e3, 1.5e3), 
#     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#mtext(latex2exp:::TeX('level of recombination (Morgans)'), 
#      1, line=1.2, cex=mtext_cex)
#mtext(latex2exp:::TeX('estimated N'), 2, line=2, cex=mtext_cex)
#dev.off()


# The original standalone estimated N figure
pdf('mom-fits-N.pdf', width=textwidth/2, height=textwidth/2)
boxplot(Nest ~ genlen, mom_fitsd, axes=FALSE, ylim=c(0, 2e3), outline=FALSE,
        pars = list(boxcol = "transparent", medlty = "blank", 
                    medpch=16, whisklty = c(1, 1),
                    medcex = 0.7,  outcex = 0, staplelty = "blank"))
abline(h=1e3, col=dashed_gray, lty=2)
axis(1, at=unique(mom_fitsd$genlen), labels=levels(mom_fitsd$genlen), 
     padj=-2.4,
     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
axis(2, at=c(0, 500, 1e3, 1.5e3, 2e3), 
     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)
mtext(latex2exp:::TeX('level of recombination (Morgans)'), 
      1, line=1.2, cex=mtext_cex)
mtext(latex2exp:::TeX('estimated N'), 2, line=2, cex=mtext_cex)
#text(0, 2010, 'B', line=2, cex=2, xpd=TRUE)
dev.off()




## Estimation Figure -- both Va and N estimation, with N boxplot grouped by Va
#pdf('mom-fits-both.pdf', width=textwidth, height=textwidth/2 * 0.9)

#layout(matrix(c(1, 2, 3, 4), ncol=2, byrow=TRUE), heights=c(0.1, 0.90))
#opar <- par(no.readonly=TRUE)
#par(mar=c(0,0,0,0), xpd=TRUE)
#plot.new()
#text(0.1, 0.3, "A", cex=1.5)

#plot.new()
#text(0.1, 0.3, "B", cex=1.5)

#data(mom_fitsd)
#wescols <- wesanderson::wes_palette('Darjeeling1')[c(4, 1, 2, 3, 5)]
#mom_fitsd$genlen <- droplevels(mom_fitsd$genlen)

#par(mar=c(4.1, 3.1, 0.1, 2.1))
#with(mom_fitsd[mom_fitsd$va0_est > 1e-4, ],
#     plot(emp_va, va0_est, col=wescols[as.factor(genlen)],
#          xlab='', ylab='',
#          pch=19, cex=0.5, axes=FALSE,
#          ylim=c(1e-4, 10),
#          xlim=c(1e-4, 1),
#          lwd=0, log='xy'))
##abline(a=0, b=1, col=dashed_gray, lty=2)
#segments(1e-4, 1e-4, 1, 1, col=dashed_gray, lty=2)
#axes_col <- 'gray12'
#axis_cex <- 1*0.7
#axis(1, at=10^(-(4:0)), labels=pretty_log(-(4:0)), 
#     padj=-1.5,
#     tck=0.01, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#axis(2, at=10^(c(-(4:0), 1)), labels=pretty_log(c(-(4:0), 1)), 
#     tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex)
#mtext(latex2exp:::TeX('empirical $V_A(1)$'), 1, line=1.3, cex=mtext_cex)
#mtext(latex2exp:::TeX('estimated $V_A(1)$'), 2, line=1.5, cex=mtext_cex)
#legend(1e-4, 5, inset=0, legend=unique(mom_fitsd$genlen), fill=wescols,
#           title=latex2exp:::TeX('level of\nrecombination (Morgans)'),
#           cex=0.6,
#           ncol=2,
#           bty='n', text.col=title_col,
#           border=0)
##text(1e-4, 11, 'A', line=2, cex=2, xpd=TRUE)

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
##text(0, 2010, 'B', line=2, cex=2, xpd=TRUE)
#dev.off()


## Figure Prop var due to linked sel figure -- conservative G
data(mom_fits_g)
pdf('estimate-g.pdf', width=textwidth, 
    height=textwidth * 0.9)

par(mar=c(4, 4, 0, 1), oma=c(0, 0, 0, 0))
wesboxcols <- wesanderson::wes_palette("Zissou1", length(unique(mom_fits$genlen)), type = "continuous")
box_va_params <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1)
yscale <- seq(-1, 1, by=1/4)
x <- boxplot(G_after ~ genlen + Va, mom_fits_g %>% filter(Va %in% box_va_params), axes=FALSE, 
        col = wesboxcols,
        outline=FALSE,
        #ylim=c(-3/4,1.1),
        ylim=c(min(yscale) - 0.1, max(yscale) + 0.1),
        at=atsfun(7, length(box_va_params)),
        ylab='', xlab='',
        pars = list(boxcol = "transparent", medlty = "blank", 
                    medpch=19, whisklty = c(1, 1),
                    whiskcol='gray61',
                    medcex = 0.5, outcex = 0, staplelty = "blank"))

abline(h=0, lty='dashed', col=adjustcolor(dashed_gray, alpha.f=0.4))
abline(h=1, lty='dashed', col=adjustcolor(dashed_gray, alpha.f=0.4))
ats <- atsfun(7, length(box_va_params), 2, ends=TRUE)
legend(34, -1/4, legend=unique(mom_fits$genlen), fill = wesboxcols,
       cex=legend_cex,
       bty='n', text.col=title_col,
       title='level of recombination\n(Morgans)',
       ncol=2,
       inset=0,
       border=0)
axis(2, at=yscale, tck=0.018, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex) 
buf <- 1
for (i in 1:nrow(ats)) {
  axis(1, at=c(ats[i,1]-buf, ats[i,2] + buf), labels=c('', ''), line=0.4, padj=-0.9, tck=0.01, 
       col=axes_col)
  axis(1, at=mean(ats[i,]), labels=box_va_params[i], line=0.5, padj=-0.9, tck=0, lwd=0,
       col=axes_col)
}
mtext_cex <- 1.1
mtext(latex2exp:::TeX('level of additive genetic variation'), 1, line=2.8, cex=mtext_cex)
mtext(latex2exp:::TeX("fraction of variance due to linked selection, $G$"), 2, line=2, cex=mtext_cex)
dev.off()


## Figure Prop var due to linked sel figure
data(mom_fits_gp)
pdf('estimate-gp.pdf', width=textwidth, 
    height=textwidth * 0.9)

par(mar=c(4, 4, 0, 1), oma=c(0, 0, 0, 0))
wesboxcols <- wesanderson::wes_palette("Zissou1", length(unique(mom_fits$genlen)), type = "continuous")
box_va_params <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1)
yscale <- seq(-1, 1, by=1/4)
x <- boxplot(Gp_after ~ genlen + Va, mom_fits_gp %>% filter(Va %in% box_va_params), axes=FALSE, 
        col = wesboxcols,
        outline=FALSE,
        #ylim=c(-3/4,1.1),
        ylim=c(min(yscale) - 0.1, max(yscale) + 0.1),
        at=atsfun(7, length(box_va_params)),
        pars = list(boxcol = "transparent", medlty = "blank", 
                    medpch=19, whisklty = c(1, 1),
                    whiskcol='gray61',
                    medcex = 0.5, outcex = 0, staplelty = "blank"))

abline(h=0, lty='dashed', col=adjustcolor(dashed_gray, alpha.f=0.4))
abline(h=1, lty='dashed', col=adjustcolor(dashed_gray, alpha.f=0.4))
legend_cex <- 0.9
ats <- atsfun(7, length(box_va_params), 2, ends=TRUE)
legend(34, -0.35, legend=unique(mom_fits$genlen), fill = wesboxcols,
       cex=legend_cex,
       bty='n', text.col=title_col,
       title='level of recombination\n(Morgans)',
       ncol=2,
       inset=0,
       border=0)
axis(2, at=yscale, tck=0.005, hadj=0.4, las=1, col=axes_col, lwd=1.2, cex.axis=axis_cex) 
buf <- 1
for (i in 1:nrow(ats)) {
  axis(1, at=c(ats[i,1]-buf, ats[i,2] + buf), labels=c('', ''), line=0.4, padj=-0.9, tck=0.007, 
       col=axes_col)
  axis(1, at=mean(ats[i,]), labels=box_va_params[i], line=0.5, padj=-0.9, tck=0, lwd=0,
       col=axes_col)
}
mtext_cex <- 1.1
mtext(latex2exp:::TeX('level of additive genetic variation'), 1, line=2.8, cex=mtext_cex)
mtext(latex2exp:::TeX("fraction of variance due to linked selection, $G'$"), 2, line=2, cex=mtext_cex)

dev.off()



## Supplementary Figure -- neutral Ne estimates
data(Ne_ests)
data(Ne_fits)
pdf('supp-Ne-est-neutral.pdf', width=textwidth, height=textwidth/phi)

panel_plot(Ne_ests, Ne_fits, t0, Ne, N, rho, group,
                  col_lab="$\\rho = $",
                  row_lab='$N = $',
                  pch_cex=1.2,
                  axis_lwd=1.2,
                  y_axis_line=3,
                  x_axis_line=3,
                  row_lab_cex=1.1,
                  col_lab_cex=1.1,
                  grp_cols=dashed_gray, lty=2,
                  ylab=latex2exp::TeX('$N_e$ estimate'), xlab='generation',
                  legend_prop=0.0001, legend_title='')

dev.off()

## Supplementary Figure -- neutral LD estimates
data(neut_het_fits)
data(neut_het_sims)
pdf('supp-het-neutral.pdf', width=textwidth, height=textwidth/phi)

panel_plot(neut_het_sims, neut_het_fits, x=gen, y=het, col=rho, row=N, groups=group,
                  row_lab="$\\rho = $",
                  col_lab='$N = $',
                  pch_cex=0.5,
                  axis_lwd=1.2,
                  y_axis_line=3,
                  x_axis_line=3,
                  row_lab_cex=1.1,
                  col_lab_cex=1.1,
                  data_col='gray40',
                  lwd=3,
                  data_no_group=TRUE,
                  point_alpha=1,
                  line_alpha=0.6,
                  grp_cols=wesanderson::wes_palette('Darjeeling1')[1],,
                  ylab=latex2exp::TeX('Heterozygosity'), xlab='generation',
                  legend_prop=0.0001, legend_title='')

dev.off()

#

## Supplementary Figure -- neutral LD estimates
data(neut_ld_sims)
data(neut_ld_fits)
pdf('supp-ld-neutral.pdf', width=textwidth, height=textwidth/phi)

neut_ld_fits  <-  neut_ld_fits %>% mutate(grp = as.factor('L'))

panel_plot(neut_ld_sims, neut_ld_fits, x=gen, y=D, col=N, row=r, groups=grp,
                  row_lab="$r = $",
                  col_lab='$N = $',
                  pch_cex=0.5,
                  axis_lwd=1.2,
                  y_axis_line=3,
                  x_axis_line=3,
                  row_lab_cex=1.1,
                  col_lab_cex=1.1,
                  data_col='gray40',
                  lwd=3,
                  data_no_group=TRUE,
                  point_alpha=1,
                  line_alpha=0.6,
                  grp_cols=wesanderson::wes_palette('Darjeeling1')[1],,
                  ylab=latex2exp::TeX('Linkage Disequilibrium (D)'), xlab='generation',
                  legend_prop=0.0001, legend_title='')

dev.off()

# the code for the fluctuating selection figure 
source('./fluct-plot.r')

# code for new method of moment plots, for finite sample,
# and relative error
source('./both-plot-alt-finite.r')
source('./both-plot-alt.r')
source('./relative-error.r')
