# fluct-plot.r -- a hellish amount of base plot code for a complex figure
# ¯\ _ (ツ) _ /¯ 

# TODO look at flip r sim file, this should be what creates cov mats

pkgload::load_all()
textwidth <- 7 # 426 pts to inches
phi <- (1 + sqrt(5))/2
axs_col <- 'gray12'

pdf('fluct-sel.pdf', width=textwidth, height=textwidth/phi * 0.8)

data(pd5_fluctf)
data(pred_fluctf) ## for predictions
data(cumd)
data(cumd_abs)
pd5_fluctf2 <- pd5_fluctf %>% filter(genlen == 0.1, Va == 0.05, L==500, N==1e3)
pred_fluctf2 <- pred_fluctf %>% filter(genlen == 0.1, Va == 0.05, L==500, N==1e3)


#ggplot(pd5_fluctf2, aes(gen, cov)) +  geom_point()

## plot
#cols <- wesanderson::wes_palette('Darjeeling1')[c(3, 5)]
cols <- wesanderson::wes_palette('Zissou1', 16, 'continuous')
cols_before <- cols[c(1, 5)]
cols_after <- cols[c(13, 9)]
subfig_col <- 'gray12'
#cols2 <- adjustcolor(wesanderson::wes_palette('Zissou1'), alpha.f=0.4)
cols2 <- adjustcolor(wesanderson::wes_palette('Darjeeling1'), alpha.f=0.3)
#red_col <- wesanderson::wes_palette('Darjeeling1')[1]


ymin <- -6e-4
ymax <- 8e-4
layout_mat <- matrix(c(1, 2, 3, 
                       1, 2, 4), ncol=3, byrow=TRUE)
layout(layout_mat, widths=c(0.64, 0.07, 0.47))
par(mar=c(1, 2.9, 0, 0), oma=c(2.5, 2.5, 1.1, 1.5))
plot(cov ~ gen, pd5_fluctf2, axes=FALSE, pch=19, cex=0.8,
     xlim=c(0, 51), ylim=c(ymin, ymax), type='n')
#segments(15, -3e-4, 15, 3e-4, 
#         col=adjustcolor(red_col, alpha.f=0.5), lty=5)
rect(-9, ymin, 0, ymax, col='gray96', border=0)
rect(0, -8e-4, 5, ymax, col='gray96', border=0)
rect(5, -8e-4, 14.5, ymax, col=adjustcolor(cols[1], alpha.f=0.4), border=0)
rect(14.5, -8e-4, 50, ymax, col=adjustcolor(cols_after[2], alpha.f=0.4), border=0)
points(cov ~ gen, pd5_fluctf2 %>% filter(gen<50), pch=19, cex=1.,
       col='gray30')

lines(pred ~ gen, pred_fluctf2 %>% filter(gen<50, gen > 5, gen < 15), pch=19, lwd=3, cex=1.,
       col=adjustcolor('gray40', alpha.f=0.7))


lines(pred ~ gen, pred_fluctf2 %>% filter(gen<50, gen >= 5, gen >= 15), pch=19, lwd=3, cex=1.,
       col=adjustcolor('gray40', alpha.f=0.7))


text_cex <- 1
text_base <- 7e-4
text(1.5, text_base, latex2exp:::TeX('$w(z_i)=1$'), cex=text_cex, col=axs_col)
text(0 + (15+5)/2, text_base, latex2exp:::TeX('$w(z_i)=e^{z_i}$'), cex=text_cex, col=axs_col)
text(15 + (20+15)/2, text_base, latex2exp:::TeX('$w(z_i)=e^{- z_i}$'), cex=text_cex, col=axs_col)

##text(1.5, -2.7e-4, latex2exp:::TeX('$w(z_i) = 1$'), cex=0.7, col=axs_col)
#text(1.5, -2.7e-4, latex2exp:::TeX('no\nsel.'), cex=text_cex, col=axs_col)
#text(0 + (15+5)/2, -2.7e-4, latex2exp:::TeX('$w(z_i) = e^{z_i}$'), cex=text_cex, col=axs_col)
#text(15 + (20+15)/2, -2.7e-4, latex2exp:::TeX('$w(z_i) = e^{- z_i}$'), cex=text_cex, col=axs_col)

axis(1, col.axis=axs_col, col=axs_col, 
     line=0,
     lwd=1.2,
     padj=-1.2,
     cex.axis=0.9,
     tck=0.01)
mtext(latex2exp::TeX("cov$(\\Delta p_5, \\; \\Delta p_s)$"), 2, line=3.5)

segments(0, 0, 50, 0, col=dashed_gray, lty=2)
#text(22.7, 2.5e-4, 'change in direction\nof selection', offset=1,
#     col=red_col)
axis(2, col.axis=axs_col, col=axs_col, las=1,
     line=0,
     hadj=0.75,
     lwd=1.2,
     cex.axis=0.9,
     tck=0.01)
mtext('generation', 1, line=1.6)

mtext('A', 3, cex=1.2, adj=-.2, line=-1.4, col=subfig_col)

plot.new()


# there's a subtle off by one issue here; temp_cov() uses Δp_{t-1}.
# 4th row/col is 5th timepoint
end <- 21
dd <- cumd %>% filter(genlen == 0.1, Va == 0.02) 
tmp = dd %>% spread(type, val) %>% mutate(tot=cov+var)
dd_mat <- tmp %>% select(gen, cov, var, tot) %>% t()
mat <- dd_mat[c(3, 2), 1:end] 
# divide through by starting time
mat <- sweep(mat, 2, dd_mat[1,1:end]-3, '/')

# this is a bunch of messy code that gets barplot() to move the negative bars down,
# not plotting them above. It does this if we have rows ordered c(2, 3) (see above),
# but we want the variance to be on the bottom and it doesn't handle it
# oh base graphics
mat_org <- mat
mat[is.na(mat)] <- 0 
mat[1, ] <- ifelse(mat[2,] <= 0, mat[1,] + mat[2,], mat[1,])
mat[2, ] <- ifelse(mat[2,] < 0, 0, mat[2,])
flip <- 10
mat2 <- cbind(rbind(mat[, 1:flip], matrix(0, ncol=flip, nrow=2)),
              rbind(matrix(0, ncol=ncol(mat)-flip, nrow=2), mat[, (flip+1):ncol(mat)]))
x <- mat2 %>% barplot(col=c(cols_before, cols_after), border=0, axes=FALSE)
mat <- mat_org
overlay <- mat[2, ]
overlay[overlay > 0] <- 0
#overlay %>% barplot(col=adjustcolor(cols[1], alpha.f=1), border=0, add=TRUE, axes=FALSE, axisnames=FALSE)
#overlay %>% barplot(col=adjustcolor(cols[2], alpha.f=0.4), border=0, add=TRUE, axes=FALSE, axisnames=FALSE)
overlay %>% barplot(col=adjustcolor(cols_after[2], alpha.f=1), border=0, add=TRUE, axes=FALSE, axisnames=FALSE)
axis(1, at=x[seq(1, 25, by=5)], labels=rep('', 5), col.axis=axs_col, col=axs_col, line=0.5, 
     lwd=1.2, padj=-0.9, tck=0.01, las=1)
segments(0, 1/2e3, 25, 1/2e3, col='gray12', lty=2)
#mtext(latex2exp::TeX('$1/2N$'), 4, xpd=TRUE, cex=0.5, line=-0, adj=0.25, padj=-1.5, srt=45)
mtext(latex2exp::TeX('$1/2N$'), 1, xpd=TRUE, cex=0.5, line=-2.6, adj=1.05)
axis(2, col.axis=axs_col, col=axs_col, las=1, line=0, hadj=0.7, lwd=1.2, tck=0.01, cex.axis=0.9)
#mtext('cumulative variance/covariance', 2, line=3.1, adj=1.4)
mtext('cumulative\ncov + var', 2, line=3.1, cex=0.9)
# mtext('cumulative var/cov', 2, line=3.1, adj=-23)
ly <- 0.0023
lx <- 16.5
legend(lx, ly, legend=c('', ''), fill=c(cols_before),
       title='',
       adj=c(0.1, 0.4),
       cex=1,
       bty='n', text.col='gray10',
       border=0)
legend(lx+1.8, ly, legend=c('var', 'cov'), fill=c(cols_after),
       title='',
       text.width=c(0, 0),
       adj=c(0.1, 0.4),
       cex=1,
       bty='n', text.col='gray10',
       border=0)
text(20.9, ly-0.0004, 'before flip', srt=45, adj=c(0, 0))
text(23, ly-0.0004, 'after flip', srt=45, adj=c(0, 0))
mtext('B', 3, cex=1.2, adj=.029, line=-1.5, col=subfig_col)


dd_abs <- cumd_abs %>% filter(genlen == 0.1, Va == 0.02)  
tmp_abs = dd_abs %>% spread(type, val) %>% mutate(tot=cov+var)
dd_abs_mat <- tmp_abs %>% select(gen, cov, var, tot) %>% t()
mat_abs <- dd_abs_mat[c(3, 2), 1:end] 
mat_abs <- sweep(mat_abs, 2, dd_abs_mat[1,1:end]-3, '/')

mat_abs2 <- cbind(rbind(mat_abs[, 1:flip], matrix(0, ncol=flip, nrow=2)),
              rbind(matrix(0, ncol=ncol(mat_abs)-flip, nrow=2), mat_abs[, (flip+1):ncol(mat_abs)]))

mat_abs2 %>% barplot(col=c(cols_before, cols_after), border=0, axes=FALSE)
axis(1, at=x[seq(1, 25, by=5)], labels=seq(5, 25, by=5), col.axis=axs_col, col=axs_col, line=0.2, 
     lwd=1.2, padj=-1.4, tck=0.01, las=1)
# abline(h=1/2e3, col='gray24', lty=2)
segments(-1, 1/2e3, 25, 1/2e3, col='gray12', lty=2)
mtext(latex2exp::TeX('$1/2N$'), 1, xpd=TRUE, cex=0.5, line=-1.7, adj=1.05)
axis(2, col.axis=axs_col, col=axs_col, las=1, 
     line=0, hadj=0.63, lwd=1.2, tck=0.01, cex.axis=0.9)
mtext('cumulative\n|cov| + var', 2, line=3.1, cex=0.9)
mtext('generation', 1, line=1.6)
mtext('C', 3, cex=1.2, adj=.029, line=-1.5, col=subfig_col)

dev.off()
