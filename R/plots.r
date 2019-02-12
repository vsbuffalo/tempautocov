## plots.r
# library(lattice)
# library(latticeExtra)

dashed_gray <- 'gray58'

# create the ats vector for boxplot(), with spacing
atsfun <- function(neach, ngroups, space=2, ends=FALSE) {
  ats <- rep(1, (neach * ngroups))
  ats[neach*(1:ngroups) + 1L] <- ats[neach*(1:ngroups) + 1L] + space
  ats <- cumsum(ats)
  if (!ends) {
    return(ats[!is.na(ats)])
  }
  a <- atsfun(neach, ngroups, ends=FALSE)[seq(1, neach*ngroups, neach)]
  b <- atsfun(neach, ngroups, ends=FALSE)[seq(neach, neach*ngroups, neach)]
  cbind(a, b)
}



finite_min <- function(x) {
  min(x[is.finite(x)], na.rm=TRUE)
}

finite_max <- function(x) {
  max(x[is.finite(x)], na.rm=TRUE)
}

pretty_rng <- function(xmin, xmax, nticks) 
  range(labeling::extended(xmin, xmax, nticks))

nunique <- function(x) length(unique(x))

panel_plot <- function(data, fits, x, y, row, col, groups, 
                       ynticks=5, xnticks=4, 
                       point_alpha=0.5,
                       line_alpha=1,
                       fontsize=12,
                       row_lab='row = ',
                       col_lab='col = ',
                       legend_title='legend',
                       xlab='x',
                       ylab='y',
                       axis_label_prop=0.05,
                       axis_cex=1,
                       mtext_cex=1,
                       legend_cex=1.1,
                       row_lab_cex=1.4,
                       col_lab_cex=1.4,
                       legend_prop=0.4,
                       y_axis_line=0.5,
                       x_axis_line=-0.9,
                       pretty_rng=TRUE,
                       lwd=2,
                       axis_lwd=1,
                       pch_cex=0.2,
                       oma=c(5, 5, 2, 1),
                       grp_cols=wesanderson::wes_palette('Darjeeling1'),
                       data_col='gray3',
                       title_col='grey12',
                       subtitle_col='grey34',
                       lty=1,
                       free_y=TRUE,
                       data_no_group=FALSE) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  axs_col <- 'gray38'
  has_group <- !missing(groups)

  # quote variables
  row <- rlang::enquo(row)
  col <- rlang::enquo(col)

  # get row/col names
  row_var <- rlang::quo_text(row)
  col_var <- rlang::quo_text(col)
  # get var names
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x_var <- rlang::quo_text(x)
  y_var <- rlang::quo_text(y)

  if (has_group) {
    group <- rlang::enquo(groups)
    grp_var <- rlang::quo_text(group)
    grp_levels <- levels(fits[[grp_var]])
  }
 
  # now, nest data into each panel 
  data_nested <- data %>% group_by(UQ(row), UQ(col)) %>% nest(.key='data')
  fits_nested <- fits %>% group_by(UQ(row), UQ(col)) %>% nest(.key='fits')

  all_data <- left_join(data_nested, fits_nested) %>% 
              mutate(data_ymin=map_dbl(data, ~ finite_min(.[[y_var]])),
                     data_ymax=map_dbl(data, ~ finite_max(.[[y_var]])),

                     data_xmin=map_dbl(data, ~ finite_min(.[[x_var]])),
                     data_xmax=map_dbl(data, ~ finite_max(.[[x_var]])),

                     fits_ymin=map_dbl(fits, ~ finite_min(.[[y_var]])),
                     fits_ymax=map_dbl(fits, ~ finite_max(.[[y_var]])),
 
                     fits_xmin=map_dbl(fits, ~ finite_min(.[[x_var]])),
                     fits_xmax=map_dbl(fits, ~ finite_max(.[[x_var]]))) 
  # add row/cols as factors
  all_data <- all_data %>% mutate(row_fct=as.factor(UQ(row)), 
                                  col_fct=as.factor(UQ(col)))
 
 
  # if free_y, the rows have same y range as determined by largest
  # in data, fits.
  if (free_y) {
    all_data <- all_data %>% group_by(UQ(row)) 
  }
  all_data <- all_data %>%
                  mutate(ymax=pmax(max(data_ymax), max(fits_ymax)), 
                         ymin=pmin(min(data_ymin), min(fits_ymin)),
                         xmax=pmax(max(data_xmax), max(fits_xmax)), 
                         xmin=pmin(min(data_xmin), min(fits_xmin)))

  all_data <- all_data %>% 
   mutate(pretty_ymin=map2_dbl(ymin, ymax, 
                      ~ pretty_rng(.x, .y, nticks=ynticks)[1]), 
          pretty_ymax=map2_dbl(ymin, ymax,
                      ~ pretty_rng(.x, .y, nticks=ynticks)[2]))
  all_data$nplot <- 1:nrow(all_data)
  dims <- c(nunique(data[[row_var]]), nunique(data[[col_var]]))

  # figure out which panels have axes labels
  all_data <- all_data %>% 
    mutate(y_axis=as.integer(col_fct) == 1, x_axis=as.integer(row_fct)==dims[1])

  # now, create all plots
  npans <- prod(dims)
  layout_vec <- 1:npans
  layout_mat <- matrix(layout_vec, nrow=dims[1], ncol=dims[2], byrow=TRUE)
  layout_widths <- rep.int(1, ncol(layout_mat))
  layout_heights <- rep.int(1, nrow(layout_mat))
  if (!is.na(legend_title)) {
    layout_mat <- cbind(layout_mat, npans+1L)
  } else {
    legend_prop <- 0
  }
  # add plots for y/x labels
  layout_mat <- rbind(layout_mat, npans+2L)
  layout_mat <- cbind(npans+3L, layout_mat)
  layout_widths <- c(axis_label_prop, layout_widths, legend_prop)
  layout_heights <- c(layout_heights, axis_label_prop)

  # create dummy panels 
  bottom_left <- row(layout_mat) == nrow(layout_mat) & col(layout_mat) == 1
  layout_mat[bottom_left] <- npans+4L
  bottom_right <- (row(layout_mat) == nrow(layout_mat) & 
                   col(layout_mat) == ncol(layout_mat))
  layout_mat[bottom_right] <- npans+5L

  layout(layout_mat, widths=layout_widths, heights=layout_heights)
  par(mar=c(0, 1, 1, 0.5),
      #oma=c(5*axis_cex, 5*axis_cex, 2, 1*legend_cex),
      oma=oma,
      ps=fontsize)
  for (i in all_data$nplot) {
    plot.new()

    # set up coordinates.
    ymin <- all_data$ymin[i]
    ymax <- all_data$ymax[i]
    xmin <- all_data$xmin[i]
    xmax <- all_data$xmax[i]
    pretty_ymax <- all_data$pretty_ymax[i]
    pretty_ymin <- all_data$pretty_ymin[i]
    fac <- 0.08
    if (pretty_rng) {
      plot.window(xlim=c(xmin*(1-fac), xmax*(1+fac)), 
                  ylim=c(pretty_ymin*(1-fac), pretty_ymax*(1+fac)))
    } else {
      plot.window(xlim=c(xmin*(1-fac), xmax*(1+fac)), 
                  ylim=c(ymin*(1-fac), ymax*(1+fac)))
    }
    panel_data <- all_data$data[[i]]
    panel_fits <- all_data$fits[[i]]
    
    # group data
    data_grp <- panel_data[[grp_var]]
    fit_grp <- panel_fits[[grp_var]]
    segments(xmin, 0, xmax, 0, col=dashed_gray, lty=2)

    # data
    if (!data_no_group) {
    points(panel_data[[x_var]], panel_data[[y_var]], 
           col=adjustcolor(grp_cols[data_grp], alpha.f=point_alpha),
           pch=19, cex=pch_cex, lwd=0)
    } else {
      points(panel_data[[x_var]], panel_data[[y_var]], 
           col=adjustcolor(data_col, alpha.f=point_alpha),
           pch=19, cex=pch_cex, lwd=0)
    }

    for (gi in seq_along(levels(fit_grp))) {
      grp_i <- which(panel_fits[[grp_var]] == levels(fit_grp)[gi])
      lines(panel_fits[[x_var]][grp_i], 
            panel_fits[[y_var]][grp_i],
            col=adjustcolor(grp_cols[gi], alpha.f=line_alpha), lwd=lwd, lty=lty)
    }
    panel_row <- all_data[[row_var]][[i]]
    panel_col <- all_data[[col_var]][[i]]

    corners <- par('usr')
    if (i <= dims[2]) {
      text(x=mean(corners[1:2]), y=corners[4]*(1+0.04),
           latex2exp:::TeX(paste0(col_lab, prettyNum(panel_col))),
           font=2, cex=col_lab_cex, xpd=NA, col=title_col)
    }
    if (i %% dims[2] == 0) {
      corners = par("usr")
      text(x = corners[2], y = mean(corners[3:4]), 
           latex2exp:::TeX(paste0(row_lab, prettyNum(panel_row))),
           srt = 270, font=2, cex=row_lab_cex, xpd=NA, col=title_col)
      # mtext(
           # latex2exp:::TeX(paste0(row_lab, prettyNum(panel_row))),
           # side=4)
    }
    # text(0.55*corners[2], 0.8*corners[4], adj=c(0, 0),
         # latex2exp::TeX(paste0('$V_A = $', panel_row, ", $r = $", panel_col))
         # )
    # title(latex2exp::TeX(paste0('$V_A = $', panel_row, ", $r = $", panel_col)), 
          # line=-0.5)

    # axes with labels
    labeled_x <- all_data$x_axis[[i]]
    labeled_y <- all_data$y_axis[[i]]
    xticks <- labeling::extended(xmin, xmax, xnticks)
    yticks <- labeling::extended(ymin, ymax, ynticks)
    axis(1, at=xticks, labels=if (labeled_x) prettyNum(xticks) else FALSE,
         col.axis=axs_col, col=axs_col, 
         cex.axis=axis_cex,
         lwd=axis_lwd,
         line=0.4,
         padj=-0.9,
         tck=0.018)
    axis(2, at=yticks, labels=if (labeled_y) prettyNum(yticks) else FALSE,
         col.axis=axs_col, col=axs_col, las=1,
         cex.axis=axis_cex,
         lwd=axis_lwd,
         line=0.5,
         hadj=0.9,
         tck=0.018)

  }
  par(mar=c(0, 0, 0, 0))
  plot.new()
  if (legend_title != '') {
    # dirty hack :| 
  if (has_group && length(grp_levels) > 1) {
    legend('center', inset=0, legend=grp_levels, fill=grp_cols,
           title=legend_title,
           # adj=c(0, 0.45),
           cex=legend_cex,
           bty='n', text.col=title_col,
           border=0)
    }
  } 
  plot.new()
  mtext(xlab, 1, line=x_axis_line, col=title_col, cex=mtext_cex)
  plot.new()
  mtext(ylab, 2, line=y_axis_line, col=title_col, cex=mtext_cex)

  # two dummy plots
  plot.new(); plot.new()

}


# hacked version of previous function to handle one case ¯\_(ツ)_/¯ 
panel_plot2 <- function(data, fits, x, y, row, col, groups_data, 
                        groups_fits,
                       ynticks=5, xnticks=4, 
                       point_alpha=0.5,
                       line_alpha=1,
                       fontsize=12,
                       row_lab='row = ',
                       col_lab='col = ',
                       legend_title1='legend',
                       legend_title2='legend',
                       lg1y=0.1, 
                       lg2y=lg1y + 0.5,
                       xlab='x',
                       ylab='y',
                       axis_label_prop=0.05,
                       axis_cex=1,
                       mtext_cex=1,
                       legend_cex=1.1,
                       row_lab_cex=1.4,
                       col_lab_cex=1.4,
                       legend_prop=0.4,
                       y_axis_line=0.5,
                       x_axis_line=-0.9,
                       pretty_rng=TRUE,
                       lwd=2,
                       axis_lwd=1,
                       pch_cex=0.2,
                       oma=c(5, 5, 2, 1),
                       grp_cols=wesanderson::wes_palette('Darjeeling1'),
                       data_col='gray28',
                       title_col='grey12',
                       subtitle_col='grey34',
                       free_y=TRUE,
                       data_no_group=FALSE) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  axs_col <- 'gray38'
  # stopifnot(!missing(group_fits) && !missing(group_data))

  # quote variables
  row <- rlang::enquo(row)
  col <- rlang::enquo(col)

  # get row/col names
  row_var <- rlang::quo_text(row)
  col_var <- rlang::quo_text(col)
  # get var names
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x_var <- rlang::quo_text(x)
  y_var <- rlang::quo_text(y)

  group_fits <- rlang::enquo(groups_fits)
  group_data <- rlang::enquo(groups_data)
  grp_fits_var <- rlang::quo_text(group_fits)
  grp_data_var <- rlang::quo_text(group_data)
  grp_data_levels <- levels(data[[grp_data_var]])
  grp_fits_levels <- levels(fits[[grp_fits_var]])

  # now, nest data into each panel 
  data_nested <- data %>% group_by(UQ(row), UQ(col)) %>% nest(.key='data')
  fits_nested <- fits %>% group_by(UQ(row), UQ(col)) %>% nest(.key='fits')

  all_data <- left_join(data_nested, fits_nested) %>% 
              mutate(data_ymin=map_dbl(data, ~ finite_min(.[[y_var]])),
                     data_ymax=map_dbl(data, ~ finite_max(.[[y_var]])),

                     data_xmin=map_dbl(data, ~ finite_min(.[[x_var]])),
                     data_xmax=map_dbl(data, ~ finite_max(.[[x_var]])),

                     fits_ymin=map_dbl(fits, ~ finite_min(.[[y_var]])),
                     fits_ymax=map_dbl(fits, ~ finite_max(.[[y_var]])),
 
                     fits_xmin=map_dbl(fits, ~ finite_min(.[[x_var]])),
                     fits_xmax=map_dbl(fits, ~ finite_max(.[[x_var]]))) 
  # add row/cols as factors
  all_data <- all_data %>% mutate(row_fct=as.factor(UQ(row)), 
                                  col_fct=as.factor(UQ(col)))
 
 
  # if free_y, the rows have same y range as determined by largest
  # in data, fits.
  if (free_y) {
    all_data <- all_data %>% group_by(UQ(row)) 
  }
  all_data <- all_data %>%
                  mutate(ymax=pmax(max(data_ymax), max(fits_ymax)), 
                         ymin=pmin(min(data_ymin), min(fits_ymin)),
                         xmax=pmax(max(data_xmax), max(fits_xmax)), 
                         xmin=pmin(min(data_xmin), min(fits_xmin)))

  all_data <- all_data %>% 
   mutate(pretty_ymin=map2_dbl(ymin, ymax, 
                      ~ pretty_rng(.x, .y, nticks=ynticks)[1]), 
          pretty_ymax=map2_dbl(ymin, ymax,
                      ~ pretty_rng(.x, .y, nticks=ynticks)[2]))
  all_data$nplot <- 1:nrow(all_data)
  dims <- c(nunique(data[[row_var]]), nunique(data[[col_var]]))

  # figure out which panels have axes labels
  all_data <- all_data %>% 
    mutate(y_axis=as.integer(col_fct) == 1, x_axis=as.integer(row_fct)==dims[1])

  # now, create all plots
  npans <- prod(dims)
  layout_vec <- 1:npans
  layout_mat <- matrix(layout_vec, nrow=dims[1], ncol=dims[2], byrow=TRUE)
  layout_widths <- rep.int(1, ncol(layout_mat))
  layout_heights <- rep.int(1, nrow(layout_mat))
  if (!is.na(legend_title1)) {
    layout_mat <- cbind(layout_mat, npans+1L)
  } else {
    legend_prop <- 0
  }
  # add plots for y/x labels
  layout_mat <- rbind(layout_mat, npans+2L)
  layout_mat <- cbind(npans+3L, layout_mat)
  layout_widths <- c(axis_label_prop, layout_widths, legend_prop)
  layout_heights <- c(layout_heights, axis_label_prop)

  # create dummy panels 
  bottom_left <- row(layout_mat) == nrow(layout_mat) & col(layout_mat) == 1
  layout_mat[bottom_left] <- npans+4L
  bottom_right <- (row(layout_mat) == nrow(layout_mat) & 
                   col(layout_mat) == ncol(layout_mat))
  layout_mat[bottom_right] <- npans+5L

  layout(layout_mat, widths=layout_widths, heights=layout_heights)
  par(mar=c(0, 1, 1, 0.5),
      #oma=c(5*axis_cex, 5*axis_cex, 2, 1*legend_cex),
      oma=oma,
      ps=fontsize)
  for (i in all_data$nplot) {
    plot.new()

    # set up coordinates.
    ymin <- all_data$ymin[i]
    ymax <- all_data$ymax[i]
    xmin <- all_data$xmin[i]
    xmax <- all_data$xmax[i]
    pretty_ymax <- all_data$pretty_ymax[i]
    pretty_ymin <- all_data$pretty_ymin[i]
    fac <- 0.08
    if (pretty_rng) {
      plot.window(xlim=c(xmin*(1-fac), xmax*(1+fac)), 
                  ylim=c(pretty_ymin*(1-fac), pretty_ymax*(1+fac)))
    } else {
      plot.window(xlim=c(xmin*(1-fac), xmax*(1+fac)), 
                  ylim=c(ymin*(1-fac), ymax*(1+fac)))
    }
    panel_data <- all_data$data[[i]]
    panel_fits <- all_data$fits[[i]]
    
    # group data
    data_grp <- panel_data[[grp_data_var]]
    fit_grp <- panel_fits[[grp_fits_var]]
    segments(xmin, 0, xmax, 0, col=dashed_gray, lty=2)

    # data
    # browser()
    points(panel_data[[x_var]], panel_data[[y_var]], 
           col=adjustcolor(grp_cols[data_grp], alpha.f=point_alpha),
           pch=19, cex=pch_cex, lwd=0)

    for (gi in seq_along(levels(fit_grp))) {
      grp_i <- which(panel_fits[[grp_fits_var]] == levels(fit_grp)[gi])
      #browser()
      lines(panel_fits[[x_var]][grp_i], 
            panel_fits[[y_var]][grp_i],
            col=adjustcolor(data_col, alpha.f=line_alpha),
            lty=ifelse(gi==3, 5, gi), lwd=lwd)
    }
    panel_row <- all_data[[row_var]][[i]]
    panel_col <- all_data[[col_var]][[i]]

    corners <- par('usr')
    if (i <= dims[2]) {
      text(x=mean(corners[1:2]), y=corners[4]*(1+0.04),
           latex2exp:::TeX(paste0(col_lab, prettyNum(panel_col))),
           font=2, cex=col_lab_cex, xpd=NA, col=title_col)
    }
    if (i %% dims[2] == 0) {
      corners = par("usr")
      text(x = corners[2], y = mean(corners[3:4]), 
           latex2exp:::TeX(paste0(row_lab, prettyNum(panel_row))),
           srt = 270, font=2, cex=row_lab_cex, xpd=NA, col=title_col)
      # mtext(
           # latex2exp:::TeX(paste0(row_lab, prettyNum(panel_row))),
           # side=4)
    }
    # text(0.55*corners[2], 0.8*corners[4], adj=c(0, 0),
         # latex2exp::TeX(paste0('$V_A = $', panel_row, ", $r = $", panel_col))
         # )
    # title(latex2exp::TeX(paste0('$V_A = $', panel_row, ", $r = $", panel_col)), 
          # line=-0.5)

    # axes with labels
    labeled_x <- all_data$x_axis[[i]]
    labeled_y <- all_data$y_axis[[i]]
    xticks <- labeling::extended(xmin, xmax, xnticks)
    yticks <- labeling::extended(ymin, ymax, ynticks)
    axis(1, at=xticks, labels=if (labeled_x) prettyNum(xticks) else FALSE,
         col.axis=axs_col, col=axs_col, 
         cex.axis=axis_cex,
         lwd=axis_lwd,
         line=0.4,
         padj=-0.9,
         tck=0.018)
    axis(2, at=yticks, labels=if (labeled_y) prettyNum(yticks) else FALSE,
         col.axis=axs_col, col=axs_col, las=1,
         cex.axis=axis_cex,
         lwd=axis_lwd,
         line=0.5,
         hadj=0.9,
         tck=0.018)

  }
  plot.new()
  par(mar=c(0, 0, 0, 0))
  # browser()
  usr <- par('usr')
  legend(0, lg2y, inset=0, legend=grp_data_levels, fill=grp_cols,
         title=legend_title1,
         # adj=c(0, 0.45),
         cex=legend_cex,
         bty='n', text.col=title_col,
         border=0)

  legend(0, lg1y, inset=0, legend=latex2exp::TeX(grp_fits_levels), 
         lty=c(1, 2, 5),
         col=adjustcolor(data_col, alpha.f=line_alpha),
         title=legend_title2,
         # adj=c(0, 0.45),
         cex=legend_cex,
         lwd=1.7,
         bty='n', text.col=title_col,
         border=0)

  plot.new()
  mtext(xlab, 1, line=x_axis_line, col=title_col, cex=mtext_cex)
  plot.new()
  mtext(ylab, 2, line=y_axis_line, col=title_col, cex=mtext_cex)

  # two dummy plots
  plot.new(); plot.new()

}



single_panel_plot <- function(data, fits, x, y, row, col, groups, 
                       ynticks=4, xnticks=4, 
                       point_alpha=0.5,
                       fontsize=12,
                       row_lab='row = ',
                       col_lab='col = ',
                       row_lab_cex=1.4,
                       col_lab_cex=1.4,
                       legend_title='legend',
                       xlab='x',
                       ylab='y',
                       axis_label_prop=0.05,
                       legend_prop=0.4,
                       pretty_rng=TRUE,
                       type='p',
                       y_axis_line=0.5,
                       x_axis_line=-0.9,
                       lwd=1,
                       axis_lwd=1.2,
                       grp_cols=wesanderson::wes_palette('Darjeeling1'),
                       N=NULL,
                       title_col='grey10',
                       subtitle_col='grey34',
                       free_y=TRUE) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  axs_col <- 'gray38'
  has_group <- !missing(groups)
  has_fits <- !missing(fits)

  # quote variables
  row <- rlang::enquo(row)
  col <- rlang::enquo(col)

  # get row/col names
  row_var <- rlang::quo_text(row)
  col_var <- rlang::quo_text(col)
  # get var names
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x_var <- rlang::quo_text(x)
  y_var <- rlang::quo_text(y)

  if (has_group) {
    group <- rlang::enquo(groups)
    grp_var <- rlang::quo_text(group)
    if (!is.factor(data[[grp_var]])) stop('grouping variable must be factor')
    grp_levels <- levels(data[[grp_var]])
  }
 
  # now, nest data into each panel 
  data_nested <- data %>% group_by(UQ(row), UQ(col)) %>% nest(.key='data')
  if (has_fits) {
    fits_nested <- fits %>% group_by(UQ(row), UQ(col)) %>% nest(.key='fits')

    all_data <- left_join(data_nested, fits_nested) %>% 
      mutate(data_ymin=map_dbl(data, ~ finite_min(.[[y_var]])),
             data_ymax=map_dbl(data, ~ finite_max(.[[y_var]])),

                     data_xmin=map_dbl(data, ~ finite_min(.[[x_var]])),
                     data_xmax=map_dbl(data, ~ finite_max(.[[x_var]])),

                     fits_ymin=map_dbl(fits, ~ finite_min(.[[y_var]])),
                     fits_ymax=map_dbl(fits, ~ finite_max(.[[y_var]])),
 
                     fits_xmin=map_dbl(fits, ~ finite_min(.[[x_var]])),
                     fits_xmax=map_dbl(fits, ~ finite_max(.[[x_var]]))) 
  } else {
    all_data <- data_nested %>%
      mutate(data_ymin=map_dbl(data, ~ finite_min(.[[y_var]])),
             data_ymax=map_dbl(data, ~ finite_max(.[[y_var]])),

             data_xmin=map_dbl(data, ~ finite_min(.[[x_var]])),
             data_xmax=map_dbl(data, ~ finite_max(.[[x_var]])))
  }
  # add row/cols as factors
  all_data <- all_data %>% mutate(row_fct=as.factor(UQ(row)), 
                                  col_fct=as.factor(UQ(col)))
 
 
  # if free_y, the rows have same y range as determined by largest
  # in data, fits.
  if (free_y) {
    all_data <- all_data %>% group_by(UQ(row)) 
  }

  if (has_fits) {
    all_data <- all_data %>%
                  mutate(ymax=pmax(max(data_ymax), max(fits_ymax)), 
                         ymin=pmin(min(data_ymin), min(fits_ymin)),
                         xmax=pmax(max(data_xmax), max(fits_xmax)), 
                         xmin=pmin(min(data_xmin), min(fits_xmin)))
  } else {
    all_data <- all_data %>%
                  mutate(ymax=max(data_ymax), 
                         ymin=min(data_ymin),
                         xmax=max(data_xmax), 
                         xmin=min(data_xmin))
  }


  all_data <- all_data %>% 
   mutate(pretty_ymin=map2_dbl(ymin, ymax, 
                      ~ pretty_rng(.x, .y, nticks=ynticks)[1]), 
          pretty_ymax=map2_dbl(ymin, ymax,
                      ~ pretty_rng(.x, .y, nticks=ynticks)[2]))
  all_data$nplot <- 1:nrow(all_data)
  dims <- c(nunique(data[[row_var]]), nunique(data[[col_var]]))

  # figure out which panels have axes labels
  all_data <- all_data %>% 
    mutate(y_axis=as.integer(col_fct) == 1, x_axis=as.integer(row_fct)==dims[1])

  # now, create all plots
  npans <- prod(dims)
  layout_vec <- 1:npans
  layout_mat <- matrix(layout_vec, nrow=dims[1], ncol=dims[2], byrow=TRUE)
  layout_widths <- rep.int(1, ncol(layout_mat))
  layout_heights <- rep.int(1, nrow(layout_mat))
  if (!is.na(legend_title)) {
    layout_mat <- cbind(layout_mat, npans+1L)
  } else {
    legend_prop <- 0
  }
  # add plots for y/x labels
  layout_mat <- rbind(layout_mat, npans+2L)
  layout_mat <- cbind(npans+3L, layout_mat)
  layout_widths <- c(axis_label_prop, layout_widths, legend_prop)
  layout_heights <- c(layout_heights, axis_label_prop)

  # create dummy panels 
  bottom_left <- row(layout_mat) == nrow(layout_mat) & col(layout_mat) == 1
  layout_mat[bottom_left] <- npans+4L
  bottom_right <- (row(layout_mat) == nrow(layout_mat) & 
                   col(layout_mat) == ncol(layout_mat))
  layout_mat[bottom_right] <- npans+5L

  layout(layout_mat, widths=layout_widths, heights=layout_heights)
  par(
      # mar=c(0, 1.5, 1.5, 0.5),
      mar=c(0, 1, 1, 0.5),
      oma=c(5, 5, 2, 1),
      ps=fontsize)
  for (i in all_data$nplot) {
    plot.new()

    # set up coordinates.
    ymin <- all_data$ymin[i]
    ymax <- all_data$ymax[i]
    xmin <- all_data$xmin[i]
    xmax <- all_data$xmax[i]
    pretty_ymax <- all_data$pretty_ymax[i]
    pretty_ymin <- all_data$pretty_ymin[i]
    fac <- 0.08
    if (pretty_rng) {
      plot.window(xlim=c(xmin*(1-fac), xmax*(1+fac)), 
                  ylim=c(pretty_ymin*(1-fac), pretty_ymax*(1+fac)))
    } else {
      plot.window(xlim=c(xmin*(1-fac), xmax*(1+fac)), 
                  ylim=c(ymin*(1-fac), ymax*(1+fac)))
    }
    panel_data <- all_data$data[[i]]
    # group data
    data_grp <- panel_data[[grp_var]]
    if (has_fits) {
      panel_fits <- all_data$fits[[i]]
      fit_grp <- panel_fits[[grp_var]]
    }
    
    segments(xmin, 0, xmax, 0, col=dashed_gray, lty=2)
    single_type <- length(type) == 1 
    for (gi in seq_along(levels(data_grp))) {
      data_grp_i <- which(panel_data[[grp_var]] == levels(data_grp)[gi])
      pty <- type
      if (!single_type)
        pty <- type[gi]
      if (pty == 'p') {
        points(panel_data[[x_var]][data_grp_i], 
               panel_data[[y_var]][data_grp_i], 
               col=adjustcolor(grp_cols[data_grp], alpha.f=point_alpha),
               pch=19, cex=0.6, lwd=0)
      } else if (pty == 'l') {
        lines(panel_data[[x_var]][data_grp_i], 
               panel_data[[y_var]][data_grp_i], 
               col=grp_cols[gi],
               lwd=lwd)
      }

       if (has_fits) {
         fit_grp_i <- which(panel_fits[[grp_var]] == levels(fit_grp)[gi])
         lines(panel_fits[[x_var]][grp_i], 
               panel_fits[[y_var]][grp_i],
               col=grp_cols[gi], lwd=2)
       }
      if (!is.null(N)) {
        # hacky way to add decay of Va due to drift
        times <- panel_data[[x_var]][data_grp_i]
        Vas <- all_data[[row_var]][[i]]*exp(-times/(2*N))
        lines(times, Vas, col=dashed_gray, lwd=1, lty=2)
      }
    }

    panel_row <- all_data[[row_var]][[i]]
    panel_col <- all_data[[col_var]][[i]]

    corners <- par('usr')
    if (i <= dims[2]) {
      text(x=mean(corners[1:2]), y=corners[4]*(1+0.04),
           latex2exp:::TeX(paste0(col_lab, prettyNum(panel_col))),
           font=2, cex=col_lab_cex, xpd=NA, col=title_col)
    }
    if (i %% dims[2] == 0) {
      corners = par("usr")
      text(x = corners[2], y = mean(corners[3:4]), 
           latex2exp:::TeX(paste0(row_lab, prettyNum(panel_row))),
           srt = 270, font=2, cex=row_lab_cex, xpd=NA, col=title_col)
      # mtext(
           # latex2exp:::TeX(paste0(row_lab, prettyNum(panel_row))),
           # side=4)
    }
    # text(0.55*corners[2], 0.8*corners[4], adj=c(0, 0),
         # latex2exp::TeX(paste0('$V_A = $', panel_row, ", $r = $", panel_col))
         # )
    # title(latex2exp::TeX(paste0('$V_A = $', panel_row, ", $r = $", panel_col)), 
          # line=-0.5)

    # axes with labels
    labeled_x <- all_data$x_axis[[i]]
    labeled_y <- all_data$y_axis[[i]]
    xticks <- labeling::extended(xmin, xmax, xnticks)
    yticks <- labeling::extended(ymin, ymax, ynticks)
    axis(1, at=xticks, labels=if (labeled_x) prettyNum(xticks) else FALSE,
         col.axis=axs_col, col=axs_col, 
         line=0.4,
         lwd=axis_lwd,
         padj=-0.9,
         tck=0.018)
    axis(2, at=yticks, labels=if (labeled_y) prettyNum(yticks) else FALSE,
         col.axis=axs_col, col=axs_col, las=1,
         line=0.4,
         hadj=0.9,
         lwd=axis_lwd,
         tck=0.018)

  }
  plot.new()
  par(mar=c(0, 0, 0, 0))
  legend('center', inset=0, legend=grp_levels, fill=grp_cols,
         title=legend_title,
         # adj=c(0, 0.45),
         cex=1.1,
         bty='n', text.col=title_col,
         border=0)
  plot.new()
  mtext(xlab, 1, line=x_axis_line, col=title_col)
  plot.new()
  mtext(ylab, 2, line=y_axis_line, col=title_col)

  # two dummy plots
  plot.new(); plot.new()

}



cumcov_panels <- function(x, include_data=FALSE, cex=1, N=1e3, 
                          axs_col='gray12',
                          mtext_cex=1,
                          lwd=2,
                          lx=2, ly=2,
                          include_pred=FALSE,
                          ymax=NULL,
                          basic=TRUE,
                          pt_rng_col='gray12',
                          ylab="Var$(p_{t} - p_{0})/t$") {
  va_levels <- sort(unique(x$Va))

  cols <- wesanderson::wes_palette("Darjeeling1")
  title_col <- "grey10"
  subtitle_col <- "grey30"  # bold font, lighter color looks better

  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  ncols <- length(unique(va_levels))
  par(mfrow=c(1, ncols),
      mar=c(0, .1, 1.5, 0),
      oma=c(3, 4, 0.5, 4.5),
      cex=cex)
  if (is.null(ymax)) {
    ymax <- x %>% mutate(tot = pred_cov+pred_var) %>% 
      pull(tot) %>% max(na.rm=TRUE)
  } 

  for (va_i in seq_along(va_levels)) {
      labs <- latex2exp::TeX(paste(rep(c('cov with', 'var with'), 2), 
                                   rep(c('V_A', 'neutral SSH'), each=2)))
    va_level <- va_levels[va_i]
    pd <- x %>% filter(Va == va_level) 
    if (basic) {
      pred <- rbind(pd$pred_var, pd$pred_cov)
      colnames(pred) <- pd$genlen
    } else {
      pred_zvar <- with(pd, rbind(pd$pred_var_zvar, pd$pred_cov_zvar))
      pred_genic <- with(pd, rbind(pd$pred_var_ssh, pd$pred_cov_ssh))
      pred_zvarb <- rbind(pred_zvar, matrix(0, ncol=ncol(pred_genic), 
                                            nrow=nrow(pred_genic)))
      pred_genicb <- rbind(matrix(0, ncol=ncol(pred_zvar), 
                                  nrow=nrow(pred_zvar)), pred_genic) 
      m <- cbind(pred_zvarb, pred_genicb)
      pred <- m[,c(1,5,2,6,3,7,4,8)]
      colnames(pred) <- paste0(rep(c('', ''), 
                                   length(pd$genlen)), pd$genlen)
    }
    barx <- barplot(pred, ylim=c(0, ymax*1.1), 
                    col=cols[c(5, 3, 2,4)],
                    cex.axis=0.7*cex,
                    col.axis=axs_col,
                    xaxt='n',
                    border=0,
                    axes=FALSE, las=2, cex.names=0.7*cex)
    mtext(latex2exp:::TeX(ylab), side=2, 
          col=title_col,      
          cex=mtext_cex,
          line=2.4, outer=TRUE)
    mtext("region length (Morgans)", side=1, 
          cex=mtext_cex,
          col=title_col,      
          line=1.3, outer=TRUE)

    if (include_data)  {
      this_dsum <- x %>% filter(Va==va_level) 
      if (!basic) 
        barx <- (barx[seq(1, length(barx), 2)] + barx[seq(2, length(barx), 2)])/2
      points(barx, this_dsum$cov_mean, pch=19, col=pt_rng_col, cex=0.4)
      points(barx, this_dsum$var_mean, pch=19, col=pt_rng_col, cex=0.4)
      segments(barx, this_dsum$cov_lower, barx, this_dsum$cov_upper, 
               col=pt_rng_col)
      segments(barx, this_dsum$var_lower, barx, this_dsum$var_upper, 
               col=pt_rng_col)
    }
    twoN <- 2*N  # if later, we want to add drift expectation
    axis(1, barx, pd$genlen, las=1, tick=FALSE, 
         line=-0.9,
         cex.axis=0.8*cex, 
         col.axis=axs_col)
    if (va_i == 1) {
      axis(2, col.axis=axs_col, cex.axis=cex*0.8, col=axs_col, 
           lwd=lwd, tck=0.018, las=1, hadj=0.6)
    }
 
    mtext(latex2exp:::TeX(sprintf('$V_a = %0.2g$', va_level)), 
          line=0,
          col=axs_col,
          side=3, cex=mtext_cex, las=1)

    if (FALSE || !basic && va_i == 1) {
      legend(lx, ly, legend=labs, fill=cols[c(3,5,4,2)], 
             #cex=0.86*cex,
             cex=0.9,
             x.intersp=0.6,
             adj=c(0, 0.45),
             bty='n', text.col=axs_col,
             border=0,
             xpd=NA) 
    }

  }
  if (basic) {
    legend(lx, ly, legend=c('covariance', 'variance'), fill=cols[c(3, 5)], 
           #cex=0.86*cex,
           cex=0.9,
           x.intersp=0.6,
           adj=c(0, 0.45),
           bty='n', text.col=axs_col,
           border=0,
           xpd=NA)
  }

   
}




# rec_params <- c(0.05, 0.5, 0.1, 1.5)
# va_params <- c(0.01, 0.02, 0.05)


# rec_params <- c(0.5, 0.1, 1.5)
# va_params <- c(0.01, 0.02, 0.05)

# pd5f <- pd5 %>% mutate(gen=t1) %>% filter(Va %in% va_params, genlen %in% rec_params)
# predf <- pred %>% filter(Va %in% va_params, genlen %in% rec_params)

# panel_plot(pd5f, predf, Va, genlen, L, free_y=FALSE, row_lab='V_A', col_lab='r')

# ggplot() +  geom_point(data=pd5f, aes(t1, cov, color=as.factor(L)), 
#                        alpha=0.5, size=0.8) +
#   geom_hline(yintercept=0, color='red') + 
#   geom_line(data=predf, aes(gen, pred/2, color=as.factor(L))) + 
#   facet_grid(Va ~ genlen, scales='free_y') 

simple_theme <- function(base_size = 11, base_family = "", ticks = TRUE) {
  ## TODO: start with theme_minimal
  ret <- theme_classic(base_family = base_family, base_size = base_size) +
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line=element_line(size = 0.3, colour = "grey38"),
          axis.ticks=element_line(size = 0.3, colour = "grey38"),
          plot.background = element_blank(),
          strip.background=element_rect(linetype=0),
          strip.text=element_text(size=12),
          panel.grid = element_blank())
  if (!ticks) {
    ret <- ret + theme(axis.ticks = element_blank())
  }
  ret
}


