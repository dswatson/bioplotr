#' Format Features
#'
#' This utility function formats categorical and continuous features for plots.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param features Vector of length equal to sample size, or several such
#'   vectors organized into a list or data frame.
#' @param var_type Character string specifying whether features are \code{
#'   "Categorical"} or \code{"Continuous"}.
#'
#' @importFrom purrr map
#'

format_features <- function(dat,
                            features,
                            var_type) {

  # Listify, add names
  if (is.data.frame(features)) {
    features <- as.list(features)
  } else if (!is.list(features)) {
    if (var_type == 'Categorical') {
      features <- list('Group' = features)
    } else if (var_type == 'Continuous') {
      features <- list('Variable' = features)
    }
  }
  if (var_type == 'Categorical') {
    if (is.null(names(features))) {
      if (length(features) == 1L) {
        names(features) <- 'Group'
      } else {
        names(features) <- paste('Factor', seq_along(features))
      }
    }
    features <- map(features, as.factor)
  } else if (var_type == 'Continuous') {
    if (is.null(names(features))) {
      if (length(features) == 1L) {
        names(features) <- 'Variable'
      } else {
        names(features) <- paste('Variable', seq_along(features))
      }
    }
  }

  # Check dimensions, invariance, object class
  for (x in seq_along(features)) {
    if (var_type == 'Categorical') {
      if (length(features[[x]]) != ncol(dat)) {
        stop(paste('Categorical feature', x, 'not equal to sample size.'))
      }
      if (length(levels(features[[x]])) == 1L) {
        stop(paste('Categorical feature', x, 'is invariant.'))
      }
    } else if (var_type == 'Continuous') {
      if (length(features[[x]]) != ncol(dat)) {
        stop(paste('Continuous feature', x, 'not equal to sample size.'))
      }
      if (var(features[[x]]) == 0L) {
        stop(paste('Continuous feature', x, 'is invariant.'))
      }
      if (!is.numeric(features[[x]])) {
        stop('All variables passed to covar must be of class numeric.')
      }
    }
  }

  # Output
  return(features)

}

#' Format Binomial Things
#'
#' This utility function formats vectors of observed and predicted values for
#' binomially distributed events.
#'
#' @param vec Vector of observed outcomes, or one or several vectors of
#'   predicted values.
#' @param vec_type Character string specifying whether vector is \code{"obs"}
#'   or \code{"pred"}.
#' @param n Number of events for which predictions are required.
#'

format_binom <- function(vec,
                         vec_type,
                         n) {

  # Format obs
  if (vec_type == 'obs') {
    if (any(!is.finite(vec))) {
      stop('Missing or infinite values detected in obs.')
    }
    if (is.character(vec)) {
      vec <- as.factor(vec)
    }
    if (is.factor(vec)) {
      if (length(levels(vec)) != 2L) {
        stop('Response must be dichotomous.')
      } else {
        warning(paste0('A positive outcome is hereby defined as obs == "',
                       levels(vec)[1], '". To change this to obs == "',
                       levels(vec)[2], '", either relevel the factor or ',
                       'recode response as numeric (1/0).'))
        vec <- ifelse(vec == levels(vec)[1], 1L, 0L)
      }
    }
    if (is.logical(vec)) {
      vec <- ifelse(vec, 1L, 0L)
    }
    if (!all(vec %in% c(0L, 1L))) {
      stop('A numeric response can only take values of 1 or 0.')
    }
    if (var(vec) == 0L) {
      stop('Response is invariant.')
    }

  # Format pred
  } else if (vec_type == 'pred') {
    if (is.data.frame(vec)) {
      vec <- as.list(vec)
    } else if (!is.list(vec)) {
      vec <- list(vec)
    }
    if (is.null(names(vec))) {
      names(vec) <- paste0('M', seq_along(vec))
    }
    for (m in seq_along(vec)) {
      if (any(!is.finite(vec[[m]]))) {
        stop('Missing or infinite values detected in pred.')
      }
      if (!is.numeric(vec[[m]])) {
        stop('pred must be a numeric vector, or several such vectors ',
             'organized into a list or data frame.')
      }
      if (n != length(vec[[m]])) {
        stop('obs and pred vectors must be of equal length.')
      }
    }
  }

  # Output
  return(vec)

}

#' Test for Valid Color Representation
#'
#' This utility function checks whether each string in a character vector
#' denotes a valid color in R.
#'
#' @param chr A character vector.
#'

is.color <- function(chr) {
  sapply(chr, function(x) {
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)
  })
}

#' Assign Colors
#'
#' This utility function applies a user selected color palette to a ggplot
#' image.
#'
#' @param pal Color palette.
#' @param var_type Character string specifying whether features are categorical
#'   or continuous.
#' @param n Number of unique groups for which colors must be assigned. Only
#'   relevant if \code{var_type = "Categorical"}.
#'
#' @importFrom scales hue_pal
#' @importFrom RColorBrewer brewer.pal
#' @import ggsci
#'

colorize <- function(pal,
                     var_type,
                     n) {

  # Preliminaries
  group_pals <- c('ggplot', 'npg', 'aaas', 'nejm', 'lancet', 'jco', 'ucscgb',
                  'd3', 'locuszoom', 'igv', 'uchicago', 'startrek', 'futurama',
                  'rickandmorty', 'simpsons', 'gsea')
  covar_pals <- c('blues', 'greens', 'purples', 'greys', 'oranges', 'reds')
  if (!all(is.color(pal))) {
    if (length(pal) > 1L) {
      stop('When passing multiple strings to define a color palette, ',
           'each must denote a valid color in R.')
    } else {
      if (var_type == 'Categorical' & !pal %in% group_pals) {
        stop(paste('Invalid palette for categorical features. Must be either a',
                   'vector of valid R colors or one of the following palettes:',
                   group_pals))
      }
      if (var_type == 'Continuous' & !pal %in% covar_pals) {
        stop(paste('Invalid palette for continuous features. Must be either a',
                   'vector of valid R colors or one of the following palettes:',
                   covar_pals))
      }
    }
  } else {
    if (var_type == 'Categorical' & length(pal) < n) {
      stop(paste('Insufficient colors in palette. Need at least', n, 'for',
                 'this plot.'))
    }
    if (var_type == 'Categorical' & length(pal) > n) {
      warning(paste('Too many colors in palette. Only the first', n, 'will',
                    'be used.'))
    }
    out <- pal
  }

  # Presets
  if (var_type == 'Categorical') {
    if (pal == 'ggplot') {
      out <- hue_pal()(n)
    } else if (pal == 'npg') {
      out <- pal_npg()(n)
    } else if (pal == 'aaas') {
      out <- pal_aaas()(n)
    } else if (pal == 'nejm') {
      out <- pal_nejm()(n)
    } else if (pal == 'lancet') {
      out <- pal_lancet()(n)
    } else if (pal == 'jco') {
      out <- pal_jco()(n)
    } else if (pal == 'ucscgb') {
      out <- pal_ucscgb()(n)
    } else if (pal == 'd3') {
      out <- pal_d3()(n)
    } else if (pal == 'locuszoom') {
      out <- pal_locuszoom()(n)
    } else if (pal == 'igv') {
      out <- pal_igv()(n)
    } else if (pal == 'uchicago') {
      out <- pal_uchicago()(n)
    } else if (pal == 'startrek') {
      out <- pal_startrek()(n)
    } else if (pal == 'futurama') {
      out <- pal_futurama()(n)
    } else if (pal == 'rickandmorty') {
      out <- pal_rickandmorty()(n)
    } else if (pal == 'simpsons') {
      out <- pal_simpsons()(n)
    }
  } else if (var_type == 'Continuous') {
    if (pal == 'blues') {
      out <- colorRampPalette(brewer.pal(9, 'Blues'))(n = 256)
    } else if (pal == 'greens') {
      out <- colorRampPalette(brewer.pal(9, 'Greens'))(n = 256)
    } else if (pal == 'purples') {
      out <- colorRampPalette(brewer.pal(9, 'Purples'))(n = 256)
    } else if (pal == 'greys') {
      out <- colorRampPalette(brewer.pal(9, 'Greys'))(n = 256)
    } else if (pal == 'oranges') {
      out <- colorRampPalette(brewer.pal(9, 'Oranges'))(n = 256)
    } else if (pal == 'reds') {
      out <- colorRampPalette(brewer.pal(9, 'Reds'))(n = 256)
    }
  }

  # Output
  return(out)

}

#' Standardize Matrix
#'
#' This utility function takes data objects from \code{limma}, \code{edgeR}, or
#' \code{DESeq2} pipelines and outputs a standard probe by sample matrix. Raw
#' counts are log2-CPM transformed with a warning.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#'
#' @importFrom edgeR cpm calcNormFactors
#' @importFrom limma getEAWP
#'

matrixize <- function(dat) {

  # Transform, extract
  if (is(dat, 'DGEList')) {
    keep <- rowSums(dat$counts) > 1L             # Minimal count filter
    dat <- dat[keep, ]
    if (is.null(dat$samples$norm.factors) |      # Calculate size factors
        all(dat$samples$norm.factors == 1L)) {
      dat <- calcNormFactors(dat)
    }
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqDataSet')) {
    require(DESeq2)
    if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
      dat <- estimateSizeFactors(dat)            # Normalize counts
    }
    dat <- counts(dat, normalized = TRUE)
    keep <- rowMeans(dat) > 0L                   # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqTransform')) {
    require(SummarizedExperiment)
    dat <- assay(dat)
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }

  # Output
  return(dat)

}

#' Create Distance Matrix
#'
#' This utility function calculates distance based on a user defined measure
#' and optionally filters probes by leading log fold change.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for distance calculations.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}.
#' @param dist Distance measure to be used. Currently supports \code{
#'   "euclidean", "pearson", "MI",} or \code{"KLD"}.
#'
#' @references
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., & Smyth, G.K.
#' (2015).
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/25605792}{limma powers differential
#' expression analyses for RNA-sequencing and microarray studies}. \emph{Nucleic
#' Acids Res.}, emph{43}(7): e47.
#'
#' @importFrom matrixStats rowVars
#' @importFrom wordspace dist.matrix
#' @importFrom bioDist MIdist KLdist.matrix
#' @importFrom KernSmooth dpih
#'

dist_mat <- function(dat,
                     top,
                     filter_method,
                     dist) {

  # Preliminaries
  if (!is.null(top)) {
    if (top > 1L) {
      if (top > nrow(dat)) {
        warning(paste(
          'top exceeds nrow(dat), at least after removing probes with missing',
          'values and/or applying a minimal expression filter. Proceeding with',
          'the complete', nrow(dat), 'x', ncol(dat), 'matrix.'
          ))
        top <- NULL
      }
    } else {
      top <- round(top * nrow(dat))
    }
  }

  # Variance filter?
  if (!is.null(top) & filter_method == 'common') {
    vars <- rowVars(dat)
    keep <- order(vars, decreasing = TRUE)[seq_len(min(top, nrow(dat)))]
    dat <- dat[keep, , drop = FALSE]
  }

  # Median center
  dat <- sweep(dat, 1, apply(dat, 1, median))

  # Create distance matrix
  if (is.null(top) | filter_method == 'common') {
    if (dist == 'euclidean') {
      dm <- dist.matrix(t(dat), method = 'euclidean')
    } else if (dist == 'pearson') {
      dm <- 1 - cor(dat)
    } else if (dist == 'MI') {
      dm <- as.matrix(MIdist(t(dat)))
    } else if (dist == 'KLD') {
      dm <- as.matrix(KLdist.matrix(t(dat)))
    }
  } else if (!is.null(top)) {
    dm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat))
    for (i in 2L:ncol(dat)) {
      for (j in 1L:(i - 1L)) {
        if (dist == 'euclidean') {
          top_idx <- nrow(dat) - top + 1L
          dm[i, j] <- sqrt(sum(sort.int((dat[, i] - dat[, j])^2L,
                                        partial = top_idx)[top_idx:nrow(dat)]))
        } else {
          tops <- order((dat[, i] - dat[, j])^2, decreasing = TRUE)[seq_len(top)]
          if (dist == 'pearson') {
            dm[i, j] <- 1 - cor(dat[tops, i], dat[tops, j])
          } else if (dist == 'MI') {
            dm[i, j] <- max(as.matrix(MIdist(t(dat[tops, c(i, j)]))))
          } else if (dist == 'KLD') {
            dm[i, j] <- max(as.matrix(KLdist.matrix(t(dat[tops, c(i, j)]))))
          }
        }
      }
    }
    dm <- pmax(dm, t(dm))
  }

  # Output
  return(dm)

}

#' Size Data Points
#'
#' This utility function assigns a reasonable point size parameter based on data
#' dimensionality.
#'
#' @param df Data frame to be passed to \code{ggplot}.
#'

pt_size <- function(df) {
  if (nrow(df) <= 10L) {
    out <- 3L
  } else if (nrow(df) <= 20L) {
    out <- 2L
  } else if (nrow(df) <= 100L) {
    out <- 1.5
  } else if (nrow(df) <= 1e3L) {
    out <- 1L
  } else if (nrow(df) <= 1e4L) {
    out <- 0.75
  } else if (nrow(df) <= 1e5L) {
    out <- 0.5
  } else {
    out <- 0.25
  }
  return(out)
}

#' Set Data Point Transparency
#'
#' This utility function assigns a reasonable alpha parameter based on data
#' dimensionality.
#'
#' @param df Data frame to be passed to \code{ggplot}.
#'

pt_alpha <- function(df) {
  if (nrow(df) <= 20L) {
    out <- 1L
  } else if (nrow(df) <= 100L) {
    out <- 0.85
  } else if (nrow(df) <= 1e5L) {
    out <- 0.5
  } else {
    out <- 0.25
  }
  return(out)
}

#' Output Image
#'
#' This utility function locates the legend (if supplied) and prints a ggplot or
#' ggplotly figure.
#'
#' @param p A \code{ggplot2} object.
#' @param hover Add text tooltip using plotly?
#' @param loc String specifying legend location.
#'
#' @importFrom ggplot2 theme
#'

gg_out <- function(p,
                   hover,
                   loc = NULL) {

  # Locate legend
  if (loc == 'right') {
    p <- p + theme(legend.position = 'right')
  } else if (loc == 'left') {
    p <- p + theme(legend.position = 'left')
  } else if (loc == 'top') {
    p <- p + theme(legend.position = 'top', legend.box = 'horizontal')
  } else if (loc == 'bottom') {
    p <- p + theme(legend.position = 'bottom', legend.box = 'horizontal')
  } else if (loc == 'bottomright') {
    p <- p + theme(legend.justification = c(0.99, 0.01),
                   legend.position = c(0.99, 0.01))
  } else if (loc == 'bottomleft') {
    p <- p + theme(legend.justification = c(0.01, 0.01),
                   legend.position = c(0.01, 0.01))
  } else if (loc == 'topright') {
    p <- p + theme(legend.justification = c(0.99, 0.99),
                   legend.position = c(0.99, 0.99))
  } else if (loc == 'topleft') {
    p <- p + theme(legend.justification = c(0.01, 0.99),
                   legend.position = c(0.01, 0.99))
  }

  # Output
  if (!hover) {
    print(p)
  } else {
    require(plotly)
    if (loc %in% c('right', 'left')) {
      p <- ggplotly(p, tooltip = 'text', height = 600, width = 650)
    } else if (loc %in% c('top', 'bottom')) {
      p <- ggplotly(p, tooltip = 'text', height = 650, width = 600)
    } else {
      p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    }
    print(p)
  }

}

#' Embed Manifold
#'
#' This utility function formats the output for \code{bioplotr}'s projection
#' functions (\code{plot_pca}, \code{plot_mds}, and \code{plot_tsne}).
#'
#' @param df Data frame formatted for \code{ggplot}ing.
#' @param group Optional factor variable formatted by \code{format_features}.
#' @param covar Optional continuous variable formatted by \code{format_features}.
#' @param group_cols Palette for \code{group} formatted by \code{colorize}.
#' @param covar_cols Palette for \code{covar} formatted by \code{colorize}.
#' @param feature_names Feature names.
#' @param label Label data points by sample name?
#' @param title Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param legend Legend location.
#' @param hover Show sample name by hovering mouse over data point?
#' @param D3 Render plot in three dimensions?
#'
#' @import ggplot2
#'

embed <- function(df,
                  group,
                  covar,
                  group_cols,
                  covar_cols,
                  feature_names,
                  label,
                  title,
                  xlab,
                  ylab,
                  legend,
                  hover,
                  D3) {

  size <- pt_size(df)
  alpha <- pt_alpha(df)
  if (!D3) {
    p <- ggplot(df, aes(PC1, PC2)) +
      geom_hline(yintercept = 0L, color = 'grey') +
      geom_vline(xintercept = 0L, color = 'grey') +
      labs(title = title, x = xlab, y = ylab) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (is.null(c(group, covar))) {
      p <- p + geom_text(aes(label = Sample), alpha = alpha,
                         hjust = 'inward', vjust = 'inward')
    } else if (length(c(group, covar)) == 1L) {
      if (label) {
        p <- p + geom_text(aes(label = Sample, color = Feature1), alpha = alpha,
                           hjust = 'inward', vjust = 'inward') +
          labs(color = feature_names[1L])
      } else {
        if (!is.null(covar)) {
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1),
                                size = size, alpha = alpha) +
              labs(color = feature_names[1L])
          )
        } else {
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature1),
                                size = size, alpha = alpha) +
              labs(color = feature_names[1L], shape = feature_names[1L])
          )
        }
      }
    } else {
      suppressWarnings(
        p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature2),
                            size = size, alpha = alpha) +
          labs(color = feature_names[1L], shape = feature_names[2L])
      )
      if (is.null(covar)) {
        p <- p + guides(color = guide_legend(order = 1L),
                        shape = guide_legend(order = 2L))
      } else {
        p <- p + guides(color = guide_colorbar(order = 1L),
                        shape = guide_legend(order = 2L))
      }
    }
    if (is.null(covar)) {
      p <- p + scale_color_manual(values = group_cols)
    } else {
      p <- p + scale_color_gradientn(colors = covar_cols)
    }
    gg_out(p, hover, legend)
  } else {
    ### REWRITE ###
    # symbls <- c(16, 17, 15, 3, 7, 8)           # This would be right if plotly worked
    symbls <- c(16, 18, 15, 3, 7, 8)
    require(plotly)
    p <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
                 text = ~Sample, color = ~Group, symbol = ~Group,
                 colors = pal_d3()(length(unique(df$Group))),
                 symbols = symbls[1:length(unique(df$Group))],
                 type = 'scatter3d', mode = 'markers',
                 alpha = 0.85, hoverinfo = 'text', marker = list(size = 5)) %>%
      layout(hovermode = 'closest', title = title, scene = list(
        xaxis = list(title = pve[min(pcs)]),
        yaxis = list(title = pve[other]),
        zaxis = list(title = pve[max(pcs)])))
    print(p)
  }

}


#' Format Annotation Track Colors
#'
#' This utility function formats the colors for plotting annotation tracks atop
#' heatmaps.
#'
#' @param features List output by \code{\link{format_features}}.
#' @param var_type Character string specifying whether features are categorical
#'   or continuous.
#'
#' @importFrom purrr map_dbl map
#' @importFrom RColorBrewer brewer.pal
#'

track_cols <- function(features,
                       pal,
                       var_type) {

  # Create color list
  if (var_type == 'Categorical') {
    n.cols <- map_dbl(seq_along(features), function(x) {
      length(levels(features[[x]]))
    })
    pal <- colorize(pal, 'Categorical', sum(n.cols))
    cols <- split(pal, rep(seq_along(features), n.cols))
  } else if (var_type == 'Continuous') {
    if (length(pal) != length(features)) {
      stop('Number of palettes passed to pal_covar must match number of
           continuous features passed to covar.')
    } else if (is.list(pal)) {
      cols <- map(seq_along(features), function(x) {
        colorize(pal[[x]], 'Continuous')
      })
    } else {
      cols <- map(seq_along(features), function(x) {
        colorize(pal[x], 'Continuous')
      })
    }
  }
  names(cols) <- NULL

  # Output
  return(cols)

}

#' Format String
#'
#' This utility function formats a character vector into a string representing
#' a well formed English list, with quotations around each element and an "and"
#' between the penultimate and final elements.
#'
#' @param x Character vector.
#'

stringify <- function(x) {
  n <- length(x)
  x <- paste0('"', x)
  x <- paste0(x, '"')
  x <- c(x[seq_len(n - 1L)], 'and', x[n])
  if (n > 2L) {
    x <- paste(x, sep = '', collapse = ', ')
    x <- gsub('and,', 'and', x)
  } else {
    x <- paste(x, sep = '', collapse = ' ')
  }
  return(x)
}


