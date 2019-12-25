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
#' @import dplyr
#'

format_features <- function(dat,
                            features,
                            var_type) {

  # Listify, add names
  if (features %>% is.data.frame) {
    features <- as.list(features)
  } else if (!(features %>% is.list)) {
    if (var_type == 'Categorical') {
      features <- list('Group' = features)
    } else if (var_type == 'Continuous') {
      features <- list('Variable' = features)
    }
  }
  if (var_type == 'Categorical') {
    if (names(features) %>% is.null) {
      if (length(features) == 1L) {
        names(features) <- 'Group'
      } else {
        names(features) <- paste('Factor', seq_along(features))
      }
    }
    features <- features %>% map(as.factor)
  } else if (var_type == 'Continuous') {
    if (names(features) %>% is.null) {
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
        stop('Categorical feature ', x, ' not equal to sample size.')
      }
      if (length(levels(features[[x]])) == 1L) {
        stop('Categorical feature ', x, ' is invariant.')
      }
    } else if (var_type == 'Continuous') {
      if (length(features[[x]]) != ncol(dat)) {
        stop('Continuous feature ', x, ' not equal to sample size.')
      }
      if (var(features[[x]]) == 0L) {
        stop('Continuous feature ', x, ' is invariant.')
      }
      if (!(features[[x]] %>% is.numeric)) {
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
#' binomially distributed data.
#'
#' @param vec Vector of observed outcomes, or one or several vectors of
#'   predicted values.
#' @param vec_type Character string specifying whether vector is \code{"obs"}
#'   or \code{"pred"}.
#' @param n Number of events for which predictions are required.
#'
#' @importFrom purrr every

format_binom <- function(vec,
                         vec_type,
                         n) {

  # Format obs
  if (vec_type == 'obs') {
    if (!(every(vec, is.finite))) {
      stop('Missing or infinite values detected in obs.')
    }
    if (vec %>% is.character) {
      vec <- vec %>% as.factor(.)
    }
    if (is.factor(vec)) {
      if (length(levels(vec)) != 2L) {
        stop('Response must be dichotomous.')
      } else {
        warning('A positive outcome is hereby defined as obs == "',
                levels(vec)[1], '". To change this to obs == "', levels(vec)[2],
                '", either relevel the factor or recode response as numeric ',
                '(1/0).')
        vec <- ifelse(vec == levels(vec)[1L], 1L, 0L)
      }
    }
    if (vec %>% is.logical) {
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
    if (vec %>% is.data.frame) {
      vec <- vec %>% as.list(.)
    } else if (!(vec %>% is.list)) {
      vec <- list(vec)
    }
    if (names(vec) %>% is.null) {
      names(vec) <- paste0('M', seq_along(vec))
    }
    for (m in seq_along(vec)) {
      if (!every(vec[[m]], is.finite)) {
        stop('Missing or infinite values detected in pred.')
      }
      if (!(vec[[m]] %>% is.numeric)) {
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
#' @param pal Palette name or constituents.
#' @param var_type Character string specifying whether features are \code{
#'   "Categorical"} or \code{"Continuous"}.
#' @param n Number of unique groups for which colors must be assigned. Only
#'   relevant if \code{var_type = "Categorical"}.
#'
#' @importFrom purrr map_lgl
#' @importFrom scales hue_pal
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis.map
#' @import ggsci
#' @import dplyr
#'

colorize <- function(pal,
                     var_type,
                     n) {

  # Preliminaries
  qual_pals <- c('ggplot', 'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2',
                 'Set1', 'Set2', 'Set3', 'npg', 'aaas', 'nejm', 'lancet', 'jco',
                 'ucscgb', 'd3', 'locuszoom', 'igv', 'uchicago', 'startrek',
                 'futurama', 'rickandmorty', 'simpsons', 'gsea')
  seq_pals <- c('Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
                'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
                'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
  div_pals <- c('BrBg', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu',
                'RdYlGn', 'Spectral', 'viridis', 'magma', 'plasma', 'inferno',
                'cividis')
  if (!(all(map_lgl(pal, is.color)))) {
    if (length(pal) > 1L) {
      stop('When passing individual colors to define a palette, each must ',
           'denote a valid color in R.')
    } else {
      if (var_type == 'Categorical' && !pal %in% qual_pals) {
        stop('Invalid palette for categorical features. Must be either a ',
             'vector of valid R colors or one of the following palettes: ',
             stringify(group_pals, 'or'), '.')
      }
      if (var_type == 'Continuous' && !pal %in% c(seq_pals, div_pals)) {
        stop('Invalid palette for continuous features. Must be either a ',
             'vector of valid R colors or one of the following palettes: ',
             stringify(c(seq_pals, div_pals), 'or'), '.')
      }
    }
  } else {
    if (var_type == 'Categorical') {
      if (length(pal) < n) {
        stop('Insufficient colors in palette. ',
             'Need at least ', n, ' for this plot.')
      } else if (length(pal) < n) {
        warning('Too many colors in palette. ',
                'Only the first ', n, ' will be used.')
      }
      out <- pal
    } else {
      out <- colorRampPalette(pal)(256)
    }
  }

  # Presets
  if (pal %in% qual_pals) {
    if (pal == 'ggplot') {
      out <- hue_pal()(n)
    } else if (pal %in% c('Accent', 'Dark2', 'Paired', 'Pastel1',
                          'Pastel2', 'Set1', 'Set2', 'Set3')) {
      out <- brewer.pal(n, pal)
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
  } else if (pal %in% seq_pals) {
    out <- colorRampPalette(brewer.pal(9, pal))(256)
  } else if (pal %in% div_pals) {
    if (pal %in% c('viridis', 'magma', 'plasma', 'inferno', 'cividis')) {
      if (pal == 'viridis') {
        df <- viridis.map %>% filter(opt == 'A')
      } else if (pal == 'magma') {
        df <- viridis.map %>% filter(opt == 'B')
      } else if (pal == 'plasma') {
        df <- viridis.map %>% filter(opt == 'C')
      } else if (pal == 'inferno') {
        df <- viridis.map %>% filter(opt == 'D')
      } else if (pal == 'cividis') {
        df <- viridis.map %>% filter(opt == 'E')
      }
      out <- rgb(df$R, df$G, df$B)
    } else {
      out <- colorRampPalette(brewer.pal(11, pal))(256) %>% rev(.)
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
#' @import dplyr
#'

matrixize <- function(dat) {

  # Transform, extract
  if (dat %>% is('DGEList')) {
    keep <- rowSums(dat$counts) > 1L             # Minimal count filter
    dat <- dat[keep, ]
    nf <- dat$samples$norm.factors               # Calculate size factors
    if (nf %>% is.null || all(nf == 1L)) {
      dat <- calcNormFactors(dat)
    }
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (dat %>% is('DESeqDataSet')) {
    require(DESeq2)
    if (sizeFactors(dat) %>% is.null && 
        normalizationFactors(dat) %>% is.null) {
      dat <- estimateSizeFactors(dat)            # Normalize counts
    }
    dat <- counts(dat, normalized = TRUE)
    keep <- rowMeans(dat) > 0L                   # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (dat %>% is('DESeqTransform')) {
    require(SummarizedExperiment)
    dat <- assay(dat)
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }
  if (rownames(dat) %>% is.null) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (colnames(dat) %>% is.null) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }

  # Output
  return(dat)

}

#' Filter by Variance
#'
#' This utility function filters a matrix by probewise variance or median
#' absolute deviation.
#'
#' @param dat Omic data matrix.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for distance calculations.
#' @param robust Use robust variance statistic?
#'
#' @importFrom Rfast rowMads rowVars
#'

var_filt <- function(dat,
                     top,
                     robust) {

  # Preliminaries
  if (top > 1L) {
    if (top > nrow(dat)) {
      warning('top exceeds nrow(dat), at least after removing probes with ',
              'missing values and/or applying a minimal expression filter. ',
              'Proceeding with the complete ', nrow(dat), ' x ', ncol(dat),
              ' matrix.')
    }
  } else {
    top <- round(top * nrow(dat))
  }

  # Filter
  if (robust) {
    stats <- rowMads(dat)
  } else {
    stats <- rowVars(dat)
  }
  keep <- order(stats, decreasing = TRUE)[seq_len(min(top, nrow(dat)))]
  dat <- dat[keep, , drop = FALSE]

  # Output
  return(dat)

}

#' Create Distance Matrix
#'
#' This utility function calculates distance based on a user defined measure
#' and optionally filters probes by leading log fold change.
#'
#' @param dat Omic data matrix.
#' @param dist Distance measure to be used.
#' @param p Power of the Minkowski distance.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for distance calculations.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}.
#' @param center Center each probe prior to computing distances?
#' @param robust Use robust probe centering?
#'
#' @details
#' Data are optionally centered by probe and samplewise distance calculated 
#' using one of the following methods:
#'
#' \itemize{
#'   \item \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, and \code{
#'     "minkowski"} are all documented in the \code{\link[stats]{dist}}
#'     function. \code{bioplotr} relies on a lower level implementation via
#'     \code{Rfast::\link[Rfast]{Dist}} to speed up computations.
#'   \item \code{"cosine"} and \code{"canberra"} are implemented and documented
#'     in \code{wordspace::\link[wordspace]{dist.matrix}}.
#'   \item \code{"pearson"}, \code{"kendall"}, and \code{"spearman"} correspond
#'     to various forms of correlation distance, generally defined as 1 -- the
#'     absolute value of the correlation coefficient. See 
#'     \code{\link[stats]{cor}} for more details.
#'   \item \code{"bray"}, \code{"kulczynski"}, \code{"jaccard"}, \code{"gower"},
#'     \code{"altGower"}, \code{"morisita"}, \code{"horn"}, \code{"mountford"},
#'     \code{"raup"}, \code{"binomial"}, \code{"chao"}, \code{"cao"}, and
#'     \code{"mahalanobis"} are all available and documented in the \code{
#'     vegan::\link[vegan]{vegdist}} function. These are designed for use with
#'     ecological data, e.g. a matrix of microbial OTU counts.
#'   \item \code{"bhattacharyya"}, \code{"hellinger"}, \code{"kullback_leibler"},
#'     and \code{"MI"} are information theoretic distance metrics. The former
#'     three are implemented in \code{Rfast::\link[Rfast]{Dist}}. See
#'     \code{bioDist::\link[bioDist]{MIdist}} for details on the latter.
#' }
#'
#' @importFrom Rfast rowMedians rowmeans Dist
#' @importFrom wordspace dist.matrix
#' @importFrom vegan vegdist
#' @import dplyr
#'

dist_mat <- function(dat,
                     dist,
                     p,
                     top,
                     filter_method,
                     center,
                     robust = FALSE) {

  # Preliminaries
  pow <- p
  n <- ncol(dat)
  p <- nrow(dat)
  if (!(top %>% is.null)) {
    dat <- var_filt(dat, top, robust = FALSE)
  }

  # Center probes?
  if (center) {
    if (robust) {
      dat <- dat - rowMedians(dat)
    } else {
      dat <- dat - rowmeans(dat)
    }
  }

  # Create distance matrix
  if (top %>% is.null || filter_method == 'common') {
    if (dist %in% c('euclidean', 'manhattan', 'minimum', 'maximum', 'minkowski',
                    'bhattacharyya', 'hellinger', 'kullback_leibler')) {
      dm <- Dist(t(dat), method = dist, p = pow)
    } else if (dist %in% c('canberra', 'cosine')) {
      dm <- dist.matrix(dat, method = dist, byrow = FALSE)
    } else if (dist %in% c('bray', 'kulczynski', 'jaccard', 'gower', 'altGower',
                           'morisita', 'horn', 'mountford', 'raup' , 'binomial',
                           'chao', 'cao', 'mahalanobis')) {
      dm <- as.matrix(vegdist(t(dat), method = dist))
    } else if (dist %in% c('pearson', 'kendall', 'spearman')) {
      dm <- 1L - abs(cor(dat, method = dist))
    } else if (dist == 'MI') {
      require(bioDist)
      dm <- as.matrix(MIdist(t(dat)))
    }
  } else if (!(top %>% is.null)) {
    dm <- matrix(nrow = n, ncol = n)
    for (i in 2L:ncol(dat)) {
      for (j in 1L:(i - 1L)) {
        top_idx <- nrow(dat) - top + 1L
        if (dist == 'euclidean') {
          dm[i, j] <- sqrt(sum(sort.int((dat[, i] - dat[, j])^2L,
                                        partial = top_idx)[top_idx:p]))
        } else if (dist == 'maximum') {
          dm[i, j] <- max(sort.int(abs(dat[, i] - dat[, j]),
                                   partial = top_idx)[top_idx:p])
        } else if (dist == 'manhattan') {
          dm[i, j] <- sum(sort.int(abs(dat[, i] - dat[, j]),
                                   partial = top_idx)[top_idx:p])
        } else if (dist == 'minkowski') {
          dm[i, j] <- (sum(sort.int(abs(dat[, i] - dat[, j]),
                                    partial = top_idx)[top_idx:p]^pow))^(1L / pow)
        } else {
          tops <- order(abs(dat[, i] - dat[, j]), decreasing = TRUE)[seq_len(top)]
          m <- dat[tops, c(i, j)]
          if (dist %in% c('pearson', 'kendall', 'spearman')) {
            dm[i, j] <- max(1L - abs(cor(m, method = dist)))
          } else if (dist %in% c('bhattacharyya', 'hellinger',
                                 'total_variation', 'kullback_leibler')) {
            dm[i, j] <- max(Dist(t(m), method = dist))
          } else if (dist %in% c('canberra', 'cosine')) {
            dm[i, j] <- max(dist.matrix(m, method = dist, byrow = FALSE))
          } else if (dist %in% c('bray', 'kulczynski', 'jaccard', 'gower',
                                 'altGower', 'morisita', 'horn', 'mountford',
                                 'raup', 'binomial', 'chao', 'cao',
                                 'mahalanobis')) {
            dm[i, j] <- max(vegdist(t(m), method = dist))
          } else if (dist == 'MI') {
            require(bioDist)
            dm[i, j] <- max(MIdist(t(m)))
          }
        }
      }
    }
    dm <- pmax(dm, t(dm), na.rm = TRUE)
    diag(dm) <- 0L
  }

  # Output
  return(dm)

}

#' Size Data Points
#'
#' This utility function assigns a reasonable default point size parameter based 
#' on data dimensionality.
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
    out <- 0.5
  } else {
    out <- 0.25
  }
  return(out)
}

#' Set Data Point Transparency
#'
#' This utility function assigns a reasonable default alpha parameter based on 
#' data dimensionality.
#'
#' @param df Data frame to be passed to \code{ggplot}.
#'

pt_alpha <- function(df) {
  if (nrow(df) <= 20L) {
    out <- 1L
  } else if (nrow(df) <= 100L) {
    out <- 0.9
  } else if (nrow(df) <= 1e4L) {
    out <- 0.75
  } else {
    out <- 0.5
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
  if (!(loc %>% is.null)) {
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
  }


  # Output
  if (!hover) {
    print(p)
  } else {
    suppressPackageStartupMessages(require(plotly))
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
                  size,
                  alpha,
                  title,
                  xlab,
                  ylab,
                  legend,
                  hover,
                  D3) {

  if (size %>% is.null) {
    size <- pt_size(df)
  }
  if (alpha %>% is.null) {
    alpha <- pt_alpha(df)
  }
  if (!D3) {
    p <- ggplot(df, aes(PC1, PC2)) +
      geom_hline(yintercept = 0L, color = 'grey') +
      geom_vline(xintercept = 0L, color = 'grey') +
      labs(title = title, x = xlab, y = ylab) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (c(group, covar) %>% is.null) {
      p <- p + geom_text(aes(label = Sample), alpha = alpha,
                         hjust = 'inward', vjust = 'inward')
    } else if (length(c(group, covar)) == 1L) {
      if (label) {
        p <- p + geom_text(aes(label = Sample, color = Feature1), alpha = alpha,
                           hjust = 'inward', vjust = 'inward') +
          labs(color = feature_names[1L])
      } else {
        if (!(covar %>% is.null)) {
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
      if (covar %>% is.null) {
        p <- p + guides(color = guide_legend(order = 1L),
                        shape = guide_legend(order = 2L))
      } else {
        p <- p + guides(color = guide_colorbar(order = 1L),
                        shape = guide_legend(order = 2L))
      }
    }
    if (covar %>% is.null) {
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
#' @param pal Palette name or constituents.
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
    n.cols <- seq_along(features) %>%
      map_dbl(~ length(levels(features[[.x]])))
    pal <- colorize(pal, 'Categorical', sum(n.cols))
    cols <- split(pal, rep(seq_along(features), n.cols))
  } else if (var_type == 'Continuous') {
    if (length(pal) != length(features)) {
      stop('Number of palettes passed to pal_covar must match number of ',
           'continuous features passed to covar.')
    } else if (pal %>% is.list) {
      cols <- seq_along(features) %>%
        map(~ colorize(pal[[.x]], 'Continuous'))
    } else {
      cols <- seq_along(features) %>%
        map(~ colorize(pal[.x], 'Continuous'))
    }
  }
  names(cols) <- NULL

  # Output
  return(cols)

}

#' Format String
#'
#' This utility function formats a character vector into a string representing
#' a well formed English list, with quotations around each element and the
#' appropriate conjunction between the penultimate and final elements.
#'
#' @param x Character vector.
#'

stringify <- function(x,
                      conjunct) {
  n <- length(x)
  if (n > 1) {
    x <- paste0('"', x)
    x <- paste0(x, '"')
    x <- c(x[seq_len(n - 1L)], conjunct, x[n])
    if (n > 2L) {
      x <- paste(x, sep = '', collapse = ', ')
      if (conjunct == 'and') {
        x <- gsub('and,', 'and', x)
      } else {
        x <- gsub('or,', 'or', x)
      }
    } else {
      x <- paste(x, sep = '', collapse = ' ')
    }
  }
  return(x)
}


