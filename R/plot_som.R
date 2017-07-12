#' SOM
#'
#' This function plots the output of a self-organizing map.
#'
#' @param dat Either a probe by sample omic data matrix or an object of class
#'   \code{kohonen}.
#' @param type What should the plot visualize? Options include \code{"model"},
#'   \code{"sample"}, \code{"distance"}, \code{"counts"}, \code{"quality"}, and
#'   \code{"train"}. See Details.
#' @param design Design matrix for linear model with rows corresponding to
#'   samples and columns to model coefficients. Only relevant if \code{type =
#'   "model"}.
#' @param coef Column number or name specifying which coefficient or contrast
#'   of the linear model is of interest. Only relevant if \code{type = "model"}.
#' @param stat Which nodewise summary statistic should be plotted in the SOM?
#'   Options are \code{"lfc"} or \code{"t"}. Only relevant if \code{type =
#'   "model"}.
#' @param sample Column number or name specifying which sample should be
#'   plotted. Only relevant if \code{type = "sample"}.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable
#'   probes to be used for SOM.
#' @param grid_dims Vector of length two specifying x- and y-axis dimensions for
#'   the SOM grid. If \code{NULL}, values are chosen so as to create a square
#'   grid with approximately 10 probes in each node, presuming a uniform
#'   distribution across the map.
#' @param topo Train SOM using \code{"hexagonal"} or \code{"rectangular"}
#'   topology?
#' @param neighb Train SOM using a \code{"gaussian"} or \code{"bubble"}
#'   neighborhood function?
#' @param rlen Run length for SOM training.
#' @param parallel Allow for parallel computation of SOM?
#' @param pal_tiles String specifying the color palette to use for heatmap
#'   nodes. Defaults vary by \code{type} (see Details). Any divergent or
#'   sequential palette in \code{RColorBrewer} is acceptable. Alternatively,
#'   any vector of recognized colors may be supplied.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show constituent probes by hovering mouse over module tile? If
#'   \code{TRUE}, the plot is rendered in HTML and will either open in your
#'   browser's graphic display or appear in the RStudio viewer.
#' @param export Return fitted SOM model?
#' @param ... Additional arguments to be passed to \code{\link[limma]{lmFit}} if
#'   \code{type = "model"}.
#'
#' @details
#' A self-organizing map (SOM) is a type of neural network used for
#' dimensionality reduction and data visualization. SOMs come in many flavors,
#' but the kind used here clusters probes together in an unsupervised manner to
#' identify functionally relevant gene sets, examine patterns across samples,
#' and explore hypotheses. See the references below for more details on SOM
#' algorithms and their applications for omic research.
#'
#' If \code{type = "model"}, then a linear model is fit to the codebook vectors.
#' Units are colored either by each node's log fold change (if \code{stat =
#' "lfc"}) or moderated \emph{t}-statistic as calculated by \code{limma} (if
#' \code{stat = "t"}) for a given model coefficient. Default color palette is
#' \code{"Spectral"}.
#'
#' If \code{type = "sample"}, then the plot depicts that sample's expression
#' profile across the SOM. Default color palette is \code{"Spectral"}.
#'
#' If \code{type = "distance"}, then the figure renders the SOM's U-matrix,
#' which represents the Euclidean distance between codebook vectors for
#' neighboring units. This can be used to inspect for clusters and borders
#' in the SOM space. Default color palette is \code{"Greys"}.
#'
#' If \code{type = "counts"}, then units are colored by each node's probe
#' count. This distribution should ideally be nearly uniform. A large number of
#' empty units suggests that \code{grid_dim}s should be reduced; an uneven
#' spread across the units suggests that \code{grid_dim}s should be increased.
#' Default color palette is \code{"Blues"}.
#'
#' If \code{type = "quality"}, then nodes are colored by the mean distance of
#' intra-unit probes from one another. Very large values or uneven distributions
#' suggest the SOM may need more training iterations to reach a stable solution.
#' Default color palette is \code{"Reds"}.
#'
#' If \code{type = "train"}, then the plot displays the SOM's learning curve
#' over its \code{rlen} training iterations.
#'
#' Depending on the data's size, the model's parameters, and your computational
#' resources, it may take a considerable amount of time to fit a SOM to an
#' omic matrix. That is why \code{plot_som} returns the SOM object by default,
#' which may in turn be reused in subsequent calls to the function to explore
#' alternative aspects of the map. SOMs are trained using the batch algorithm
#' of \code{kohonen::\link[kohonen]{som}}. This may be executed in parallel for
#' faster mapping. See the package documentation for more details.
#'
#' @return
#' If \code{dat} is an omic data matrix or matrix-like object and \code{export
#' = TRUE}, then \code{plot_som} returns an object of class \code{kohonen}. This
#' is a fitted SOM as generated by the \code{kohonen} package. See \code{
#' \link[kohonen]{som}} for more details.
#'
#' @references
#' Kohonen, T. (1995). \emph{Self-Organizing Maps}. Berlin: Springer-Verlag.
#'
#' Nikkilä, J., Törönen, P., Kaski, S., Venna, J., Castrén, E., & Wong, G.
#' (2002).
#' \href{http://www.sciencedirect.com/science/article/pii/S0893608002000709}{
#' Analysis and visualization of gene expression data using Self-Organizing
#' Maps}. \emph{Neural Networks}, \emph{15}(8-9): 953-966.
#'
#' Wehrens, R. & Buydens, L.M.C. (2007). \href{http://bit.ly/2tBFE6R}{Self- and
#' Super-organizing Maps in R: The kohonen Package} \emph{J. Stat. Softw.},
#' \emph{21}(5).
#'
#' Wirth, H. (2012). \href{http://bit.ly/2tGq6iY}{Analysis of large-scale
#' molecular biological data using self-organizing maps}. Dissertation thesis,
#' University of Leipzig.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_som(mat, type = "sample", sample = 1)
#'
#' @seealso
#' \code{\link[kohonen]{som}}
#'
#' @export
#' @importFrom kohonen somgrid som unit.distances object.distances
#' @importFrom limma lmFit eBayes topTable
#' @importFrom purrr map_dbl map_chr
#' @import dplyr
#' @import ggplot2
#'

plot_som <- function(dat,
                     type = 'model',
                   design = NULL,
                     coef = NULL,
                     stat = 'lfc',
                   sample = 1,
                      top = 0.5,
                grid_dims = NULL,
                     topo = 'hexagonal',
                   neighb = 'gaussian',
                     rlen = 1000,
                 parallel = TRUE,
                pal_tiles = NULL,
                    title = NULL,
                   legend = 'right',
                    hover = FALSE,
                   export = TRUE, ...) {

  # Preliminaries
  if (dat %>% is('kohonen')) {
    y <- dat
    grid_dims <- c(y$grid$xdim, y$grid$ydim)
    export <- FALSE
  } else {
    if (ncol(dat) < 3L) {
      stop('dat includes only ', ncol(dat), ' samples; ',
           'need at least 3 to train a SOM.')
    }
  }
  if (type == 'model') {
    if (design %>% is.null || coef %>% is.null) {
      stop('design and coef must be supplied when type = "model".')
    }
    if (coef %>% is.character && !coef %in% colnames(design)) {
      stop('No coefficient named "', coef, '" found in design matrix.')
    }
    if (coef %>% is.numeric && coef > ncol(design)) {
      stop('No coef number ', coef, ' found in design matrix.')
    }
    if (!stat %in% c('lfc', 't')) {
      stop('stat must be "lfc" or "t".')
    }
  } else if (type == 'sample') {
    if (sample %>% is.null) {
      stop('sample must be specified when type = "sample".')
    }
    if (sample %>% is.numeric) {
      if ((dat %>% is('kohonen') && sample > ncol(dat$codes[[1]])) ||
          (!dat %>% is('kohonen') && sample > ncol(dat))) {
        stop('sample number exceeds sample size.')
      }
    } else {
      if ((dat %>% is('kohonen') && !sample %in% colnames(dat$codes[[1]])) ||
          (!dat %>% is('kohonen') && !sample %in% colnames(dat))) {
        stop('sample not found.')
      }
    }
  }
  if (!pal_tiles %>% is.null) {
    cols <- colorize(pal_tiles, var_type = 'Continuous')
  }

  # Tidy data
  if (!dat %>% is('kohonen')) {                  # Fit SOM
    dat <- matrixize(dat)
    if (!(top %>% is.null)) {                    # Filter by variance?
      dat <- var_filt(dat, top, robust = FALSE)
    }
    dat <- scale(dat, scale = FALSE)             # Mean center by sample
    if (grid_dims %>% is.null) {                 # ~10 probes/node
      grid_dims <- sqrt(nrow(dat) / 10L) %>% round(.)
    }
    som_grid <- kohonen::somgrid(grid_dims, grid_dims,
                                 neighbourhood.fct = neighb, topo = topo)
    if (parallel) {
      mode <- 'pbatch'
    } else {
      mode <- 'batch'
    }
    y <- kohonen::som(dat, grid = som_grid, rlen = rlen, mode = mode)
  }

  # Plot type?
  if (type == 'train') {
    df <- data_frame(Iteration = seq_len(rlen),
                      Distance = as.numeric(y$changes))
  } else {
    n_nodes <- nrow(y$grid$pts)
    if (type == 'model') {
      top <- lmFit(y$codes[[1]], design, ...) %>%
        eBayes(.) %>%
        topTable(coef = coef, number = Inf, sort.by = 'none')
      if (stat == 'lfc') {
        value <- top$logFC
        leg.txt <- expression(log[2]~'FC')
      } else if (stat == 't') {
        value <- top$t
        leg.txt <- expression('Moderated'~italic(t))
      }
      if (pal_tiles %>% is.null) {
        cols <- colorize('Spectral', var_type = 'Continuous')
      }
      if (title %>% is.null) {
        if (coef %>% is.character) {
          title <- paste('SOM Model:', coef)
        } else {
          title <- paste('SOM Model: Coefficient', coef)
        }
      }
    } else if (type == 'sample') {
      value <- y$codes[[1]][, sample]
      if (pal_tiles %>% is.null) {
        cols <- colorize('Spectral', var_type = 'Continuous')
      }
      leg.txt <- 'Expression'
      if (title %>% is.null) {
        title <- paste('SOM: Sample', sample)
      }
    } else if (type == 'quality') {
      value <- seq_len(n_nodes) %>%
        map_dbl(~ mean(y$distances[y$unit.classif == .x]))
      if (pal_tiles %>% is.null) {
        cols <- colorize('Reds', var_type = 'Continuous')
      }
      leg.txt <- 'Mean\nIntra-Unit\nDistance'
      if (title %>% is.null) {
        title <- 'SOM Quality'
      }
    } else if (type == 'distance') {
      node_dists <- unit.distances(y$grid)
      code_dists <- y %>%
        object.distances(type = 'codes') %>%
        as.matrix(.)
      code_dists[abs(node_dists - 1L) > .001] <- NA_real_
      value <- colMeans(code_dists, na.rm = TRUE)
      if (pal_tiles %>% is.null) {
        cols <- colorize('Greys', var_type = 'Continuous')
      }
      leg.txt <- 'Distance to\nNearest\nNeighbors'
      if (title %>% is.null) {
        title <- 'SOM U-Matrix'
      }
    } else if (type == 'counts') {
      value <- seq_len(n_nodes) %>%
        map_dbl(~ sum(y$unit.classif == .x))
      if (pal_tiles %>% is.null) {
        cols <- colorize('Blues', var_type = 'Continuous')
      }
      leg.txt <- 'Probes'
      if (title %>% is.null) {
        title <- 'SOM Node Size'
      }
    }
    probes <- seq_len(n_nodes) %>%
      map_chr(~ rownames(y$data[[1]])[y$unit.classif == .x] %>%
                paste(collapse = '\n'))
    probes[probes == ''] <- NA_character_
    df <- as_tibble(y$grid$pts) %>%
      mutate(Value = value,
            Probes = probes)
  }

  # Build Plot
  if (type == 'train') {
    p <- ggplot(df, aes(Iteration, Distance)) +
      geom_path() +
      labs(title = 'SOM Learning Curve',
               y = 'Mean Distance to Nearest Unit') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    p <- ggplot(df, aes(x, y, fill = Value)) +
      scale_fill_gradientn(name = leg.txt, colors = cols) +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
           axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
           axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
    if (y$grid$topo == 'hexagonal') {
      if (hover) {
        suppressWarnings(
          p <- p + geom_hex(aes(text = Probes), stat = 'identity')
        )
      } else {
        p <- p + geom_hex(stat = 'identity')
      }
    } else if (y$grid$topo == 'rectangular') {
      if (hover) {
        suppressWarnings(
          p <- p + geom_raster(aes(text = Probes))
        )
      } else {
        p <- p + geom_raster()
      }
    }
  }

  # Output
  gg_out(p, hover, legend)
  if (export) {
    return(y)
  }

}


