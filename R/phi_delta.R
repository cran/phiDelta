#' @title Raw Confusion Statistics
#' @description reformarts the raw file data to c_statistics data
#' so it can be used for most of the functions in this package.
#' it can be used directly after loading data from a file like .csv
#' @param file raw data from a file, for example the output of \link{read.csv}.
#' the file must be formarted as follows:
#' The first column contains tho output of the classifier.  It should only be 1 or 0
#' The other columns represent the features.
#' The names of the columns 2.. are considered as the names of the features
#' @return dataframe, first column are the labels, 0 is a negative sample, 1 a positve
#' the other columns contain the
#' @author rothe
#' @examples
#' data("climate_data")
#' x <- c_statistics(climate_data)
#' @export
c_statistics <- function(file) {
  # check if file format is correct
  if (ncol(file) < 2) {
    return(NULL)
  }

  # get lables
  labels <- file[, 1]
  # get data
  raw_data <- file[, 2:ncol(file)]
  # normalize data
  normalize_Feature <- function(x) {
    x <- sweep(x, 2, apply(x, 2, min))
    x <- sweep(x, 2, apply(x, 2, max), "/")
    2 * x - 1
  }

  data <- data.frame(normalize_Feature(raw_data))
  # return dataframe
  return(data.frame(labels, data))
}

#############################################################################

#' @title calculate ratio
#' @description calculates the ratio between positive and negative samples
#' @param stats c_statistics
#' @return ratio
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' ratio <- calculate_ratio(x)
#' @export
calculate_ratio <- function(stats) {
  # number of 1 / number of 0
  return(sum(stats[1]) / (nrow(stats) - sum(stats[1])))
}

#############################################################################

#' @title confusion matrices
#' @description calculates the confusion matrices from the c_statistics
#' @param stats c_statistics
#' @return a matrix.  Each column represents a feature.
#' Each row describes in this order: true negative, FALSE negative, true positive, FALSE negative
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' cmat <- c_matrices(x)
#' @export
c_matrices <- function(stats) {
  # initialize a matrix with one col for every feature and TN, FN, TP and FP
  matrix <- matrix(nrow = 4, ncol = ncol(stats) - 1)

  # naming the rows and cols
  rownames(matrix) <- c("TN", "FN", "TP", "FP")
  colnames(matrix) <- colnames(stats[2:ncol(stats)])
  kron_neg <- stats$labels > 0.5
  kron_pos <- stats$labels < 0.5

  # calculating the values for every col
  for (x in 1:ncol(matrix)) {
    current_col <- stats[, x + 1]
    prod <- (current_col * (stats$labels * 2 - 1) + 1) / 2
    # TN
    matrix[1, x] <- sum(prod[kron_pos])
    # FN
    matrix[3, x] <- sum(prod[kron_neg])
    # TP
    matrix[2, x] <- sum(kron_neg) - matrix[3, x]
    # FP
    matrix[4, x] <- sum(kron_pos) - matrix[1, x]
  }
  return(matrix)
}

#############################################################################

#' @title normalized confusion matrices
#' @description normalizes the confusion matrices
#' @param c_matrices confusion matrices
#' @return a matrix.  Each column represents a feature.
#' Each row describes in this order: true negative rate, FALSE negative rate,
#' true positive rate, FALSE negative rate
#' @export
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' cmat <- c_matrices(x)
#' nmat <- n_matrices(cmat)
n_matrices <- function(c_matrices) {
  # TN + FP
  cm_neg <- c_matrices[1, ] + c_matrices[4, ]
  # FN + TP
  cm_pos <- c_matrices[2, ] + c_matrices[3, ]
  matrix1 <- c_matrices
  # naming
  rownames(matrix1) <- c("tn", "fn", "tp", "fp")

  matrix1[1, ] <- matrix1[1, ] / cm_neg
  matrix1[2, ] <- matrix1[2, ] / cm_pos
  matrix1[3, ] <- matrix1[3, ] / cm_pos
  matrix1[4, ] <- matrix1[4, ] / cm_neg
  return(matrix1)
}

#############################################################################

#' @title calculate phi
#' @description calculates phi out of specificity and sensitivity
#' depending on the ratio
#' @param spec is the specificity, the true negative rate
#' @param sens is the sensitivity, the true positive rate
#' @param ratio is the ratio of positive and negative of the data. The default is 1
#' @return phi
#' @export
#' @author rothe
#' @examples
#' calculate_phi(1,0)
#' calculate_phi(0.5,0.3)
calculate_phi <- function(spec, sens, ratio = 1) {
  # formula from Prof. Armano
  return(2 * (sens - ratio * spec + ratio - 1) / (ratio + 1))
}

#' @title calculate delta
#' @description calculates delta out of specificity and sensitivity
#' depending on the ratio
#' @param spec is the specificity, the true negative rate
#' @param sens is the sensitivity, the true positive rate
#' @param ratio is the ratio of positive and negative of the data. The default is 1
#' @return delta
#' @author rothe
#' @examples
#' calculate_delta(1,0)
#' calculate_delta(0.5,0.3)
#' @export
calculate_delta <- function(spec, sens, ratio = 1) {
  # formula from Prof. Armano
  return(2 * (ratio * spec + sens) / (ratio + 1)  - 1)
}
#############################################################################

#' @title Convertion of specificity and sensitivity to phi and delta
#' @description converts specificity and sensitivity to phi and delta
#' depending on the ratio
#' @param spec is the specificity, the true negative rate
#' @param sens is the sensitivity, the true positive rate
#' @param ratio is the ratio of positive and negative of the data. The default is 1
#' @return List with phi and delta vectors
#' @author neumann
#' @examples
#' phiDelta.convert(1,0)
#' phiDelta.convert(0.5,0.3, ratio = 0.8)
#' @export
phiDelta.convert <- function(spec, sens, ratio = 1) {
  phi <- calculate_phi(spec, sens, ratio)
  delta <- calculate_delta(spec, sens, ratio)
  return(list(phi, delta))
}
#############################################################################

#' @title phi delta matrix
#' @description calculates phi and delta directly from the stats
#' and generates a matrix with the names of the features, their phi
#' and their delta value
#' @param stats c_statistics
#' @param ratio_corrected locigal, if true phi and delta will be calculated in respect to the ratio
#' of positive and negative samples
#' @return dataframe, first column are the names of the features
#' second column the phi values
#' third column the delta values
#' @importFrom utils tail
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' phiDelta <- phiDelta_from_data(x, ratio_corrected = FALSE)
#' with_ratio <- phiDelta_from_data(x)
#' @export
phiDelta_from_data <- function(stats, ratio_corrected = TRUE) {
  c_mat <- c_matrices(stats)
  n_mat <- n_matrices(c_mat)
  features <- tail(colnames(stats),-1)
  ratio <- ifelse(ratio_corrected,calculate_ratio(stats),1)
  phi <- calculate_phi(n_mat[1, ], n_mat[3, ], ratio)
  delta <- calculate_delta(n_mat[1, ], n_mat[3, ], ratio)
  return(data.frame(features, phi, delta))
}
#############################################################################

#' @title Phi delta statistics from dataframe
#' @description calculates phi, delta and the ratio directly from the dataframe
#' with provided information and generates a list with the names of the features,
#' their phi and delta value and the ratio
#' @param data dataframe without labels
#' @param labels vector of labels
#' @param ratio_corrected locigal, if true phi and delta will be calculated in respect to the ratio
#' of positive and negative samples
#' @return dataframe, first column are the names of the features
#' second column the phi values
#' third column the delta values
#' @importFrom utils tail
#' @author rothe
#' @examples
#' x <- climate_data
#' phiDelta <- phiDelta.stats(x[,-1],x[,1], ratio_corrected = FALSE)
#' with_ratio <- phiDelta.stats(x[,-1],x[,1])
#' @export
phiDelta.stats <- function(data, labels, ratio_corrected = TRUE) {
  stats <- c_statistics(cbind(labels, data))
  c_mat <- c_matrices(stats)
  n_mat <- n_matrices(c_mat)
  features <- tail(colnames(stats),-1)
  r <-calculate_ratio(stats)
  ratio <- ifelse(ratio_corrected,r,1)
  phi <- calculate_phi(n_mat[1, ], n_mat[3, ], ratio)
  delta <- calculate_delta(n_mat[1, ], n_mat[3, ], ratio)
  return(list("Names" = features, "phi" = phi, "delta" = delta, "ratio" = r))
}
#############################################################################

#' @title Plot of phi delta diagram
#' @description Plots delta against phi within the phi delta diagram shape
#' @param phi numeric value or vector of phi
#' @param delta numeric value or vector of delta
#' @param names string with feature names
#' @param ratio numeric, is the ratio of positive and negative of the data
#' @param border the color of the border of the shape.  NA for no border
#' @param filling the color to fill the shape with
#' @param crossing logical, if the crossing should be drawn
#' @param iso_specificity logical, if isometric lines of the specificity should be drawn
#' @param iso_sensitivity logical, if isometric lines of the sensitivity should be drawn
#' @param iso_neg_predictive_value logical, if isometric lines of the negative predictive value should be drawn
#' @param iso_precision logical, if isometric lines of the precision should be drawn
#' @param iso_accuracy logical, if isometric lines of the accuracy should be drawn
#' @param highlighted numeric vector, indices of the points to higlight
#' highlighted points will be orange
#' @importFrom graphics plot polygon points
#' @author rothe
#' @examples
#' x <- climate_data
#' phiDelta <- phiDelta.stats(x[,-1],x[,1])
#' phiDelta.plot(phiDelta$phi, phiDelta$delta)
#' phiDelta.plot(phiDelta$phi, phiDelta$delta,
#'   ratio = phiDelta$ratio,
#'   border = "green",
#'   iso_neg_predictive_value = TRUE,
#'   crossing = FALSE)
#' @export
phiDelta.plot <-
  function(phi,
           delta,
           ratio = 1,
           names = NULL,
           border = "red",
           filling = "grey",
           crossing = TRUE,
           iso_specificity = FALSE,
           iso_sensitivity = FALSE,
           iso_neg_predictive_value = FALSE,
           iso_precision = FALSE,
           iso_accuracy = FALSE,
           highlighted = NULL) {
    limits <- borders(ratio)
    xlim <- c(limits[4, 1], limits[2, 1])
    ylim <- c(limits[3, 2], limits[1, 2])

    plot(
      phi,
      delta,
      main = "Phi Delta Diagram",
      type = "p",
      asp = 1,
      ylab = "Discriminant Capability",
      xlab = "Characteristic Capability",
      xlim=xlim,
      ylim=ylim
    )

    # draw the background shape
    polygon(limits[, 1], limits[, 2], col = filling, border = border)

    # draw additional lines if selected
    if (crossing)
      crossings(ratio)
    if (iso_specificity)
      iso_specificity(ratio)
    if (iso_sensitivity)
      iso_sensitivity(ratio)
    if (iso_neg_predictive_value)
      iso_negative_predictive_value(ratio)
    if (iso_precision)
      iso_precision(ratio)
    if (iso_accuracy) {
      iso_accuracy(ratio)
    }

    # draw border of the phi-delta space
    polygon(limits[, 1], limits[, 2], border = border)

    # set color of the points
    point_colors <-
      replace(rep("black", length(phi)), highlighted, "magenta")

    # insert points
    points(phi, delta, pch = 19, col = point_colors)
  }

#############################################################################

#' @title phi delta plot of raw statistic data
#' @description this will create a basic plot directly out of the statistic data (c_statistics)
#' @param stats matrix of the statistic data of the features and the classifier
#' @param names vector with feature names
#' @param ratio_corrected logical, if true the plot will concider the ratio
#' of the positive and negative data samples
#' @param ... further parameters for the diagram see \link{phiDelta.plot}
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' phiDelta_plot_from_data(x)
#' phiDelta_plot_from_data(x, ratio_corrected = FALSE, iso_spec = TRUE, iso_sens = TRUE)
#' @export
phiDelta_plot_from_data <-
  function(stats, names = NULL, ratio_corrected = TRUE, ...) {
    # load information from file
    data <- phiDelta_from_data(stats, ratio_corrected)
    ratio <- ifelse(ratio_corrected,calculate_ratio(stats),1)
    # plot the information
    phiDelta.plot(data$phi, data$delta, names = names, ratio = ratio, ...)
  }

#############################################################################

#' @title borders of the phi delta space
#' @description calculates the corners of the phi delta space
#' @param ratio is the ratio of positive and negative of the data.
#' The default is 1
#' @return a matrix. Each row represents a corner in the following order:
#' top, right, bottom, left
#' @export
#' @author rothe
#' @examples
#' borders(1.0)
#' borders(0.5)
#' borders(2)
borders <- function(ratio) {
  spec <- c(1, 0, 0, 1)
  sens <- c(1, 1, 0, 0)
  x <- (calculate_phi(spec, sens, ratio))
  y <- (calculate_delta(spec, sens, ratio))
  return(cbind(x, y))
}

###############################################################################

# additional Lines

#' @title Diagram crossings
#' @description adds crossings to the plot depending on the ratio
#' @param ratio is the ratio of positive and negative of the data
#' @param col the color of the lines.  Default is darkblue
#' @param ... further graphical parameters, see \link{par}
#' @importFrom graphics segments par
#' @author Neumann
#' @examples
#' x <- c_statistics(climate_data)
#' ratio <- calculate_ratio(x)
#' phiDelta_plot_from_data(x, crossing = FALSE)
#' crossings(ratio, col = "green")
#' @export
crossings <- function(ratio, col = "darkblue", ...) {
  limits <- borders(ratio)
  segments(limits[1:2, 1], limits[1:2, 2], limits[3:4, 1], limits[3:4, 2], col = col, ...)
}

#' @title isometric sensitivity lines
#' @description adds isometric lines for the sensitivity to the plot depending on the ratio
#' @param ratio numeric value for the ratio of positive and negative of the data
#' @param granularity numeric value between 0 and 1 for the granularity of the lines.
#' It is a value for the distance between 2 lines
#' @param col the color of the lines
#' @param lty the type of line, see \link{par}
#' @param ... further graphical parameters, see \link{par}
#' @importFrom graphics segments par
#' @export
#' @author Neumann
#' @examples
#' x <- c_statistics(climate_data)
#' ratio <- calculate_ratio(x)
#' phiDelta_plot_from_data(x)
#' iso_sensitivity(ratio, col = "green")
iso_sensitivity <-
  function(ratio = 1,
           granularity = 0.25,
           col = "blue",
           lty = "longdash",
           ...) {
    a <- top_left(ratio, granularity)
    b <- down_right(ratio, granularity)
    segments(a[1, ], a[2, ], b[1, ], b[2, ], lty = lty, col = col, ...)
  }

#' @title isometric specificity lines
#' @description adds isometric lines for the specificity to the plot depending on the ratio
#' @inheritParams iso_sensitivity
#' @export
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' ratio <- calculate_ratio(x)
#' phiDelta_plot_from_data(x)
#' iso_specificity(ratio, col = "green")
iso_specificity <-
  function(ratio = 1,
           granularity = 0.25,
           col = "blue",
           lty = "longdash",
           ...) {
    a <- top_right(ratio, granularity)
    b <- down_left(ratio, granularity)
    segments(a[1, ], a[2, ], b[1, ], b[2, ], lty = lty, col = col, ...)
  }

#' @title isometric precision lines
#' @description adds isometric lines for the precision to the plot depending on the ratio
#' @inheritParams iso_sensitivity
#' @export
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' ratio <- calculate_ratio(x)
#' phiDelta_plot_from_data(x)
#' iso_precision(ratio, col = "green")
iso_precision <-
  function(ratio = 1,
           granularity = 0.25,
           lty = "longdash",
           col = "blue",
           ...) {
    a <- top_right(ratio, granularity)
    a <- cbind(a, down_right(ratio, granularity))
    corner <- borders(ratio)[4, ]
    segments(corner[1], corner[2], a[1, ], a[2, ], lty = lty, col = col, ...)
  }

#' @title isometric negative predictive value lines
#' @description adds isometric lines for the negative predictive value to the plot depending on the ratio
#' @inheritParams iso_sensitivity
#' @export
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' ratio <- calculate_ratio(x)
#' phiDelta_plot_from_data(x)
#' iso_negative_predictive_value(ratio, col = "green")
iso_negative_predictive_value <-
  function(ratio = 1,
           granularity = 0.25,
           lty = "longdash",
           col = "blue",
           ...) {
    a <- top_left(ratio, granularity)
    a <- cbind(a, down_left(ratio, granularity))
    corner <- borders(ratio)[2, ]
    segments(corner[1], corner[2], a[1, ], a[2, ], lty = lty, col = col, ...)
  }

#' @title isometric accuracy lines
#' @description adds isometric lines for the accuracy to the plot depending on the ratio
#' @inheritParams iso_sensitivity
#' @importFrom utils head tail 
#' @importFrom graphics text 
#' @export
#' @author rothe
#' @examples
#' x <- c_statistics(climate_data)
#' ratio <- calculate_ratio(x)
#' phiDelta_plot_from_data(x)
#' iso_accuracy(ratio, col = "green")
iso_accuracy <-
  function(ratio = 1,
           granularity = 0.25,
           lty = "longdash",
           col = "blue",
           ...) {
    # last line needn't be drawn
    a <- head(seq(0, 1, granularity),-1)
    a <- append(a, tail(-a,-1))

    get_acc_points <- function(point, ratio) {
      # top corner
      corner1 <- borders(ratio)[1, ]
      # bot corner
      corner2 <- borders(ratio)[3, ]

      # values on the lines from the top corner
      y_diff <- corner1[2] - point
      x1 <- rbind(corner1[1] + y_diff, corner1[1] - y_diff)

      # values on the lines from the bot corner
      y_diff <- corner2[2] - point
      x2 <- rbind(corner2[1] - y_diff, corner2[1] + y_diff)

      # compare the values, always take those closer to the middle
      return(rbind(pmin(x1[1, ], x2[1, ]), pmax(x1[2, ], x2[2, ])))
    }

    # get the phi values of the shape
    x <- get_acc_points(a, ratio)
    segments(x[1, ], a, x[2, ], a, lty = lty, col = col, ...)
    text(max(x[1, ]) + 0.2,
         sort(a),
         labels = seq(granularity / 2, 1 - granularity / 2, granularity / 2))
  }

###############################################################################

#' @title X symmetric distance of a point
#' @description calculates the Distance from the positive anchor and the negative anchor
#' to the point and returns the smaller one.  That means, if y is positive the distance to
#' the positive anchor will be return, if it is negative,
#' the negative anchor distance will be calculated
#' @param x,y numerical, in this case phi and delta but in general the input coordinates
#' @param anchor vector (x,y) the anchor for the calculation of the distance
#' @return the smaller distance of (x,y) to eather the positive or negative anchor
#' @export
#' @examples
#' symmetric_distance(0.5,0.5,c(0,0))
symmetric_distance <- function(x,y, anchor){
  top <- sqrt((x - anchor[1]) * (x - anchor[1])
              + (y - anchor[2]) * (y - anchor[2]))
  bot <- sqrt((x - anchor[1]) * (x - anchor[1])
              + (y + anchor[2]) * (y + anchor[2]))
  return(pmin(top, bot))
}

#' @title distance to the middle of the space
#' @description calculates the euclidic distance of a phi delta tuple
#' to the middle of the phi delta space.  This could be used for a rating of the features
#' @param phi numeric value or vector of phi
#' @param delta numeric value or vector of delta
#' @param ratio is the ratio of positive and negative of the data. The default is 1
#' @return the euclidic distance of the tuple to the middle
#' @export
#' @author rothe
#' @examples
#' dist_to_middle(1,0,1)
#' dist_to_middle(0.5,0.3,1)
dist_to_middle <- function(phi, delta, ratio) {
  # the middle of the rectangle
  middle <- (borders(ratio)[1, ] + borders(ratio)[3, ]) / 2
  dist <- sqrt((phi - middle[1]) * (phi - middle[1])
              + (delta - middle[2]) * (delta - middle[2]))
  return(dist)
}

#' @title distance to top or bottom
#' @description calculates the distance of the tuple to the closer
#' corner of top and bottom of the phi delta space with ratio 1.
#' This can be used for a ranking of the features
#' @param phi numeric value or vector of phi
#' @param delta numeric value or vector of delta
#' @return distance to the top or the bottom corner
#' @export
#' @author rothe
#' @examples
#' dist_to_top(1,0)
#' dist_to_top(0.5,0.3)
dist_to_top <- function(phi, delta) {
  return(symmetric_distance(phi,delta, c(0,1)))
}

#' @title calculate entropy
#' @description calculates the entropy of a specificity and sensitivity tuple
#' considering the ratio
#' @param spec numeric, is the specificity, the true negative rate
#' @param sens numeric, is the sensitivity, the true positive rate
#' @param ratio numeric, is the ratio of positive and negative of the data
#' @return entropy of the tuple
#' @export
#' @author rothe
#' @examples
#' calculate_entropy(1,0)
#' calculate_entropy(0.5,0.6,0.7)
calculate_entropy <- function(spec, sens, ratio = 1) {

  # formular from Prof. Armano
  #calculate n and p
  n <- 1 / (ratio + 1)
  p <- 1 / (1 / ratio + 1)

  # change the behavior in case of zero input
  log2 <- function(x){
    ifelse(x == 0,0,log(x,base = 2))
  }

  # inner function to calculate H
  pos_neg_H <- function(x,y,ratio){
    result <- ifelse(x == 0,0,{
      k <- (1 - y) / x
      prec <- 1 / (1 + ratio * k)
      (log2(prec) + k * log2(k * prec / ratio) / ratio)
    })
    return(result)
  }

  negH <- pos_neg_H(spec,sens,ratio)
  posH <- pos_neg_H(sens,spec,ratio)

  # calculate entropy
  h = -n * spec * negH - p * sens * posH
  return(h)
}

calculate_entropy_from_stats <- function(stats, ratio_corrected = FALSE){
  cmat <- c_matrices(stats)
  nmat <- n_matrices(cmat)
  ratio <- ifelse(ratio_corrected,calculate_ratio(stats),1)
  return(calculate_entropy(nmat[1,],nmat[3,],ratio))
}

rank_delta <- function(delta) {
  rank(abs(delta))
}

#' @title ranking of the features
#' @description this function puts together a number of rankings of the features
#' @param stats c_statistics, the data input
#' @param ratio_corrected logical, true if ratio shoud be considerd
#' @param delta_dist, numeric, the delta value of the anchor for the geometrical ranking
#' see \link{symmetric_distance}
#' @author rothe
#' @export
rank_stats <- function(stats, ratio_corrected = FALSE, delta_dist = 1){
  # get ratio
  ratio <- ifelse(ratio_corrected,calculate_ratio(stats),1)
  # phi, delta
  result <- phiDelta_from_data(stats, ratio_corrected = ratio_corrected)
  # absolut delta
  delta_absolut <- rank_delta(result$delta)
  # geometric ranking
  geom_ranking <- rank(symmetric_distance(result$phi, result$delta, c(0,delta_dist)))
  # entropy ranking
  entropy_ranking <- rank(calculate_entropy_from_stats(stats, ratio_corrected = ratio_corrected))
  # distance to middle, not used
  distance_to_center <- rank(dist_to_middle(result$phi, result$delta, ratio))
  names <- ("")
  # create data frame
  result <- data.frame(result,delta_absolut,
                       geom_ranking,
                       entropy_ranking
                       # , distance_to_center
                       )
  return(result)
}
#############################################################################

# graphical helper functions to get points on each site of the phi-delta-space

down_left <- function(ratio, granularity) {
  a <- seq(0, 1, granularity)
  return(rbind(calculate_phi(a, 0, ratio), calculate_delta(a, 0, ratio)))
}

down_right <- function(ratio, granularity) {
  a <- seq(0, 1, granularity)
  return(rbind(calculate_phi(0, a, ratio), calculate_delta(0, a, ratio)))
}

top_left <- function(ratio, granularity) {
  a <- seq(0, 1, granularity)
  return(rbind(calculate_phi(1, a, ratio), calculate_delta(1, a, ratio)))
}

top_right <- function(ratio, granularity) {
  a <- seq(0, 1, granularity)
  return(rbind(calculate_phi(a, 1, ratio), calculate_delta(a, 1, ratio)))
}

# generates a grid of points
generate_grid <- function(granularity){
  x <- seq(0,1,granularity)
  spec <- rep(x,length(x))
  sens <- rep(x,each = length(x))
  return(cbind(spec,sens))
}

#' @title isometric entropy
#' @description draws isometric curves for the entropy by calculating the entropy
#' for all points in a grid and connecting those within a epsilon enviroment of the value
#' @param x numeric, is the offset for the points
#' @param ratio numeric, is the ratio
#' @param eps numeric, the epsilon for entropies to be selected
#' @param grid_granularity numeric between 0 and 1, defines the granularity of the grid
#' @export
#' @author Neumann
#' @importFrom graphics points
iso_entropy_curve <- function(x, ratio = 1, eps = 0.001, grid_granularity = 0.01){
  grid1 <- generate_grid(grid_granularity)
  entropy <- calculate_entropy(grid1[,1],grid1[,2],ratio)
  phi <- calculate_phi(grid1[,1],grid1[,2],ratio)
  delta <- calculate_delta(grid1[,1],grid1[,2],ratio)
  a <- which(entropy < x + eps & entropy > x - eps)
  pos <- a[delta[a] > 0]
  neg <- a[delta[a] < 0]
  points(phi[pos],delta[pos],type = "l")
  points(phi[neg],delta[neg],type = "l")
}
