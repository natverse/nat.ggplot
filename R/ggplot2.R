#' Convert Neuron Objects to ggplot2-Compatible Data
#'
#' @description
#' This function converts 'neuron', 'mesh3d', or 'neuronlist' objects,
#' which represent 3D points linked by lines in space, into data frames
#' that describe paths compatible with ggplot2's geom_path, or geom_polygon
#' for mesh3d objects. Neuron objects are fundamental data structures in the
#' natverse ecosystem for representing neuronal morphology (see Jefferis et al., 2007;
#' and https://natverse.org/nat/articles/neurons-intro.html).
#'
#' @param x A 'neuron', 'neuronlist', or 'mesh3d' object to be converted.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame with columns X, Y, Z, and group, where each group
#' represents a continuous path in the neuron or a polygon in the mesh.
#'
#' @seealso
#' \code{\link{rgl_view}} for a way to obtain `rotation_matrix`
#'
#' @examples
#' \dontrun{
#' library(nat.ggplot)
#'
#' # Convert a single neuron to ggplot2-compatible format
#' neuron_data <- ggplot2_neuron_path(banc.skels[[1]])
#' head(neuron_data)
#'
#' # Plot with g.anat base
#' g.anat +
#'   geom_path(data = neuron_data,
#'             aes(x = X, y = Y, group = group))
#'
#' # Apply rotation matrix
#' neuron_rotated <- ggplot2_neuron_path(banc.skels[[1]],
#'                                        rotation_matrix = banc_view)
#'
#' # Convert neuronlist
#' neuronlist_data <- ggplot2_neuron_path(banc.skels)
#'
#' # Convert mesh3d object
#' mesh_data <- ggplot2_neuron_path(banc.brain_neuropil)
#' }
#'
#' @export
ggplot2_neuron_path <- function(x, rotation_matrix = NULL, ...) UseMethod('ggplot2_neuron_path')

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path neuron
#' @export
ggplot2_neuron_path.neuron <- function(x, rotation_matrix = NULL, ...){
  x$d <- x$d[order(x$d$Parent),]
  x$d <- x$d[order(x$d$PointNo),]
  npoints <- as.data.frame(nat::xyzmatrix(x))
  if(!is.null(rotation_matrix)){
    npoints <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(npoints))))
    npoints <- npoints[,-4]
    colnames(npoints) <- c("X","Y","Z")
  }
  ss <-nat::seglist(x)
  seglist <- ss[[1]]$SegList
  edges_df <- data.frame()
  for(s in 1:length(seglist)){
    g <- npoints[seglist[[s]],]
    g$group <- s
    edges_df <- rbind(edges_df,g)
  }
  edges_df
}

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path neuronlist
#' @export
ggplot2_neuron_path.neuronlist <- function(x, rotation_matrix = NULL, ...){
  ll <- lapply(x, ggplot2_neuron_path, rotation_matrix = rotation_matrix, ...)
  max.group <- 0
  for(i in 1:length(ll)){
    ll[[i]]$group <- ll[[i]]$group + max.group
    ll[[i]]$id <- i
    max.group <- max(ll[[i]]$group, na.rm = TRUE)
  }
  do.call(rbind, ll)
}

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path mesh3d
#' @export
ggplot2_neuron_path.mesh3d <- function(x, rotation_matrix = NULL, ...) {

  # Extract vertices
  vertices <- as.data.frame(t(x$vb[-4,]))
  if(!nrow(vertices)){
    warning("ggplot2_neuron_path.mesh3d given an invalid mesh3d object")
    return(data.frame(X=NA,Y=NA,Z=0,group=NA))
  }else if(nrow(vertices)<3){
    warning("ggplot2_neuron_path.mesh3d given an invalid mesh3d object")
    return(data.frame(X=NA,Y=NA,Z=0,group=NA))
  }
  colnames(vertices) <- c("X","Y","Z")

  # Apply rotation if specified
  if(!is.null(rotation_matrix)){
    vertices <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(vertices))))
    vertices <- vertices[,-4]
    colnames(vertices) <- c("X","Y","Z")
  }

  # Extract faces
  faces.matrix <- t(x$it)
  faces <- as.vector(t(faces.matrix))

  # Create a data frame of triangles
  triangles <- data.frame(
    X = vertices$X[faces],
    Y = vertices$Y[faces],
    Z = vertices$Z[faces],
    group = rep(1:nrow(faces.matrix), each = 3)
  )
  triangles <- triangles %>%
    dplyr::arrange(dplyr::desc(.data$Z))

  # return: triangles
  triangles
}

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path NULL
#' @export
ggplot2_neuron_path.NULL <- function(x, rotation_matrix = NULL, ...) {
  NULL
}

#' Create ggplot2 Geom Layer for Neuron Visualisation
#'
#' @description
#' This function creates a ggplot2 geom layer for visualising neuron objects.
#' It supports 'neuron', 'neuronlist', and 'mesh3d' objects. For split neurons
#' (created using flow centrality analysis from Schneider-Mizell et al., 2016),
#' it colours axonal and dendritic compartments differently.
#'
#' @param x A 'neuron', 'neuronlist', or 'mesh3d' object to be visualised.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param root Numeric, if >0 and x is or contains `neuron` objects,
#' then the root node is plotted as a dot of size `root`. If `FALSE` or `0` no root node is plotted.
#' @param size Numeric, the line width for neuron skeleton paths. Default is 0.5.
#' @param cols The colour to plot the neurons in. If \code{length(cols)==length(x)} each neuron will be coloured
#' by its index in `x` applied to `cols`.
#' @param stat The statistical transformation to use on the data for this layer.
#' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' @param na.rm If FALSE, the default, missing values are removed with a warning. If TRUE, missing values are silently removed.
#' @param show.legend logical. Should this layer be included in the legends? NA, the default, includes if any aesthetics are mapped.
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them.
#' @param threshold the minimum threshold for healing neuron skeletons for visualisation, using `nat::stitch_neurons_mst`.
#' @param ... Other arguments passed on to layer().
#'
#' @return A list of ggplot2 geom layers for visualising the neuron.
#'
#' @examples
#' \dontrun{
#' library(nat.ggplot)
#'
#' # Plot a single neuron
#' g.anat +
#'   geom_neuron(banc.skels[[1]], rotation_matrix = banc_view)
#'
#' # Plot all neurons with custom colours
#' g.anat +
#'   geom_neuron(banc.skels,
#'               rotation_matrix = banc_view,
#'               cols = c("purple", "magenta"))
#'
#' # Plot brain mesh as context
#' g.anat +
#'   geom_neuron(banc.brain_neuropil,
#'               rotation_matrix = banc_view,
#'               cols = c("grey75", "grey50"),
#'               alpha = 0.3)
#'
#' # Plot split neurons showing axon/dendrite
#' g.anat +
#'   geom_neuron(banc.neurons.flow[[1]],
#'               rotation_matrix = banc_view)
#'
#' # Plot synapses as points
#' g.anat +
#'   geom_neuron(as.matrix(banc.syns[, c("X", "Y", "Z")]),
#'               rotation_matrix = banc_view,
#'               root = 0.5,
#'               cols = c("navy", "red"))
#' }
#'
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
geom_neuron <- function(x = NULL,
                        rotation_matrix = NULL,
                        root = 3,
                        size = 0.5,
                        cols = c("navy", "turquoise"),
                        stat = "identity",
                        position = "identity",
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = FALSE,
                        threshold = Inf,
                        ...) UseMethod('geom_neuron')

#' @rdname geom_neuron
#' @method geom_neuron neuron
#' @export
geom_neuron.neuron <- function(x = NULL,
                               rotation_matrix = NULL,
                               root = 3,
                               size = 0.5,
                               cols = c("navy", "turquoise"),
                               stat = "identity",
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = FALSE,
                               threshold = Inf,
                               ...) {
  if(root){
    x$tags$soma <- nat::rootpoints(x)
  }
  # Try to get soma using catmaid if available, otherwise use rootpoints
  soma <- if (requireNamespace("catmaid", quietly = TRUE)) {
    catmaid::soma(x)
  } else {
    nat::xyzmatrix(x)[nat::rootpoints(x), , drop = FALSE]
  }
  if(is.null(soma)|ncol(soma)!=3){
    soma <- nat::xyzmatrix(x)[nat::rootpoints(x), , drop = FALSE]
  }
  soma <- t(as.data.frame(soma, ncol = 3))
  if(!is.null(rotation_matrix)){
    soma <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(soma))), ncol = 3)
    soma <- soma[,-4]
    colnames(soma) <- c("X","Y","Z")
  }
  x <- ggplot2_neuron_path.neuron(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, color = .data$Z, group = .data$group),
                       data = x,
                       linewidth = size,
                       stat = stat,
                       position = position,
                       na.rm = na.rm,
                       show.legend = show.legend,
                       inherit.aes = inherit.aes,
                       ...),
    ggplot2::geom_point(mapping = ggplot2::aes(x = .data$X, y = .data$Y),
                        data = soma,
                        color = "black", #cols[1],
                        size = root,
                        ...),
    ggplot2::scale_color_gradient(low = cols[1],
                                  high = cols[length(cols)]),
    ggnewscale::new_scale_colour()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron neuronlist
#' @export
geom_neuron.neuronlist <- function(x = NULL,
                                   rotation_matrix = NULL,
                                   root = 3,
                                   size = 0.5,
                                   cols = c("navy", "turquoise"),
                                   stat = "identity",
                                   position = "identity",
                                   na.rm = FALSE,
                                   show.legend = NA,
                                   inherit.aes = FALSE,
                                   threshold = Inf,
                                   ...) {
  glist <- list()

  # Handle empty neuronlist
  if(length(x) == 0){
    return(list(ggplot2::geom_blank()))
  }

  if(length(x)!=1){
    if(cols[1]=="rainbow"){
      cols <- grDevices::rainbow(length(x))
    }else if(length(cols)==1){
      cols <- rep(cols, length(x))
    }else if(length(cols)==length(x)){
      cols <- cols
    }else{
      cols <- grDevices::colorRampPalette(c(cols))(length(x))
    }
    for(i in 1:length(x)){
      glist[[i]] <- geom_neuron(x = x[[i]],
                                rotation_matrix = rotation_matrix,
                                cols = cols[i],
                                root = root,
                                size = size,
                                stat = stat,
                                position = position,
                                na.rm = na.rm,
                                show.legend = show.legend,
                                inherit.aes = FALSE,
                                threshold = threshold,
                                ...)
    }
  }else{
    if(cols[1]=="rainbow"){
      cols <-c("#6D2C7B", "#FF1493")
    }
    for(i in 1:length(x)){
      glist[[i]] <- geom_neuron(x = x[[i]], rotation_matrix = rotation_matrix, cols = cols,
                                root = root, size = size,
                                stat = stat, position = position, na.rm = na.rm, show.legend = show.legend,
                                inherit.aes = FALSE, ...)
    }
  }
  glist
}

#' @rdname geom_neuron
#' @method geom_neuron mesh3d
#' @export
geom_neuron.mesh3d <- function(x = NULL,
                               rotation_matrix = NULL,
                               root = 3,
                               size = 0.5,
                               cols = c("navy", "turquoise"),
                               stat = "identity",
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = FALSE,
                               threshold = Inf,
                               ...) {
  x <- ggplot2_neuron_path.mesh3d(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_polygon(data = x, mapping = ggplot2::aes(x = .data$X, y = .data$Y, fill = .data$Z, group = .data$group),
                          color = NA,
                          stat = stat, position = position, na.rm = na.rm,
                          show.legend = show.legend, inherit.aes = inherit.aes, ...),
    ggplot2::scale_fill_gradient(low = cols[1], high = cols[length(cols)]),
    ggnewscale::new_scale_fill()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron hxsurf
#' @export
geom_neuron.hxsurf <- function(x = NULL,
                               rotation_matrix = NULL,
                               root = 3,
                               size = 0.5,
                               cols = c("navy", "turquoise"),
                               stat = "identity",
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = FALSE,
                               threshold = Inf,
                               ...) {
  x <- rgl::as.mesh3d(x)
  geom_neuron.mesh3d(x=x, rotation_matrix=rotation_matrix, cols=cols,
                     size=size, stat=stat, position=position, na.rm=na.rm, show.legend=show.legend, inherit.aes=inherit.aes,
                     ...)
}

#' @rdname geom_neuron
#' @method geom_neuron NULL
#' @export
geom_neuron.NULL <- function(x = NULL,
                             rotation_matrix = NULL,
                             root = 3,
                             size = 0.5,
                             cols = c("navy", "turquoise"),
                             stat = "identity",
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = FALSE,
                             threshold = Inf,
                             ...) {
  list(
    ggplot2::geom_polygon(...)
  )
}

#' @rdname geom_neuron
#' @method geom_neuron list
#' @export
geom_neuron.list <- function(x = NULL,
                             rotation_matrix = NULL,
                             root = 3,
                             size = 0.5,
                             cols = c("navy", "turquoise"),
                             stat = "identity",
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = FALSE,
                             threshold = Inf,
                             ...) {
  if(is_named_all(cols) & is_named_all(x)){
    if(all(names(x)%in%names(cols))){
      cols <- cols[names(x)]
    }
  }
  if(length(x)){
    geom_neuron.neuronlist(x=x,
                           rotation_matrix=rotation_matrix,
                           cols=cols,
                           size=size,
                           root=root,
                           stat=stat,
                           position=position,
                           na.rm=na.rm,
                           show.legend=show.legend,
                           inherit.aes=inherit.aes,
                           ...)
  }else{
    geom_neuron.NULL(x = x, ...)
  }
}

#' @rdname geom_neuron
#' @method geom_neuron matrix
#' @export
geom_neuron.matrix <- function(x = NULL,
                               rotation_matrix = NULL,
                               root = 3,
                               size = 0.5,
                               cols = c("navy", "turquoise"),
                               stat = "identity",
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = FALSE,
                               threshold = Inf,
                               ...) {
  x<-as.data.frame(nat::xyzmatrix(x))
  if(!is.null(rotation_matrix)){
    x <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(x))))
    x <- x[,-4]
    colnames(x) <- c("X","Y","Z")
  }
  list(
    ggplot2::geom_point(data = x,
                        mapping = ggplot2::aes(x = .data$X, y = .data$Y, color = .data$Z),
                        size = size,
                        ...),
    ggplot2::scale_color_gradient(low = cols[1], high = cols[length(cols)]),
    ggnewscale::new_scale_colour()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron data.frame
#' @export
geom_neuron.data.frame <- function(x = NULL,
                                   rotation_matrix = NULL,
                                   root = 3,
                                   size = 0.5,
                                   cols = c("navy", "turquoise"),
                                   stat = "identity",
                                   position = "identity",
                                   na.rm = FALSE,
                                   show.legend = NA,
                                   inherit.aes = FALSE,
                                   threshold = Inf,
                                   ...) {
  x <- nat::xyzmatrix(x)
  geom_neuron.matrix(x,
                     rotation_matrix = rotation_matrix,
                     root = root,
                     size = size,
                     cols = cols,
                     stat = stat,
                     position = position,
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = FALSE,
                     ...)
}

#' @rdname geom_neuron
#' @method geom_neuron dotprops
#' @export
geom_neuron.dotprops <- function(x = NULL,
                                 rotation_matrix = NULL,
                                 root = 3,
                                 size = 0.5,
                                 cols = c("navy", "turquoise"),
                                 stat = "identity",
                                 position = "identity",
                                 na.rm = FALSE,
                                 show.legend = NA,
                                 inherit.aes = FALSE,
                                 threshold = Inf,
                                 ...) {
  x<-as.data.frame(nat::xyzmatrix(x))
  geom_neuron.data.frame(x,
                         rotation_matrix = rotation_matrix,
                         root = root,
                         size = size,
                         cols = cols,
                         stat = stat,
                         position = position,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         threshold = Inf,
                         ...)
}

#' @rdname geom_neuron
#' @method geom_neuron synapticneuron
#' @export
geom_neuron.synapticneuron <- function(x = NULL,
                                       rotation_matrix = NULL,
                                       root = 3,
                                       size = 0.5,
                                       cols = c("navy", "turquoise"),
                                       stat = "identity",
                                       position = "identity",
                                       na.rm = FALSE,
                                       show.legend = NA,
                                       inherit.aes = FALSE,
                                       threshold = Inf,
                                       ...) {
  geomneuron<-geom_neuron.neuron(x = x,
                                 rotation_matrix = rotation_matrix,
                                 root = root,
                                 size = size,
                                 cols = cols,
                                 stat = stat,
                                 position = position,
                                 na.rm = na.rm,
                                 show.legend = show.legend,
                                 inherit.aes = inherit.aes,
                                 ...)
  if(!is.null(x$connectors)){
    syns.in <- nat::xyzmatrix(x$connectors[x$connectors$prepost == 1,])
    syns.out <- nat::xyzmatrix(x$connectors[x$connectors$prepost == 0,])
    if(!is.null(rotation_matrix)){
      syns.in <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.in)))
      syns.in <- syns.in[,-4]
      colnames(syns.in) <- c("X","Y","Z")
      syns.out <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.out)))
      syns.out <- syns.out[,-4]
      colnames(syns.out) <- c("X","Y","Z")
    }
    glist <- list(
      ggplot2::geom_point(data = syns.in,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#132157",
                          size = root/50,
                          alpha = 0.25),
      ggplot2::geom_point(data = syns.out,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#D72000",
                          size = root/50,
                          alpha = 0.25)
    )
    c(geomneuron,glist)
  }else{
    geomneuron
  }
}

#' @rdname geom_neuron
#' @method geom_neuron splitneuron
#' @export
geom_neuron.splitneuron <- function(x = NULL,
                                    rotation_matrix = NULL,
                                    root = 3,
                                    size = 0.5,
                                    cols = c("navy", "turquoise"),
                                    stat = "identity",
                                    position = "identity",
                                    na.rm = FALSE,
                                    show.legend = NA,
                                    inherit.aes = FALSE,
                                    threshold = Inf,
                                    ...) {

  # Get parts
  if(root){
    x$tags$soma <- nat::rootpoints(x)
  }
  soma <- if (requireNamespace("catmaid", quietly = TRUE)) {
    catmaid::soma(x)
  } else {
    nat::xyzmatrix(x)[nat::rootpoints(x), , drop = FALSE]
  }
  if(is.null(soma)|ncol(soma)!=3){
    soma <- nat::xyzmatrix(x)[nat::rootpoints(x), , drop = FALSE]
  }
  soma <- t(as.data.frame(soma, ncol = 3))
  if(!is.null(rotation_matrix)){
    soma <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(soma))))
    soma <- soma[,-4]
    colnames(soma) <- c("X","Y","Z")
  }
  rownames(x$d) <- 1:nrow(x$d)
  dendrites.v = subset(rownames(x$d), x$d$Label == 3)
  axon.v = subset(rownames(x$d), x$d$Label == 2)
  p.d.v = subset(rownames(x$d), x$d$Label == 4)
  p.n.v = subset(rownames(x$d), x$d$Label == 7)
  null.v = subset(rownames(x$d), x$d$Label ==  0 | is.na(x$d$Label))

  # Get cable
  dendrites <- tryCatch(nat::prune_vertices(x,
                                            .progress = "none",
                                            verticestoprune = as.integer(setdiff(rownames(x$d),dendrites.v))),
                        error = function(e) NULL)
  axon <- tryCatch(nat::prune_vertices(x,
                                       .progress = "none",
                                       verticestoprune = as.integer(setdiff(rownames(x$d),axon.v))),
                   error = function(e) NULL)
  p.d <- tryCatch(nat::prune_vertices(x,
                                      .progress = "none",
                                      verticestoprune = as.integer(setdiff(rownames(x$d),p.d.v))),
                  error = function(e) NULL)
  p.n <- tryCatch(nat::prune_vertices(x,
                                      .progress = "none",
                                      verticestoprune = as.integer(setdiff(rownames(x$d),p.n.v))),
                  error = function(e) NULL)
  nulls <- tryCatch(nat::prune_vertices(x,
                                        .progress = "none",
                                        verticestoprune = as.integer(setdiff(rownames(x$d),null.v))),
                    error = function(e) NULL)

  # Stitch subtree
  dendrites <- tryCatch(nat::stitch_neurons_mst(dendrites, threshold=threshold, .progress = "none"), error = function(e) NULL)
  axon <- tryCatch(nat::stitch_neurons_mst(axon, threshold=threshold, .progress = "none"), error = function(e) NULL)
  p.d <- tryCatch(nat::stitch_neurons_mst(p.d, threshold=threshold), .progress = "none", error = function(e) NULL)
  p.n <- tryCatch(nat::stitch_neurons_mst(p.n, threshold=threshold, .progress = "none"), error = function(e) NULL)

  # Make into a multi-segment neuronlist
  nulls <- nat::nlapply(1:length(nulls$SubTrees), function(subt) tryCatch(nat::prune_vertices(nulls,
                                                                                              verticestoprune = unlist(nulls$SubTrees[[subt]]),
                                                                                              invert = TRUE,
                                                                                              .progress = "none"),
                                                                          error = function(e) NULL),
                        .progress = "none")
  nulls <- nulls[unlist(lapply(nulls, length))>0]
  if(!length(nulls)){
    nulls <- NULL
  }

  # Make ggplot2 objects
  g.dendrites <- ggplot2_neuron_path(dendrites, rotation_matrix = rotation_matrix)
  g.axon <- ggplot2_neuron_path(axon, rotation_matrix = rotation_matrix)
  g.p.d <- ggplot2_neuron_path(p.d, rotation_matrix = rotation_matrix)
  g.p.n <- ggplot2_neuron_path(p.n, rotation_matrix = rotation_matrix)
  g.nulls <- ggplot2_neuron_path(nulls, rotation_matrix = rotation_matrix)

  # Make geom objects
  glist <- list(
    if(length(g.dendrites)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.dendrites, col = "#54BCD1", na.rm = TRUE,
                         linewidth = size,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes, alpha = 1)
    }else{
      NULL
    },
    if(length(g.axon)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.axon, col = "#EF7C12", na.rm = TRUE,
                         linewidth = size,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    }else{
      NULL
    },
    if(length(g.p.d)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.p.d, col = "#8FDA04", na.rm = TRUE,
                         linewidth = size,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    }else{
      NULL
    },
    if(length(g.p.n)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.p.n, col = "#C70E7B", na.rm = TRUE,
                         linewidth = size,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    }else{
      NULL
    },
    # if(length(g.nulls)){
    #   ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
    #                      data = g.nulls, col = "#B3B3B3", na.rm = TRUE,
    #                      linewidth = size,
    #                      stat = stat, position = position,
    #                      show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    # }
    # else{
    #   NULL
    # },
    ggplot2::geom_point(mapping = ggplot2::aes(x = .data$X, y = .data$Y),
                        data = soma,
                        color = "black", #cols[1],
                        alpha = 0.75,
                        size = root)
  )

  # And synapses?
  if(!is.null(x$connectors)){
    syns.in <- nat::xyzmatrix(x$connectors[x$connectors$prepost == 1,])
    syns.out <- nat::xyzmatrix(x$connectors[x$connectors$prepost == 0,])
    if(!is.null(rotation_matrix)){
      syns.in <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.in)))
      syns.in <- syns.in[,-4]
      colnames(syns.in) <- c("X","Y","Z")
      syns.out <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.out)))
      syns.out <- syns.out[,-4]
      colnames(syns.out) <- c("X","Y","Z")
    }
    syn.glist <- list(
      ggplot2::geom_point(data = syns.in,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#132157",
                          size = root/100,
                          alpha = 0.25),
      ggplot2::geom_point(data = syns.out,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#D72000",
                          size = root/100,
                          alpha = 0.25)
    )
    c(glist,syn.glist)
  }else{
    glist
  }
}

#' Create a ggplot2 Visualisation of Neuron Objects
#'
#' @description
#' This function creates a complete ggplot2 visualisation for neuron objects,
#' including 'neuron', 'neuronlist', 'mesh3d', and 'hxsurf' objects. It sets up
#' a minimal theme and applies consistent styling to the plot.
#'
#' @param x A 'neuron', 'neuronlist', 'mesh3d', or 'hxsurf' object to be visualised.
#' @param volume a brain/neuropil volume to be plotted in grey, for context.
#' Defaults to NULL, no volume plotted.
#' @param info Optional. A string to be used as the plot title.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param cols1 Colour for the lowest Z values. Default is "turquoise".
#' @param cols2 Colour for the highest Z values. Default is "navy".
#' @param alpha Transparency of the neuron visualisation. Default is 0.5.
#' @param title.col Colour of the plot title. Default is "darkgrey".
#' @param ... Additional arguments passed to geom_neuron().
#'
#' @return A ggplot2 object representing the neuron visualisation.
#'
#' @details
#' This function wraps around geom_neuron() to create a complete plot with a
#' consistent, minimal theme. It removes axes, legends, and other extraneous
#' elements to focus on the neuron visualisation itself.
#'
#' @examples
#' \dontrun{
#' library(nat.ggplot)
#'
#' # Visualise neurons with brain volume as context
#' ggneuron(banc.skels,
#'          volume = banc.brain_neuropil,
#'          rotation_matrix = banc_view)
#'
#' # Visualise the brain neuropil alone
#' ggneuron(banc.brain_neuropil,
#'          rotation_matrix = banc_view,
#'          cols1 = c("lightblue", "darkblue"))
#'
#' # Visualise split neurons with custom colours
#' ggneuron(banc.neurons.flow,
#'          volume = banc.brain_neuropil,
#'          rotation_matrix = banc_view,
#'          info = "LHPD2a1 neurons with axon/dendrite split")
#'
#' # Visualise neuron meshes
#' ggneuron(banc.meshes[[1]],
#'          rotation_matrix = banc_view,
#'          cols1 = c("purple", "magenta"),
#'          alpha = 0.8)
#' }
#'
#' @references
#' Jefferis, G. S. X. E., Potter, C. J., Chan, A. M., Marin, E. C., Rohlfing, T.,
#' Maurer, C. R., & Luo, L. (2007). Comprehensive maps of Drosophila higher
#' olfactory centers: Spatially segregated fruit and pheromone representation.
#' \emph{Cell}, 128(6), 1187-1203. \doi{10.1016/j.cell.2007.01.040}
#'
#' Schneider-Mizell, C. M., Gerhard, S., Longair, M., Kazimiers, T., Li, F.,
#' Zwart, M. F., Champion, A., Midgley, F. M., Fetter, R. D., Saalfeld, S.,
#' & Cardona, A. (2016). Quantitative neuroanatomy for connectomics in Drosophila.
#' \emph{eLife}, 5, e12059. \doi{10.7554/eLife.12059}
#'
#' @seealso
#' \code{\link{geom_neuron}} for the underlying geom used in this function.
#'
#' @export
ggneuron <- function(x,
                     volume = NULL,
                     info = NULL,
                     rotation_matrix = NULL,
                     cols1 = c("turquoise","navy"),
                     cols2 =  c("grey75", "grey50"),
                     alpha = 0.5,
                     title.col = "darkgrey",
                     ...){
  g.anat +
    {if(!is.null(volume)){
      geom_neuron(x = volume, rotation_matrix = rotation_matrix, alpha = max(alpha-0.25,0.01), cols = cols2)
    }} +
    geom_neuron(x = x, rotation_matrix = rotation_matrix, cols = cols1, alpha = alpha, ...) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, size = 8, face = "bold", colour = title.col)) +
    ggplot2::labs(title = info)
}

# hidden
prune_vertices.synapticneuron <- function (x, verticestoprune, invert = FALSE, ...){
  if(length(verticestoprune)==nrow(x$d)){
    warning('no points left after pruning')
    return(NULL)
  }
  soma <- if (requireNamespace("catmaid", quietly = TRUE)) {
    catmaid::somaid(x)
  } else {
    nat::rootpoints(x)[1]
  }
  if(!is.null(soma)&&!is.na(soma)){
    if (requireNamespace("catmaid", quietly = TRUE)) {
      x$d[catmaid::somaindex(x),"Label"] <- 1
    } else {
      x$d[soma,"Label"] <- 1
    }
  }
  pruned <- nat::prune_vertices(x, verticestoprune, invert = invert, ...)
  root <- nat::rootpoints(x)
  if(!is.null(root)){
    root <- nat::xyzmatrix(x)[root,]
    pruned <- nat::reroot(x = pruned, point = c(root))
  }
  pruned$connectors <- x$connectors[x$connectors$treenode_id %in%
                                      pruned$d$PointNo, ]
  relevant.points <- subset(x$d, x$d$PointNo %in% pruned$d$PointNo)
  y <- pruned
  y$d <- relevant.points[match(pruned$d$PointNo, relevant.points$PointNo),]
  y$d$Parent <-  pruned$d$Parent
  y$tags <- lapply(x$tags, function(t) t[t %in% pruned$d$PointNo])
  y$url <- x$url
  y$headers <- x$headers
  y$AD.segregation.index = x$AD.segregation.index
  smid <- which(y$d$Label==1)[1]
  if(length(smid)){
    y$tags$soma <- y$d$PointNo[smid]
  }
  class(y) <- class(x)
  y
}
