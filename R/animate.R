# Animate a sequence of neuron/mesh states as a ggplot2 GIF. This is the general
# "morph -> GIF" primitive: geodesic warps, developmental time series, rotations, any
# per-timepoint list of nat objects. It reuses geom_neuron() for every frame and
# assembles the frames with gifski (or magick) into an optionally ping-ponged GIF.

#' Animate neuron/mesh objects across timepoints as a ggplot2 GIF
#'
#' Renders a set of objects that change over time (each a per-timepoint list of
#' `neuron`/`neuronlist`/`mesh3d`/`hxsurf`/point objects accepted by [geom_neuron()])
#' to a GIF, each object in its own colour, over an optional fixed reference `volume`
#' and an optional set of static `targets` drawn in a transparent greyscale beneath
#' the moving objects (e.g. the fixed shapes a warp should land on).
#'
#' The simplest use is a **turntable**: pass a single spatial object as `x` (a
#' `neuron`/`neuronlist`/`mesh3d`/`hxsurf`/`dotprops`/matrix) and, with `flows` left
#' `NULL`, it is spun about its centroid (or `centre`) to a looping rotation GIF.
#'
#' @param x Either a single spatial object to turn into a turntable (when `flows` is
#'   `NULL`), or — for back-compatibility — a `flows` list passed positionally.
#' @param flows A named list; each element is the per-timepoint list of objects for
#'   one structure (all the same length, one entry per frame). Colours are per
#'   structure. If `NULL` (default) it is built from `x`: a turntable if `x` is a
#'   single spatial object, otherwise `x` is taken to already be the flows list.
#' @param centre Rotation centre for the turntable (length-3). Defaults to the
#'   `volume`'s centroid if a volume is given, else `x`'s centroid.
#' @param turntable_frames Number of frames in an auto-built turntable.
#' @param spin_axis Axis the turntable rotates about: `"y"` (default), `"x"` or `"z"`.
#' @param cols Named vector of colours (one per `flows` entry); recycled from a
#'   default palette if `NULL`.
#' @param volume Optional fixed reference object drawn under every frame (e.g. a brain
#'   surface / hull).
#' @param volume_col Colour of the reference volume (default `"grey80"`; pass e.g. a
#'   light pink to render a brain hull as a soft envelope).
#' @param targets Optional named list of fixed reference objects drawn *above* the
#'   volume but *below* the flow (e.g. the fixed targets each moving object should land
#'   on). Static across frames.
#' @param target_cols Named vector of colours for `targets`; defaults to a transparent
#'   greyscale ramp so the coloured flow reads clearly on top.
#' @param volume_alpha,target_alpha Alphas for the reference volume and target objects.
#' @param alpha Alpha for the moving (flow) objects: a single value for all, or a
#'   per-structure named vector (matched to `flows` by name) so e.g. a rotating brain
#'   mesh can be translucent while neurons stay opaque.
#' @param rotation_matrix Optional 4x4 view matrix (as [geom_neuron()]/rgl use),
#'   applied to every object so all layers share one projection.
#' @param file Output GIF path. Frames are written next to it. If `NULL`, frames are
#'   written to a temporary directory and their paths returned (no GIF assembled).
#' @param width,height,delay,dpi GIF frame size (px), per-frame delay (s) and dpi.
#' @param pingpong If `TRUE` (default) the frame sequence plays forward then back for
#'   a seamless loop.
#' @param fixed_limits If `TRUE` (default) all frames share one set of x/y limits
#'   (spanning every frame plus the volume and targets) so moving objects don't make
#'   the camera appear to zoom; set `FALSE` to let each frame auto-scale.
#' @param points Optional named list of point sets drawn as circles on top of the
#'   neurons — e.g. neuron root points / somata. Each entry is either a static N x 3
#'   matrix (same every frame) or, like `flows`, a per-timepoint list of matrices
#'   (circles that move with the warp). Anything [nat::xyzmatrix()] accepts works.
#' @param point_cols Named vector of fill colours for `points` (one per entry);
#'   recycled from a default palette if `NULL`.
#' @param point_size,point_alpha,point_stroke Circle size, alpha and outline width.
#' @return The GIF path if written (needs `gifski`, or falls back to `magick`), else
#'   the vector of frame PNG paths.
#' @seealso [geom_neuron()], [ggneuron()].
#' @examples
#' \dontrun{
#' # Turntable of a single object (spun about its centroid):
#' ggneuron_gif(banc.brain_neuropil, file = "spin.gif")
#'
#' # Turntable of neurons about the brain-hull centroid, brain shown translucent:
#' ggneuron_gif(banc.skels[1:6], volume = banc.brain_neuropil,
#'              rotation_matrix = banc_view, file = "spin.gif")
#'
#' # Explicit animation: `warp` is a named list of per-timepoint mesh states.
#' ggneuron_gif(warp, volume = brain, volume_col = "lightpink",
#'              targets = target_meshes, file = "morph.gif")
#' }
#' @export
ggneuron_gif <- function(x, flows = NULL, cols = NULL, volume = NULL, volume_col = "grey80",
                         volume_alpha = 0.12, targets = NULL, target_cols = NULL,
                         target_alpha = 0.18, alpha = 0.6, rotation_matrix = NULL,
                         file = NULL, width = 900, height = 800, delay = 0.14,
                         dpi = 96, pingpong = TRUE, fixed_limits = TRUE,
                         centre = NULL, turntable_frames = 36L, spin_axis = c("z", "y", "x"),
                         points = NULL, point_cols = NULL, point_size = 3,
                         point_alpha = 1, point_stroke = 0.8) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggneuron_gif() needs the 'ggplot2' package.", call. = FALSE)
  spin_axis <- match.arg(spin_axis)
  spatial <- function(o) inherits(o, c("neuron", "neuronlist", "mesh3d", "hxsurf", "dotprops")) || is.matrix(o)

  # Build `flows` if not given. A single spatial `x` (neuron/neuronlist/mesh/hxsurf/...)
  # becomes a TURNTABLE: it (and the `volume`, if any) is spun about `centre` (default
  # the volume's, else x's, centroid) over `turntable_frames`, seamlessly looping. A
  # flows-list passed as `x` is used as-is (back-compatible).
  if (is.null(flows)) {
    if (spatial(x)) {
      ctr <- if (!is.null(centre)) centre
             else if (!is.null(volume)) colMeans(nat::xyzmatrix(volume))
             else colMeans(nat::xyzmatrix(x))
      angs <- utils::head(seq(0, 2 * pi, length.out = turntable_frames + 1L), -1L)
      # Spin about a world axis (default Z, the usual dorsal-ventral / up axis for
      # natverse brain data) via Rodrigues' formula; the view is set separately by
      # rotation_matrix. Pick spin_axis to suit your object's up direction.
      ax <- switch(spin_axis, x = c(1, 0, 0), y = c(0, 1, 0), z = c(0, 0, 1))
      ax <- ax / sqrt(sum(ax^2))
      K  <- matrix(c(0, -ax[3], ax[2], ax[3], 0, -ax[1], -ax[2], ax[1], 0), 3, byrow = TRUE)
      Rmat <- function(th) diag(3) + sin(th) * K + (1 - cos(th)) * (K %*% K)
      spin <- function(o, th) { V <- nat::xyzmatrix(o)
        nat::xyzmatrix(o) <- sweep(sweep(V, 2, ctr) %*% t(Rmat(th)), 2, -ctr); o }
      flows <- list()
      if (!is.null(volume)) { flows[["volume"]] <- lapply(angs, function(th) spin(volume, th)); volume <- NULL }
      flows[["object"]] <- lapply(angs, function(th) spin(x, th))
      if (is.null(cols))
        cols <- if ("volume" %in% names(flows)) list(volume = volume_col, object = "navy") else list(object = "navy")
      if (length(alpha) == 1L && is.null(names(alpha)))
        alpha <- if ("volume" %in% names(flows)) c(volume = volume_alpha, object = alpha) else c(object = alpha)
      pingpong <- FALSE                       # a full rotation already loops seamlessly
    } else flows <- x                         # `x` is already a flows-list
  }

  # Allow a single object per structure (wrap as a one-frame flow).
  flows <- lapply(flows, function(f) if (is.list(f) && !inherits(f, c("neuron", "mesh3d"))) f else list(f))
  nstruct <- length(flows)
  if (is.null(cols)) {
    pal <- c("#E4572E", "#2E86AB", "#F3A712", "#3B7A57", "#8E44AD",
             "#17A398", "#C1121F", "#5B8C5A", "#E27396", "#7D4F50")
    cols <- stats::setNames(rep(pal, length.out = nstruct), names(flows))
  }
  if (!is.null(targets) && is.null(target_cols))
    target_cols <- stats::setNames(
      grDevices::grey.colors(length(targets), start = 0.4, end = 0.75), names(targets))
  # per-structure alpha: expand a scalar, or match a named/positional vector to flows
  if (length(alpha) == 1L && is.null(names(alpha)))
    alpha <- stats::setNames(rep(alpha, nstruct), names(flows))
  else if (is.null(names(alpha)))
    alpha <- stats::setNames(alpha, names(flows))

  nt <- length(flows[[1]])
  # Normalise `points` to one N x 3 matrix per frame per structure (static -> repeated).
  if (!is.null(points)) {
    points <- lapply(points, function(pp) {
      per_frame <- if (is.list(pp) && !is.matrix(pp) &&
                       !inherits(pp, c("neuron", "mesh3d", "dotprops"))) pp else rep(list(pp), nt)
      lapply(per_frame, nat::xyzmatrix)
    })
    if (is.null(point_cols))
      point_cols <- stats::setNames(rep(c("#C1121F", "#2E86AB", "#F3A712", "#3B7A57"),
                                        length.out = length(points)), names(points))
  }
  frames <- character(nt)
  fdir <- if (is.null(file)) tempdir() else dirname(file)
  dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  file.remove(list.files(fdir, "^flowframe_[0-9]+\\.png$", full.names = TRUE))  # stale frames

  # Fixed camera: compute one set of x/y limits spanning EVERY frame (plus volume and
  # targets) so moving objects don't make each frame auto-rescale (which reads as an
  # unwanted zoom). Projection matches geom_neuron: rows 1:2 of rotation_matrix[,1:3].
  lims <- NULL
  if (isTRUE(fixed_limits)) {
    proj_xy <- function(o) {
      if (is.null(o)) return(NULL)
      X <- tryCatch(nat::xyzmatrix(o), error = function(e) NULL); if (is.null(X)) return(NULL)
      if (!is.null(rotation_matrix)) X <- X %*% t(rotation_matrix[1:2, 1:3]) else X <- X[, 1:2, drop = FALSE]
      X
    }
    allobj <- c(list(volume), targets, unlist(flows, recursive = FALSE),
                if (!is.null(points)) unlist(points, recursive = FALSE))
    xy <- do.call(rbind, lapply(allobj, proj_xy))
    if (!is.null(xy) && nrow(xy)) {
      rx <- range(xy[, 1]); ry <- range(xy[, 2])
      lims <- list(x = rx + c(-1, 1) * diff(rx) * 0.03, y = ry + c(-1, 1) * diff(ry) * 0.03)
    }
  }
  for (t in seq_len(nt)) {
    p <- ggplot2::ggplot()
    if (!is.null(volume))
      p <- p + geom_neuron(volume, rotation_matrix = rotation_matrix,
                           cols = rep(volume_col, 2), alpha = volume_alpha)
    if (!is.null(targets))
      for (nm in names(targets))
        p <- p + geom_neuron(targets[[nm]], rotation_matrix = rotation_matrix,
                             cols = rep(target_cols[[nm]], 2), alpha = target_alpha)
    for (nm in names(flows))
      p <- p + geom_neuron(flows[[nm]][[t]], rotation_matrix = rotation_matrix,
                           cols = rep(cols[[nm]], 2), alpha = alpha[[nm]])
    # Root points / somata as filled circles on top.
    if (!is.null(points))
      for (nm in names(points)) {
        P <- points[[nm]][[t]]
        P <- if (!is.null(rotation_matrix)) P %*% t(rotation_matrix[1:2, 1:3]) else P[, 1:2, drop = FALSE]
        p <- p + ggplot2::geom_point(
          data = data.frame(px = P[, 1], py = P[, 2]), ggplot2::aes(x = .data$px, y = .data$py),
          shape = 21, fill = point_cols[[nm]], colour = "grey15",
          size = point_size, alpha = point_alpha, stroke = point_stroke)
      }
    p <- p + ggplot2::theme_void() +
      (if (is.null(lims)) ggplot2::coord_fixed()
       else ggplot2::coord_fixed(xlim = lims$x, ylim = lims$y)) +
      ggplot2::theme(legend.position = "none")
    frames[t] <- file.path(fdir, sprintf("flowframe_%03d.png", t))
    ggplot2::ggsave(frames[t], p, width = width / dpi, height = height / dpi,
                    dpi = dpi, bg = "white")
  }
  if (!is.null(file)) {
    seq <- if (pingpong && nt > 1) c(frames, rev(frames)[-c(1, nt)]) else frames
    # the per-frame PNGs are intermediate once a GIF is requested; remove them after.
    if (requireNamespace("gifski", quietly = TRUE)) {
      gifski::gifski(seq, file, width = width, height = height, delay = delay)
      file.remove(frames); return(file)
    }
    # gifski needs a Rust toolchain to build; fall back to magick if available.
    if (requireNamespace("magick", quietly = TRUE)) {
      valid <- c(1, 2, 4, 5, 10, 20, 25, 50, 100)     # magick fps must divide 100
      fps <- valid[which.min(abs(valid - 1 / delay))]
      anim <- magick::image_animate(magick::image_read(seq), fps = fps, optimize = TRUE)
      magick::image_write(anim, file)
      file.remove(frames); return(file)
    }
    message("ggneuron_gif(): install 'gifski' or 'magick' to assemble the GIF; ",
            "returning frame paths instead.")
  }
  frames
}
