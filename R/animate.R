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
#' @param flows A named list; each element is the per-timepoint list of objects for
#'   one structure (all the same length, one entry per frame). Colours are per
#'   structure. A single object (not a list) is treated as a one-frame "flow".
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
#' @return The GIF path if written (needs `gifski`, or falls back to `magick`), else
#'   the vector of frame PNG paths.
#' @seealso [geom_neuron()], [ggneuron()].
#' @examples
#' \dontrun{
#' # `warp` is a named list of per-timepoint mesh states from some morph.
#' ggneuron_gif(warp, volume = brain, volume_col = "lightpink",
#'              targets = target_meshes, file = "morph.gif")
#' }
#' @export
ggneuron_gif <- function(flows, cols = NULL, volume = NULL, volume_col = "grey80",
                         volume_alpha = 0.12, targets = NULL, target_cols = NULL,
                         target_alpha = 0.18, alpha = 0.6, rotation_matrix = NULL,
                         file = NULL, width = 900, height = 800, delay = 0.14,
                         dpi = 96, pingpong = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggneuron_gif() needs the 'ggplot2' package.", call. = FALSE)
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
  frames <- character(nt)
  fdir <- if (is.null(file)) tempdir() else dirname(file)
  dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  file.remove(list.files(fdir, "^flowframe_[0-9]+\\.png$", full.names = TRUE))  # stale frames
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
    p <- p + ggplot2::theme_void() + ggplot2::coord_fixed() +
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
