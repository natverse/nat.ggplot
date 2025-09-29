#' Check if a Package is Available and Prompt Installation
#'
#' @description
#' Internal function to check if a required package is available.
#' If not, it provides a helpful error message with installation instructions.
#'
#' @param pkg Character string naming the package to check
#' @param reason Optional character string explaining why the package is needed
#'
#' @return Invisible TRUE if package is available, otherwise throws an error
#'
#' @keywords internal
check_package_available <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required but not installed.")
    if (!is.null(reason)) {
      msg <- paste0(msg, " It is needed for ", reason, ".")
    }

    # Provide installation instructions based on the package
    if (pkg == "ggnewscale") {
      msg <- paste0(msg, "\n\nPlease install it with:\n",
                   "  remotes::install_github('eliocamp/ggnewscale')")
    } else if (pkg == "catmaid") {
      msg <- paste0(msg, "\n\nPlease install it with:\n",
                   "  remotes::install_github('natverse/rcatmaid')")
    } else {
      msg <- paste0(msg, "\n\nPlease install it with:\n",
                   "  install.packages('", pkg, "')")
    }

    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

#' Get Current RGL View Parameters for Use as Rotation Matrix
#'
#' @description
#' Captures the current viewing parameters from an active rgl window and returns
#' them in a format suitable for use as a rotation matrix with nat.ggplot functions.
#' This allows users to interactively find their desired viewing angle using
#' \code{plot3d()} and then apply the same view to ggplot2 visualisations.
#'
#' @return A list containing:
#' \describe{
#'   \item{userMatrix}{A 4x4 rotation matrix that can be used with the
#'                     \code{rotation_matrix} parameter in nat.ggplot functions}
#'   \item{zoom}{The current zoom level}
#'   \item{windowRect}{The window dimensions}
#' }
#'
#' @details
#' The typical workflow is:
#' 1. Plot neurons or brain meshes using \code{plot3d()} from the nat package
#' 2. Interactively rotate the view using the mouse to find the desired angle
#' 3. Call \code{rgl_view()} to capture the current view parameters
#' 4. Use the \code{$userMatrix} component as the \code{rotation_matrix} argument
#'    in nat.ggplot functions like \code{geom_neuron()} or \code{ggneuron()}
#'
#' @examples
#' \dontrun{
#' library(nat)
#' library(nat.ggplot)
#' library(ggplot2)
#'
#' # First, plot neurons in 3D and rotate to desired view
#' plot3d(banc.brain_neuropil, alpha = 0.3)
#' plot3d(banc.skels)
#'
#' # Capture the view after rotating with mouse
#' my_view <- rgl_view()
#'
#' # Use the captured view in ggplot2
#' ggplot() +
#'   geom_neuron(banc.brain_neuropil,
#'               rotation_matrix = my_view$userMatrix,
#'               cols = c("grey75", "grey50"),
#'               alpha = 0.3) +
#'   geom_neuron(banc.skels,
#'               rotation_matrix = my_view$userMatrix)
#' }
#'
#' @seealso \code{\link{geom_neuron}}, \code{\link{ggneuron}}
#'
#' @export
rgl_view <- function () {
  if (!requireNamespace("rgl", quietly = TRUE))
    stop("Please install rgl to use rgl_view()", call. = FALSE)
  dput(list(
    userMatrix = rgl::par3d()$userMatrix,
    zoom       = rgl::par3d()$zoom,
    windowRect = rgl::par3d()$windowRect
  ))
}

# hidden
is_named_all <- function(x, require_unique = FALSE) {
  nm <- names(x)
  ok <- !is.null(nm) &&
    length(nm) == length(x) &&
    all(!is.na(nm)) &&
    all(nzchar(nm))
  if (require_unique) {
    ok <- ok && length(unique(nm)) == length(nm)
  }
  ok
}

