# hidden
.build_gganat <- function() {
  ggplot2::ggplot() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::guides(fill = "none", colour = "none") +
    ggplot2::theme(
      legend.position   = "none",
      plot.title        = ggplot2::element_text(hjust = 0, size = 8, face = "bold", colour = "black"),
      axis.title.x      = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_blank(),
      axis.ticks.x      = ggplot2::element_blank(),
      axis.title.y      = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_blank(),
      axis.ticks.y      = ggplot2::element_blank(),
      axis.line         = ggplot2::element_blank(),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(0, 0, 0, 0),
      panel.spacing     = ggplot2::unit(0, "cm"),
      panel.border      = ggplot2::element_blank(),
      panel.background  = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill='transparent', color=NA)
    ) +
    ggplot2::labs(title = "")
}

#' Base ggplot2 Template for Neuroanatomy Plots
#'
#' @description
#' A pre-configured ggplot2 template with a minimal theme optimised for
#' neuroanatomy visualisations. This object provides a clean base with no axes,
#' grids, or extraneous elements, allowing the focus to remain on the
#' neuroanatomical structures. Plot background is transparent
#' (which is really annoying to do with an rgl plot, via with rgl.postscript)
#'
#' @format A ggplot2 object with:
#' \describe{
#'   \item{\code{coord_fixed()}}{Maintains aspect ratio for accurate morphology}
#'   \item{\code{theme_void()}}{Removes standard plot elements}
#'   \item{No legends}{Guides for fill and colour are hidden}
#'   \item{No axes}{All axis elements removed}
#'   \item{No grids}{Panel grids removed}
#'   \item{Zero margins}{Plot margins set to 0}
#'   \item{Transparent background}{Panel and plot backgrounds are transparent}
#' }
#'
#' @details
#' Use \code{gganat} as the base for neuroanatomy plots instead of \code{ggplot()}.
#' Simply add your \code{geom_neuron()} layers to this base object.
#'
#' @examples
#' \dontrun{
#' gganat + nat.ggplot::geom_neuron(banc.skels, rotation_matrix = banc_view)
#' }
#'
#' @seealso \code{\link{geom_neuron}}, \code{\link{ggneuron}}
#' @export
gganat <- .build_gganat()

# Hidden
.onLoad <- function(libname, pkgname) {
  # Re-assign canonical object at load-time to prevent any accidental redefinition
  assign("gganat", .build_gganat(), envir = asNamespace(pkgname))
}
