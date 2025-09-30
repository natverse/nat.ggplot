#' Base ggplot2 Template for Neuroanatomy Plots
#'
#' @description
#' A pre-configured ggplot2 template with a minimal theme optimised for
#' neuroanatomy visualisations. This object provides a clean base with no axes,
#' grids, or extraneous elements, allowing the focus to remain on the
#' neuroanatomical structures.
#'
#' @format A ggplot2 object with:
#' \describe{
#'   \item{\code{coord_fixed()}}{Maintains aspect ratio for accurate morphology}
#'   \item{\code{theme_void()}}{Removes all standard plot elements}
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
#' gganat +
#'   nat.ggplot::geom_neuron(banc.skels, rotation_matrix = banc_view)
#'
#' gganat +
#'   nat.ggplot::geom_neuron(
#'     banc.brain_neuropil, rotation_matrix = banc_view,
#'     cols = c("grey95","grey85"), alpha = 0.3
#'   ) +
#'   nat.ggplot::geom_neuron(
#'     banc.skels, rotation_matrix = banc_view,
#'     cols = c("purple","magenta")
#'   )
#' }
#'
#' @seealso \code{\link{geom_neuron}}, \code{\link{ggneuron}}
#' @export
gganat <- .build_gganat()
