#' @noRd
.build_gganat <- function() {
  ggplot2::ggplot() +
    ggplot2::coord_fixed() +                 # required: test expects CoordFixed
    ggplot2::theme_void() +                  # blank canvas (ensures a theme is present)
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
      plot.background   = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = "")
}
