<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.github.io)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](https://natverse.github.io/nat.ggplot/reference/)
[![R-CMD-check](https://github.com/natverse/nat.ggplot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/natverse/nat.ggplot/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

[Example circuit](https://github.com/natverse/nat.ggplot/blob/main/inst/images/banner.png?raw=true)

# nat.ggplot

`nat.ggplot` enables the [Neuroanatomy Toolbox (nat)](https://jefferis.github.io/nat/) suite to create publication-quality 2D visualisations of neurons and brain meshes using [ggplot2](https://ggplot2.tidyverse.org/). This package bridges the gap between nat's 3D neuroanatomy capabilities and ggplot2's 2D plotting framework, allowing users to create figures of neuronal morphology data. This README uses sample data from the BANC connectome, but the package is valid for any neuroanatomical data from any organism, that can be represented as a neuron object, .swc, mesh3d, .surf, .obj file, etc. 

See [nat](https://natverse.org/) for details on reading neuroanatomy data into R.

## About Neuron Objects

In the [natverse ecosystem](https://natverse.org/), neurons are represented as tree structures containing 3D coordinates and connectivity information. The main data structures include:

- **neuron**: A single neuron with XYZ coordinates and parent-child connectivity defining the tree structure, `neuron$d` provides data in the style of a .swc file, where `neuron$d$Label` can give the arbour type, e.g. dendrite versus axon.
- **neuronlist**: A collection of neurons that can be manipulated together, bundled as a list of class 'neuronlist'
- **mesh3d**: 3D surface meshes representing brain regions or neuron surfaces
- **Synaptic information**: Pre- and postsynaptic site locations that can be overlaid on morphology

For more details, see the [natverse neurons introduction](https://natverse.org/nat/articles/neurons-intro.html)

[figure](https://github.com/natverse/nat.ggplot/blob/main/inst/images/figure.png?raw=true)

## Installation

You can install the development version of nat.ggplot from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github('natverse/nat.ggplot')
```

## Quick Start

``` r
library(nat.ggplot)

# Set output directory for saving figures (change as needed)
output_dir <- "inst/images/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Visualise neurons with default settings
p <- ggneuron(banc.meshes, rotation_matrix = banc_view)
p
ggsave(file.path(output_dir, "quickstart_basic.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
```

![Basic visualisation](https://github.com/natverse/nat.ggplot/blob/main/inst/images/quickstart_basic.png?raw=true)

Here is a more complicated plot, with a brain mesh, neuron meshes, their skeletons and their root points highlighted with a small black circle:

``` r
# Smooth neurons for better presentation
banc.skels_smoothed <- nat::nlapply(banc.skels,nat::smooth_neuron,sigma = 5000)

# Customise visualisation with g.anat base
p <- g.anat + 
  geom_neuron(banc.brain_neuropil, 
              rotation_matrix = banc_view,
              cols = c("grey90", "grey60"),
              alpha = 0.3) +
  geom_neuron(banc.meshes,
              rotation_matrix = banc_view,
              cols = c("#6D2C7B", "#FF1493", "green", "blue"),
              alpha = 0.3) +
  geom_neuron(banc.skels_smoothed,
              size = 0.01,
              alpa = 0.5,
              root = 1,
              rotation_matrix = banc_view,
              cols = "black")
p
ggsave(file.path(output_dir, "quickstart_custom.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
```

![Customised visualisation](https://github.com/natverse/nat.ggplot/blob/main/inst/images/quickstart_custom.png?raw=true)

### Finding the Right View

Compared with a simple `rgl` plot for neurons (e.g. `plot3d(banc.skels)`), the challenge with `nat.ggplot2` is that you need to transform your data, for the view you wish, so that the X-axis is the horziontal extent of your view and the Y-axis is the vertical.

To do this, you can call `plot3d(banc.skels)`, position neurons into the view you like, then use the helper function `rgl_view`, to extract the view as a matrix:

``` r
library(nat)  # for plot3d
library(rgl)  # for 3D interaction

# Step 1: Plot neurons in 3D and rotate to desired view
plot3d(banc.brain_neuropil, alpha = 0.3, col = "grey")
plot3d(banc.skels_smoothed)

# Step 2: Rotate with mouse to find the desired angle
# Step 3: Capture the view
my_view <- rgl_view()

# Step 4: Use the captured view in ggplot2
ggneuron(banc.skels_smoothed, rotation_matrix = my_view$userMatrix)
```

## Example Gallery

### Setting a view for neuron data

Here are some view of the BANC brain neuropil mesh, extracted in that way:

``` r
# Define the different view matrices from bancr
views <- list(
  side = structure(c(0.000764884985983372, 0.0153511334210634, 
    -0.99988180398941, 0, -0.940421104431152, -0.339961022138596, 
    -0.00593886896967888, 0, -0.340011894702911, 0.94031423330307, 
    0.0141764245927334, 0, -401395.405539944, -128785.809090088, 
    -5607.3408203126, 1), dim = c(4L, 4L)),
    front = structure(c(0.99931389093399, 0.0139970388263464, 
      -0.0342894680798054, 0, -0.0321401171386242, -0.132316529750824, 
      -0.990686297416687, 0, -0.0184037387371063, 0.991108655929565, 
      -0.131775915622711, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
      dorsal = structure(c(0.840221107006073, 0.537661552429199, 
  -0.0703445374965668, 0, -0.541468024253845, 0.838866174221039, 
  -0.0558210015296936, 0, 0.0289968233555555, 0.0849913582205772, 
  0.9959596991539, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
  ventral = structure(c(0.945645987987518, -0.325197845697403, 
  3.18996608257294e-05, 0, -0.30071958899498, -0.874427616596222, 
  0.380715191364288, 0, -0.123779848217964, -0.360031485557556, 
  -0.924692392349243, 0, 0, 0, 0, 1), dim = c(4L, 4L))
)

# Create plots for each view
plots <- list()
for(view_name in names(views)) {
  plots[[view_name]] <- ggneuron(
    banc.meshes,
    volume = NULL, # banc.brain_neuropil,
    rotation_matrix = views[[view_name]],
    cols1 = c("#8B1A89", "#FF69B4"),
    cols2 = c("grey95", "grey85"),
    alpha = 0.8,
    info = paste(substring(view_name, 1, 1), 
                 substring(view_name, 2), " view", sep="")
  )
}

# Display and save individual views
plots$side
ggsave(file.path(output_dir, "view_side.png"), plots$side, width = 5, height = 5, dpi = 300, bg = "white")

plots$front
ggsave(file.path(output_dir, "view_front.png"), plots$front, width = 5, height = 5, dpi = 300, bg = "white")

# Create a collage using cowplot
if (requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  # Combine plots in a 2x2 grid
  collage <- plot_grid(
    plots$front + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "none"),
    plots$side + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "none"),
    plots$dorsal + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "none"),
    plots$ventral + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "none"),
    ncol = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 12
  )
  collage
  ggsave(file.path(output_dir, "views_collage.png"), collage, width = 10, height = 10, dpi = 300, bg = "white")
}
```

![Multiple views collage](https://github.com/natverse/nat.ggplot/blob/main/inst/images/views_collage.png?raw=true)

### Visualising Brain Neuropil

You can visualise the brain neuropil mesh with different colours and transparency:

``` r
# Brain neuropil with custom colours
p <- ggneuron(banc.brain_neuropil,
              rotation_matrix = banc_view,
              cols1 = c("orange", "#EFC7E6"), 
              alpha = 0.1)
p
ggsave(file.path(output_dir, "brain_neuropil.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
```

![Brain neuropil](https://github.com/natverse/nat.ggplot/blob/main/inst/images/brain_neuropil.png?raw=true)

### Plotting Neurons with Custom Colours

Visualise multiple neurons with different colouring schemes:

``` r
# All neurons in one colour
p <- g.anat +
  geom_neuron(banc.skels_smoothed,
              rotation_matrix = banc_view,
              cols = c("purple"))
p
ggsave(file.path(output_dir, "neurons_single_colour.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
```

![Single colour](https://github.com/natverse/nat.ggplot/blob/main/inst/images/neurons_single_colour.png?raw=true)

``` r
# Colour gradient in the Z dimension:
p <- g.anat +
  geom_neuron(banc.meshes[1],
              rotation_matrix = banc_view,
              cols = c("purple", "orange"))
p
ggsave(file.path(output_dir, "neurons_gradient.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
```

![Gradient colouring](https://github.com/natverse/nat.ggplot/blob/main/inst/images/neurons_gradient.png?raw=true)

``` r
# Each neuron in a different colour
neuron_colours <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8")
p <- g.anat +
  geom_neuron(banc.brain_neuropil,
              rotation_matrix = banc_view,
              cols = c("grey95", "grey80"),
              alpha = 0.2) +
  geom_neuron(banc.skels_smoothed,
              rotation_matrix = banc_view,
              cols = neuron_colours)
p
ggsave(file.path(output_dir, "neurons_multi_colour.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
```

![Multiple colours](https://github.com/natverse/nat.ggplot/blob/main/inst/images/neurons_multi_colour.png?raw=true)

### Visualising Synapses

Display synaptic sites as coloured points:

``` r
# Create synapse visualisation
p <- g.anat +
  # Neuron meshes
  geom_neuron(banc.meshes,
              rotation_matrix = banc_view,
              cols = c("grey60", "grey40"),
              alpha = 0.8) +
  # Synapses as points
  geom_neuron(banc.syns %>% # Just output synapses
                       dplyr::filter(prepost==0),
             rotation_matrix = banc_view,
             root = 0.5,  # Point size for synapses
             cols = c("#EE4244", "#8B0000"),
             alpha = 0.6) +
  geom_neuron(banc.syns %>% # Just input synapses
                       dplyr::filter(prepost==1),
             rotation_matrix = banc_view,
             root = 0.5,  # Point size for synapses
             cols = c("#1BB6AF", "#121B56"),
             alpha = 0.6)
p
ggsave(file.path(output_dir, "synapses_visualisation.png"), p, width = 7, height = 7, dpi = 300, bg = "white")
```

![Synapse visualisation](https://github.com/natverse/nat.ggplot/blob/main/inst/images/synapses_visualisation.png?raw=true)

### Axon-Dendrite Split Visualisation

Visualise neurons with compartments coloured by their functional identity (using flow centrality from Schneider-Mizell et al., 2016). Note that this algorithm in implemented in the R package `hemibrainr`, [here](https://natverse.org/hemibrainr/reference/flow_centrality.html).

``` r
# Single neuron with axon/dendrite split
p <- ggneuron(banc.neurons.flow[[1]], 
              threshold = 20000, # Max distance over which a line is drawn between matching arbour
              rotation_matrix = banc_view,
              info = paste("neuron: ", names(banc.neurons.flow)[1]))
p
ggsave(file.path(output_dir, "split_single_neuron.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
```

![Single split neuron](https://github.com/natverse/nat.ggplot/blob/main/inst/images/split_single_neuron.png?raw=true)

``` r
# All split neurons with brain context
p <- ggneuron(banc.neurons.flow,
              #volume = banc.brain_neuropil,
              threshold = 20000,
              rotation_matrix = banc_view,
              info = "LHPD2a1 neurons:\n axon (orange), dendrite (cyan), primary neurite (purple), linker (green)\n inputs (navy), outputs (red)")
p
ggsave(file.path(output_dir, "split_all_neurons.png"), p, width = 7, height = 7, dpi = 300, bg = "white")
```

![All split neurons](https://github.com/natverse/nat.ggplot/blob/main/inst/images/split_all_neurons.png?raw=true)

### Creating Figures

Combine multiple elements for a figure:

``` r
library(ggplot2)

# Create comprehensive visualisation with zoom
# Calculate bounding box of meshes after rotation
mesh_bounds <- do.call(rbind, lapply(banc.meshes, function(mesh) {
  vertices <- t(mesh$vb[-4,])
  if(!is.null(banc_view)) {
    vertices <- t(banc_view[1:3, 1:3] %*% t(vertices))
  }
  data.frame(X = vertices[,1], Y = vertices[,2])
}))

# Add some padding (10% on each side)
x_range <- range(mesh_bounds$X)
y_range <- range(mesh_bounds$Y)
x_padding <- diff(x_range) * 0.1
y_padding <- diff(y_range) * 0.1

p <- g.anat +
  # Add brain outline
  geom_neuron(banc.brain_neuropil,
              rotation_matrix = banc_view,
              cols = c("grey95", "grey85"),
              alpha = 0.3) +
  # Add neuron meshes
  geom_neuron(banc.meshes,
              rotation_matrix = banc_view,
              cols = c("grey60", "grey40"),
              alpha = 0.8) +
  # Add split neurons
  geom_neuron(banc.neurons.flow,
              threshold = 15000,
              root = 2,
              size = 0.1,
              rotation_matrix = banc_view) +
  # Zoom to neuron mesh bounds
  coord_fixed(xlim = c(x_range[1] - x_padding, x_range[2] + x_padding),
              ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
  # Add title
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(title = "LHPD2a1",
       subtitle = "lateral horn output neurons from the BANC connectome")
p

# Save figure
# ggsave("lhpd2a1_neurons.pdf", p, width = 12, height = 12, dpi = 300)
ggsave(file.path(output_dir, "figure_comprehensive.png"), p, width = 7, height = 7, dpi = 300, bg = "white")
```

![Comprehensive figure](https://github.com/natverse/nat.ggplot/blob/main/inst/images/figure_comprehensive.png?raw=true)

## Data Attribution

This package includes example data from the BANC (Brain And Nerve Cord) connectome:

- **LHPD2a1 neurons**: First characterized by Dolan et al. (2018) in *Neuron*
- **BANC connectome data**: From Bates et al. (2025) in *bioRxiv*

## Advanced Features

### Integration with ggplot2 Ecosystem

The package works with ggplot2 extensions:

``` r
library(ggplot2)
library(cowplot)  # For combining plots

# Create multiple panels - save individually
p1 <- g.anat + 
  geom_neuron(banc.skels_smoothed[[1]], rotation_matrix = banc_view) +
  ggtitle("Neuron 1")
ggsave(file.path(output_dir, "panel_neuron1.png"), p1, width = 4, height = 4, dpi = 300, bg = "white")

p2 <- g.anat + 
  geom_neuron(banc.skels_smoothed[[2]], rotation_matrix = banc_view) +
  ggtitle("Neuron 2")
ggsave(file.path(output_dir, "panel_neuron2.png"), p2, width = 4, height = 4, dpi = 300, bg = "white")

# Alternative 1: Use cowplot for combining plots
if (requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  
  # Remove the theme elements that might conflict
  p1_clean <- p1 + theme_void() + ggtitle("Neuron 1")
  p2_clean <- p2 + theme_void() + ggtitle("Neuron 2")
  
  # Combine with cowplot
  p_combined <- plot_grid(p1_clean, p2_clean, ncol = 2, labels = c("A", "B"))
  p_combined
  ggsave(file.path(output_dir, "panels_combined_cowplot.png"), p_combined, width = 8, height = 4, dpi = 300, bg = "white")
}

# Alternative 2: Use gridExtra for combining plots
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  
  # Convert plots to grobs
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  
  # Arrange in grid
  p_combined <- grid.arrange(g1, g2, ncol = 2)
  ggsave(file.path(output_dir, "panels_combined_gridextra.png"), p_combined, width = 8, height = 4, dpi = 300, bg = "white")
}

# Alternative 3: Use base R graphics with png files
# Save individual plots first (already done above)
# Then combine the PNG files using magick or other image processing
if (requireNamespace("magick", quietly = TRUE)) {
  library(magick)
  
  # Read the individual plot images
  img1 <- image_read(file.path(output_dir, "panel_neuron1.png"))
  img2 <- image_read(file.path(output_dir, "panel_neuron2.png"))
  
  # Combine side by side
  img_combined <- image_append(c(img1, img2))
  
  # Save combined image
  image_write(img_combined, file.path(output_dir, "panels_combined_cowplot.png"))
}
```

![Combined panels](https://github.com/natverse/nat.ggplot/blob/main/inst/images/panels_combined.png?raw=true)

### Controlling Line Width

The `size` parameter controls the thickness of neuron skeleton lines:

``` r
# Create plots with different line widths
p_thin <- g.anat + 
  geom_neuron(banc.skels_smoothed, 
              rotation_matrix = banc_view,
              size = 0.2,  # Thin lines
              cols = c("#FF1493", "#8B008B")) +
  ggtitle("Thin lines (size = 0.2)")
ggsave(file.path(output_dir, "size_thin.png"), p_thin, width = 5, height = 5, dpi = 300, bg = "white")

p_normal <- g.anat + 
  geom_neuron(banc.skels_smoothed, 
              rotation_matrix = banc_view,
              size = 0.5,  # Default
              cols = c("#FF1493", "#8B008B")) +
  ggtitle("Normal lines (size = 0.5)")
ggsave(file.path(output_dir, "size_normal.png"), p_normal, width = 5, height = 5, dpi = 300, bg = "white")

p_thick <- g.anat + 
  geom_neuron(banc.skels_smoothed, 
              rotation_matrix = banc_view,
              size = 1,  # Thick lines
              cols = c("#FF1493", "#8B008B")) +
  ggtitle("Thick lines (size = 1)")
ggsave(file.path(output_dir, "size_thick.png"), p_thick, width = 5, height = 5, dpi = 300, bg = "white")

# Combine using cowplot
if (requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  size_comparison <- plot_grid(
    p_thin + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10, legend.position = "none")),
    p_normal + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10, legend.position = "none")),
    p_thick + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 10, legend.position = "none")),
    ncol = 3,
    labels = c("A", "B", "C"),
    label_size = 10
  ) +
  ggplot2::theme(legend.position = "none")
  size_comparison
  ggsave(file.path(output_dir, "size_comparison.png"), size_comparison, width = 12, height = 4, dpi = 300, bg = "white")
}
```

![Line width comparison](https://github.com/natverse/nat.ggplot/blob/main/inst/images/size_comparison.png?raw=true)

## Acknowledgements

This package was developed by Alexander Shakeel Bates while in the laboratory of Rachel I. Wilson at Harvard Medical School. The package leverages the [natverse](https://natverse.github.io) ecosystem for neuroanatomy in R.

## References

Key papers referenced in this package:

- Bates, A. S., Phelps, J. S., Kim, M., Yang, H. H., Matsliah, A., Ajabi, Z., Perlman, E., et al. (2025). Distributed control circuits across a brain-and-cord connectome. *bioRxiv*. doi: [10.1101/2025.07.31.667571](https://doi.org/10.1101/2025.07.31.667571)

- Bates, Alexander Shakeel, James D. Manton, Sridhar R. Jagannathan, Marta Costa, Philipp Schlegel, Torsten Rohlfing, and Gregory Sxe Jefferis. 2020. “The Natverse, a Versatile Toolbox for Combining and Analysing Neuroanatomical Data.” eLife 9 (April). [eLife.53350](https://doi.org/10.7554/eLife.53350)

- Dolan, M. J., et al. (2018). Communication from learned to innate olfactory processing centers is required for memory retrieval in Drosophila. *Neuron*, 100(3), 651-668. doi: [10.1016/j.neuron.2018.08.037](https://doi.org/10.1016/j.neuron.2018.08.037)

- Jefferis, G. S. X. E., et al. (2007). Comprehensive maps of Drosophila higher olfactory centers. *Cell*, 128(6), 1187-1203. doi: [10.1016/j.cell.2007.01.040](https://doi.org/10.1016/j.cell.2007.01.040)

- Schneider-Mizell, C. M., et al. (2016). Quantitative neuroanatomy for connectomics in Drosophila. *eLife*, 5, e12059. doi: [10.7554/eLife.12059](https://doi.org/10.7554/eLife.12059)

## Getting Help

- For issues and feature requests, please use the [GitHub issue tracker](https://github.com/natverse/nat.ggplot/issues)
- For general natverse questions, visit the [natverse documentation](https://natverse.github.io)
- For BANC connectome data, see the [bancr package](https://github.com/natverse/bancr)
