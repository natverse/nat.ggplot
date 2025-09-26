---
output:
  pdf_document: default
  html_document: default
---
<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.github.io)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](https://natverse.github.io/nat.ggplot/reference/)
[![R-CMD-check](https://github.com/natverse/nat.ggplot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/natverse/nat.ggplot/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# nat.ggplot

`nat.ggplot` enables the [Neuroanatomy Toolbox (nat)](https://jefferis.github.io/nat/) suite to create publication-quality 2D visualisations of neurons and brain meshes using [ggplot2](https://ggplot2.tidyverse.org/). This package bridges the gap between nat's 3D neuroanatomy capabilities and ggplot2's 2D plotting framework, allowing users to create figures of neuronal morphology data.

## About Neuron Objects

In the [natverse ecosystem](https://natverse.org/), neurons are represented as tree structures containing 3D coordinates and connectivity information. The main data structures include:

- **neuron**: A single neuron with XYZ coordinates and parent-child connectivity defining the tree structure, `neuron$d` provides data in the style of a .swc file, where `neuron$d$Label` can give the arbour type, e.g. dendrite versus axon.
- **neuronlist**: A collection of neurons that can be manipulated together, bundled as a list of class 'neuronlist'
- **mesh3d**: 3D surface meshes representing brain regions or neuron surfaces
- **Synaptic information**: Pre- and postsynaptic site locations that can be overlaid on morphology

For more details, see the [natverse neurons introduction](https://natverse.org/nat/articles/neurons-intro.html).

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
p <- ggneuron(banc.meshes_lowres, rotation_matrix = banc_view)
p
ggsave(file.path(output_dir, "quickstart_basic.png"), p, width = 6, height = 6, dpi = 150, bg = "white")
```

![Basic visualisation](https://github.com/natverse/nat.ggplot/blob/main/inst/images/quickstart_basic.png?raw=true)

``` r
# Customise visualisation with g.anat base
p <- g.anat + 
  geom_neuron(banc.brain_neuropil_lowres, 
              rotation_matrix = banc_view,
              cols = c("grey90", "grey60"),
              alpha = 0.3) +
  geom_neuron(banc.meshes_lowres,
              rotation_matrix = banc_view,
              cols = c("purple", "magenta")) +
  geom_neuron(banc.skels,
              rotation_matrix = banc_view,
              cols = c("black"))
p
ggsave(file.path(output_dir, "quickstart_custom.png"), p, width = 6, height = 6, dpi = 150, bg = "white")
```

![Customised visualisation](https://github.com/natverse/nat.ggplot/blob/main/inst/images/quickstart_custom.png?raw=true)

## Example Gallery

### Multiple Views of LHPD2a1 Neurons

The package includes several pre-defined views for BANC connectome data. Here are the LHPD2a1 neurons viewed from different angles:

``` r
# Define the different view matrices from bancr
views <- list(
  main = structure(c(0.961547076702118, 0.037275392562151, 
    0.27209860086441, 0, 0.0369537360966206, -0.999296963214874, 
    0.00630810856819153, 0, 0.272142440080643, 0.00398948788642883, 
    -0.962248742580414, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
  side = structure(c(0.188666880130768, 0.137750864028931, 
    -0.972331881523132, 0, 0.130992725491524, -0.98479551076889, 
    -0.114099271595478, 0, -0.97326534986496, -0.105841755867004, 
    -0.203842639923096, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
  front = structure(c(0.99931389093399, 0.0139970388263464, 
    -0.0342894680798054, 0, -0.0321401171386242, -0.132316529750824, 
    -0.990686297416687, 0, -0.0184037387371063, 0.991108655929565, 
    -0.131775915622711, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
  brain_side = structure(c(0.000764884985983372, 0.0153511334210634, 
    -0.99988180398941, 0, -0.940421104431152, -0.339961022138596, 
    -0.00593886896967888, 0, -0.340011894702911, 0.94031423330307, 
    0.0141764245927334, 0, -401395.405539944, -128785.809090088, 
    -5607.3408203126, 1), dim = c(4L, 4L))
)

# Create plots for each view
plots <- list()
for(view_name in names(views)) {
  plots[[view_name]] <- ggneuron(
    banc.skels,
    volume = banc.brain_neuropil_lowres,
    rotation_matrix = views[[view_name]],
    cols1 = c("#8B1A89", "#FF69B4"),
    cols2 = c("grey95", "grey85"),
    alpha = 0.8,
    info = paste(toupper(substring(view_name, 1, 1)), 
                 substring(view_name, 2), " view", sep="")
  )
}

# Display and save individual views
plots$main
ggsave(file.path(output_dir, "view_main.png"), plots$main, width = 5, height = 5, dpi = 150, bg = "white")

plots$side
ggsave(file.path(output_dir, "view_side.png"), plots$side, width = 5, height = 5, dpi = 150, bg = "white")

plots$front
ggsave(file.path(output_dir, "view_front.png"), plots$front, width = 5, height = 5, dpi = 150, bg = "white")

plots$brain_side
ggsave(file.path(output_dir, "view_brain_side.png"), plots$brain_side, width = 5, height = 5, dpi = 150, bg = "white")

# Create a collage using patchwork (if available)
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  collage <- (plots$main | plots$side) / (plots$front | plots$brain_side)
  collage
  ggsave(file.path(output_dir, "views_collage.png"), collage, width = 10, height = 10, dpi = 150, bg = "white")
}
```

![Multiple views collage](https://github.com/natverse/nat.ggplot/blob/main/inst/images/views_collage.png?raw=true)

### Finding the Right View

The package includes tools to interactively find a suitable viewing angle for your neurons:

``` r
library(nat)  # for plot3d
library(rgl)  # for 3D interaction

# Step 1: Plot neurons in 3D and rotate to desired view
plot3d(banc.brain_neuropil_lowres, alpha = 0.3, col = "grey")
plot3d(banc.skels)

# Step 2: Rotate with mouse to find the desired angle
# Step 3: Capture the view
my_view <- rgl_view()

# Step 4: Use the captured view in ggplot2
ggneuron(banc.skels, rotation_matrix = my_view$userMatrix)
```

### Visualising Brain Neuropil

You can visualise the brain neuropil mesh with different colours and transparency:

``` r
# Brain neuropil with custom colours
p <- ggneuron(banc.brain_neuropil_lowres,
              rotation_matrix = banc_view,
              cols1 = c("lightblue", "darkblue"), 
              alpha = 0.5)
p
ggsave(file.path(output_dir, "brain_neuropil.png"), p, width = 6, height = 6, dpi = 150, bg = "white")
```

![Brain neuropil](https://github.com/natverse/nat.ggplot/blob/main/inst/images/brain_neuropil.png?raw=true)

### Plotting Neurons with Custom Colours

Visualise multiple neurons with different colouring schemes:

``` r
# All neurons in one colour
p <- g.anat +
  geom_neuron(banc.skels,
              rotation_matrix = banc_view,
              cols = c("navy", "navy"))
p
ggsave(file.path(output_dir, "neurons_single_colour.png"), p, width = 6, height = 6, dpi = 150, bg = "white")
```

![Single colour](https://github.com/natverse/nat.ggplot/blob/main/inst/images/neurons_single_colour.png?raw=true)

``` r
# Colour gradient from dendrite to axon tips
p <- g.anat +
  geom_neuron(banc.skels,
              rotation_matrix = banc_view,
              cols = c("purple", "orange"))
p
ggsave(file.path(output_dir, "neurons_gradient.png"), p, width = 6, height = 6, dpi = 150, bg = "white")
```

![Gradient colouring](https://github.com/natverse/nat.ggplot/blob/main/inst/images/neurons_gradient.png?raw=true)

``` r
# Each neuron in a different colour
neuron_colours <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8")
p <- g.anat +
  geom_neuron(banc.brain_neuropil_lowres,
              rotation_matrix = banc_view,
              cols = c("grey95", "grey80"),
              alpha = 0.2) +
  geom_neuron(banc.skels,
              rotation_matrix = banc_view,
              cols = neuron_colours)
p
ggsave(file.path(output_dir, "neurons_multi_colour.png"), p, width = 6, height = 6, dpi = 150, bg = "white")
```

![Multiple colours](https://github.com/natverse/nat.ggplot/blob/main/inst/images/neurons_multi_colour.png?raw=true)

### Visualising Synapses

Display synaptic sites as coloured points:

``` r
# Create synapse visualisation
p <- g.anat +
  # Brain context
  geom_neuron(banc.brain_neuropil_lowres,
              rotation_matrix = banc_view,
              cols = c("grey95", "grey85"),
              alpha = 0.3) +
  # Neuron skeletons
  geom_neuron(banc.skels,
              rotation_matrix = banc_view,
              cols = c("grey60", "grey40"),
              alpha = 0.8) +
  # Synapses
  geom_point(data = banc.syns %>%
               mutate(rotated = as.data.frame(t(banc_view[,1:3] %*% t(.[,c("X","Y","Z")])))) %>%
               cbind(rotated),
             aes(x = V1, y = V2, colour = factor(prepost)),
             size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = c("0" = "#DC143C", "1" = "#000080"),
                     labels = c("Presynaptic", "Postsynaptic"),
                     name = "Synapse Type") +
  theme(legend.position = "bottom")
p
ggsave(file.path(output_dir, "synapses_visualisation.png"), p, width = 7, height = 7, dpi = 150, bg = "white")
```

![Synapse visualisation](https://github.com/natverse/nat.ggplot/blob/main/inst/images/synapses_visualisation.png?raw=true)

### Axon-Dendrite Split Visualisation

Visualise neurons with compartments coloured by their functional identity (using flow centrality from Schneider-Mizell et al., 2016):

``` r
# Single neuron with axon/dendrite split
p <- ggneuron(bc.neurons.flow[[1]], 
              rotation_matrix = banc_view,
              info = paste("Neuron", names(bc.neurons.flow)[1]))
p
ggsave(file.path(output_dir, "split_single_neuron.png"), p, width = 6, height = 6, dpi = 150, bg = "white")
```

![Single split neuron](https://github.com/natverse/nat.ggplot/blob/main/inst/images/split_single_neuron.png?raw=true)

``` r
# All split neurons with brain context
p <- ggneuron(bc.neurons.flow,
              volume = banc.brain_neuropil_lowres,
              rotation_matrix = banc_view,
              info = "LHPD2a1 neurons: Axon (orange) vs Dendrite (cyan)")
p
ggsave(file.path(output_dir, "split_all_neurons.png"), p, width = 7, height = 7, dpi = 150, bg = "white")
```

![All split neurons](https://github.com/natverse/nat.ggplot/blob/main/inst/images/split_all_neurons.png?raw=true)

### Creating Figures

Combine multiple elements for a figure:

``` r
library(ggplot2)

# Create comprehensive visualisation
p <- g.anat +
  # Add brain outline
  geom_neuron(banc.brain_neuropil_lowres,
              rotation_matrix = banc_view,
              cols = c("grey95", "grey85"),
              alpha = 0.3) +
  # Add neuron meshes
  geom_neuron(banc.meshes_lowres,
              rotation_matrix = banc_view,
              cols = c("#E8B4F3", "#B14A9A"),
              alpha = 0.7) +
  # Add title
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(title = "LHPD2a1 Neurons",
       subtitle = "Lateral horn local neurons from BANC connectome")

# Save figure
ggsave("lhpd2a1_neurons.pdf", p, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_dir, "figure_comprehensive.png"), p, width = 7, height = 7, dpi = 150, bg = "white")
```

![Comprehensive figure](https://github.com/natverse/nat.ggplot/blob/main/inst/images/figure_comprehensive.png?raw=true)

## Data Attribution

This package includes example data from the BANC (Brain And Nerve Cord) connectome:

- **LHPD2a1 neurons**: First characterized by Dolan et al. (2018) in *Neuron*
- **BANC connectome data**: From Bates et al. (2025) in *bioRxiv*

## Advanced Features

### Custom Rotation Matrices

Create and save custom views for consistent visualisations:

``` r
# Define multiple standard views
views <- list(
  front = banc_view,
  top = rgl_view()$userMatrix,  # After rotating in rgl
  side = rgl_view()$userMatrix   # After rotating in rgl
)

# Apply views consistently across figures
for(view_name in names(views)) {
  p <- ggneuron(banc.skels, 
                rotation_matrix = views[[view_name]],
                info = paste(view_name, "view"))
  ggsave(paste0("neurons_", view_name, ".png"), p)
}
```

### Integration with ggplot2 Ecosystem

The package works with ggplot2 extensions:

``` r
library(ggplot2)
library(patchwork)  # For combining plots

# Create multiple panels
p1 <- g.anat + 
  geom_neuron(banc.skels[[1]], rotation_matrix = banc_view) +
  ggtitle("Neuron 1")
ggsave(file.path(output_dir, "panel_neuron1.png"), p1, width = 4, height = 4, dpi = 150, bg = "white")

p2 <- g.anat + 
  geom_neuron(banc.skels[[2]], rotation_matrix = banc_view) +
  ggtitle("Neuron 2")
ggsave(file.path(output_dir, "panel_neuron2.png"), p2, width = 4, height = 4, dpi = 150, bg = "white")

# Combine with patchwork
p_combined <- p1 + p2 + plot_layout(ncol = 2)
p_combined
ggsave(file.path(output_dir, "panels_combined.png"), p_combined, width = 8, height = 4, dpi = 150, bg = "white")
```

![Combined panels](https://github.com/natverse/nat.ggplot/blob/main/inst/images/panels_combined.png?raw=true)

## Acknowledgements

This package was developed by Alexander Shakeel Bates while in the laboratory of Rachel I. Wilson at Harvard Medical School. The package leverages the [natverse](https://natverse.github.io) ecosystem for neuroanatomy in R.

## References

Key papers referenced in this package:

- Bates, A. S., Phelps, J. S., Kim, M., Yang, H. H., Matsliah, A., Ajabi, Z., Perlman, E., et al. (2025). Distributed control circuits across a brain-and-cord connectome. *bioRxiv*. doi: [10.1101/2025.07.31.667571](https://doi.org/10.1101/2025.07.31.667571)

- Dolan, M. J., et al. (2018). Communication from learned to innate olfactory processing centers is required for memory retrieval in Drosophila. *Neuron*, 100(3), 651-668. doi: [10.1016/j.neuron.2018.08.037](https://doi.org/10.1016/j.neuron.2018.08.037)

- Jefferis, G. S. X. E., et al. (2007). Comprehensive maps of Drosophila higher olfactory centers. *Cell*, 128(6), 1187-1203. doi: [10.1016/j.cell.2007.01.040](https://doi.org/10.1016/j.cell.2007.01.040)

- Schneider-Mizell, C. M., et al. (2016). Quantitative neuroanatomy for connectomics in Drosophila. *eLife*, 5, e12059. doi: [10.7554/eLife.12059](https://doi.org/10.7554/eLife.12059)

## Getting Help

- For issues and feature requests, please use the [GitHub issue tracker](https://github.com/natverse/nat.ggplot/issues)
- For general natverse questions, visit the [natverse documentation](https://natverse.github.io)
- For BANC connectome data, see the [bancr package](https://github.com/natverse/bancr)
