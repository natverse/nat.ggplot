test_that("Complete workflow produces valid plots", {
  # Test complete workflow from data to plot
  
  # 1. Simple neuron plot
  p1 <- g.anat + 
    geom_neuron(banc.skels[[1]], rotation_matrix = banc_view)
  expect_s3_class(p1, "ggplot")
  
  # 2. Multiple neurons with colours
  p2 <- g.anat +
    geom_neuron(banc.skels, 
                rotation_matrix = banc_view,
                cols = c("purple", "orange"))
  expect_s3_class(p2, "ggplot")
  
  # 3. Brain mesh with neurons
  p3 <- g.anat +
    geom_neuron(banc.brain_neuropil,
                rotation_matrix = banc_view,
                cols = c("grey90", "grey70"),
                alpha = 0.3) +
    geom_neuron(banc.skels,
                rotation_matrix = banc_view)
  expect_s3_class(p3, "ggplot")
  
  # 4. Synapses as points
  synapse_coords <- as.matrix(banc.syns[1:100, c("X", "Y", "Z")])
  p4 <- g.anat +
    geom_neuron(synapse_coords,
                rotation_matrix = banc_view,
                root = 0.5,
                cols = c("red", "blue"))
  expect_s3_class(p4, "ggplot")
  
  # 5. Split neurons
  if (length(bc.neurons.flow) > 0) {
    p5 <- ggneuron(bc.neurons.flow[[1]], 
                   rotation_matrix = banc_view)
    expect_s3_class(p5, "ggplot")
  }
})

test_that("Plots can be saved without error", {
  skip_if_not(dir.exists(tempdir()))
  
  # Create a simple plot
  p <- g.anat + geom_neuron(banc.skels[[1]])
  
  # Test saving to temporary file
  temp_file <- tempfile(fileext = ".png")
  
  # This should not error
  expect_silent(
    ggsave(temp_file, p, width = 5, height = 5, dpi = 100)
  )
  
  # Check file was created
  expect_true(file.exists(temp_file))
  
  # Clean up
  unlink(temp_file)
})

test_that("Different neuron types can be combined", {
  # Combine different object types in one plot
  p <- g.anat +
    # Mesh
    geom_neuron(banc.brain_neuropil,
                rotation_matrix = banc_view,
                cols = c("grey95", "grey85"),
                alpha = 0.3) +
    # Skeletons
    geom_neuron(banc.skels[1:2],
                rotation_matrix = banc_view,
                cols = c("blue", "green")) +
    # Points (synapses)
    geom_neuron(as.matrix(banc.syns[1:50, c("X", "Y", "Z")]),
                rotation_matrix = banc_view,
                root = 0.3,
                cols = c("red", "darkred"))
  
  expect_s3_class(p, "ggplot")
  
  # Should have multiple layers
  expect_true(length(p$layers) > 3)
})

test_that("Size parameter affects line width", {
  # Create plots with different sizes
  p_thin <- g.anat + 
    geom_neuron(banc.skels[[1]], size = 0.1)
  p_thick <- g.anat + 
    geom_neuron(banc.skels[[1]], size = 2)
  
  expect_s3_class(p_thin, "ggplot")
  expect_s3_class(p_thick, "ggplot")
  
  # Both should have layers
  expect_true(length(p_thin$layers) > 0)
  expect_true(length(p_thick$layers) > 0)
})

test_that("Coordinate limits can be set", {
  # Test zooming with coord_fixed
  p <- g.anat +
    geom_neuron(banc.skels[[1]], rotation_matrix = banc_view) +
    coord_fixed(xlim = c(-1000, 1000), ylim = c(-1000, 1000))
  
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$coordinates, "CoordFixed")
  expect_equal(p$coordinates$limits$x, c(-1000, 1000))
  expect_equal(p$coordinates$limits$y, c(-1000, 1000))
})

test_that("Empty or invalid data handled gracefully", {
  # Empty neuronlist
  empty_nl <- nat::neuronlist()
  p1 <- g.anat + geom_neuron(empty_nl)
  expect_s3_class(p1, "ggplot")
  
  # NULL input
  p2 <- g.anat + geom_neuron(NULL)
  expect_s3_class(p2, "ggplot")
  
  # Empty matrix
  empty_matrix <- matrix(ncol = 3, nrow = 0)
  colnames(empty_matrix) <- c("X", "Y", "Z")
  p3 <- g.anat + geom_neuron(empty_matrix)
  expect_s3_class(p3, "ggplot")
})
