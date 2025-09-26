# Basic tests for core functionality

test_that("Package loads and has expected functions", {
  # Check main functions exist
  expect_true(exists("geom_neuron"))
  expect_true(exists("ggneuron"))
  expect_true(exists("ggplot2_neuron_path"))
  expect_true(exists("rgl_view"))
  expect_true(exists("g.anat"))
})

test_that("Basic neuron plotting works", {
  skip_if_not(requireNamespace("nat", quietly = TRUE))
  
  # Load data
  data(banc.skels)
  data(banc_view)
  
  # Create basic plot
  p <- g.anat + geom_neuron(banc.skels[[1]], rotation_matrix = banc_view)
  expect_s3_class(p, "ggplot")
})

test_that("geom_neuron handles different input types", {
  skip_if_not(requireNamespace("nat", quietly = TRUE))
  
  # Load data
  data(banc.skels)
  data(banc.syns)
  
  # Test with neuron
  result_neuron <- geom_neuron(banc.skels[[1]])
  expect_type(result_neuron, "list")
  
  # Test with neuronlist
  result_list <- geom_neuron(banc.skels[1:2])
  expect_type(result_list, "list")
  
  # Test with matrix (synapses)
  syn_matrix <- as.matrix(banc.syns[1:50, c("X", "Y", "Z")])
  result_matrix <- geom_neuron(syn_matrix, root = 0.5)
  expect_type(result_matrix, "list")
})

test_that("Rotation matrix application works", {
  skip_if_not(requireNamespace("nat", quietly = TRUE))
  
  # Load data
  data(banc.skels)
  data(banc_view)
  
  # Test ggplot2_neuron_path with rotation
  original <- ggplot2_neuron_path(banc.skels[[1]])
  rotated <- ggplot2_neuron_path(banc.skels[[1]], rotation_matrix = banc_view)
  
  # Coordinates should differ
  expect_false(identical(original$X, rotated$X))
  expect_false(identical(original$Y, rotated$Y))
  
  # Structure should be same
  expect_equal(nrow(original), nrow(rotated))
})

test_that("ggneuron wrapper creates complete plots", {
  skip_if_not(requireNamespace("nat", quietly = TRUE))
  
  # Load data
  data(banc.skels)
  data(banc_view)
  
  # Create plot with ggneuron
  p <- ggneuron(banc.skels[[1]], rotation_matrix = banc_view)
  expect_s3_class(p, "ggplot")
  
  # Test with title
  p_titled <- ggneuron(banc.skels[[1]], info = "Test Neuron")
  expect_equal(p_titled$labels$title, "Test Neuron")
})