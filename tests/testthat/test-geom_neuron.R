test_that("geom_neuron returns ggplot2 layers", {
  # Test with single neuron
  result <- geom_neuron(banc.skels[[1]])
  
  expect_type(result, "list")
  expect_true(length(result) > 0)
  
  # Check that elements are valid ggplot2 or ggnewscale objects
  for (elem in result) {
    if (!is.null(elem)) {
      # Elements can be layers, scales, ggnewscale objects, or other ggplot2 objects
      expect_true(
        inherits(elem, "Layer") || 
        inherits(elem, "Scale") || 
        inherits(elem, "ggproto") ||
        inherits(elem, "gg") ||
        inherits(elem, "new_aes") ||  # ggnewscale objects
        is.list(elem)  # Some elements may be lists of ggplot2 objects
      )
    }
  }
})

test_that("geom_neuron works with different object types", {
  # Test with neuron
  neuron_result <- geom_neuron(banc.skels[[1]])
  expect_type(neuron_result, "list")
  
  # Test with neuronlist
  neuronlist_result <- geom_neuron(banc.skels)
  expect_type(neuronlist_result, "list")
  
  # Test with mesh3d
  mesh_result <- geom_neuron(banc.brain_neuropil)
  expect_type(mesh_result, "list")
  
  # Test with matrix (synapse points)
  matrix_result <- geom_neuron(as.matrix(banc.syns[1:100, c("X", "Y", "Z")]))
  expect_type(matrix_result, "list")
  
  # Test with data.frame
  df_result <- geom_neuron(banc.syns[1:100, c("X", "Y", "Z")])
  expect_type(df_result, "list")
  
  # Test with NULL
  null_result <- geom_neuron(NULL)
  expect_type(null_result, "list")
})

test_that("geom_neuron respects parameters", {
  # Test size parameter
  result_thin <- geom_neuron(banc.skels[[1]], size = 0.2)
  result_thick <- geom_neuron(banc.skels[[1]], size = 1.0)
  
  # Both should return valid layer lists
  expect_type(result_thin, "list")
  expect_type(result_thick, "list")
  
  # Test root parameter
  result_no_root <- geom_neuron(banc.skels[[1]], root = 0)
  result_with_root <- geom_neuron(banc.skels[[1]], root = 5)
  
  expect_type(result_no_root, "list")
  expect_type(result_with_root, "list")
  
  # Test colour parameters
  result_single_col <- geom_neuron(banc.skels, cols = "purple")
  result_gradient <- geom_neuron(banc.skels[[1]], cols = c("blue", "red"))
  result_rainbow <- geom_neuron(banc.skels, cols = "rainbow")
  
  expect_type(result_single_col, "list")
  expect_type(result_gradient, "list")
  expect_type(result_rainbow, "list")
})

test_that("geom_neuron works with rotation matrices", {
  # Test with rotation
  result <- geom_neuron(banc.skels[[1]], rotation_matrix = banc_view)
  
  expect_type(result, "list")
  expect_true(length(result) > 0)
})

test_that("geom_neuron handles split neurons", {
  # Test with split neuron
  if (length(banc.neurons.flow) > 0) {
    result <- geom_neuron(banc.neurons.flow[[1]])
    
    expect_type(result, "list")
    expect_true(length(result) > 0)
    
    # Should have multiple path elements for different compartments
    path_elements <- sapply(result, function(x) inherits(x, "Layer") && !is.null(x$geom) && inherits(x$geom, "GeomPath"))
    expect_true(sum(path_elements) > 1)  # Multiple compartments
  }
})

test_that("geom_neuron can be added to g.anat", {
  # Test that geom_neuron works with g.anat base
  p <- g.anat + geom_neuron(banc.skels[[1]])
  
  expect_s3_class(p, "ggplot")
  
  # Test multiple layers
  p2 <- g.anat + 
    geom_neuron(banc.brain_neuropil, cols = c("grey90", "grey60"), alpha = 0.3) +
    geom_neuron(banc.skels)
  
  expect_s3_class(p2, "ggplot")
})

test_that("geom_neuron handles edge cases", {
  # Empty neuronlist
  empty_nl <- nat::neuronlist()
  result <- geom_neuron(empty_nl)
  expect_type(result, "list")
  
  # Single point matrix
  single_point <- matrix(c(0, 0, 0), nrow = 1, ncol = 3)
  colnames(single_point) <- c("X", "Y", "Z")
  result <- geom_neuron(single_point)
  expect_type(result, "list")
})
