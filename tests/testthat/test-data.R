test_that("banc.skels is a valid neuronlist", {
  expect_s3_class(banc.skels, "neuronlist")
  expect_true(length(banc.skels) > 0)
  
  # Check first neuron
  expect_s3_class(banc.skels[[1]], "neuron")
  expect_true(is.data.frame(banc.skels[[1]]$d))
  expect_true(all(c("X", "Y", "Z", "Parent") %in% colnames(banc.skels[[1]]$d)))
})

test_that("banc.meshes is a valid neuronlist of mesh3d objects", {
  expect_s3_class(banc.meshes, "neuronlist")
  expect_true(length(banc.meshes) > 0)
  
  # Check first mesh
  expect_s3_class(banc.meshes[[1]], "mesh3d")
  expect_true(is.matrix(banc.meshes[[1]]$vb))
  expect_true(is.matrix(banc.meshes[[1]]$it))
  expect_equal(nrow(banc.meshes[[1]]$vb), 4)  # Homogeneous coordinates
  expect_equal(nrow(banc.meshes[[1]]$it), 3)  # Triangle indices
})

test_that("banc.brain_neuropil is a valid mesh3d", {
  expect_s3_class(banc.brain_neuropil, "mesh3d")
  expect_true(is.matrix(banc.brain_neuropil$vb))
  expect_true(is.matrix(banc.brain_neuropil$it))
  expect_equal(nrow(banc.brain_neuropil$vb), 4)
  expect_equal(nrow(banc.brain_neuropil$it), 3)
  expect_true(ncol(banc.brain_neuropil$vb) > 0)  # Has vertices
  expect_true(ncol(banc.brain_neuropil$it) > 0)  # Has faces
})

test_that("banc.syns is a valid data.frame", {
  expect_s3_class(banc.syns, "data.frame")
  expect_true(nrow(banc.syns) > 0)
  
  # Check required columns
  expect_true(all(c("X", "Y", "Z", "prepost") %in% colnames(banc.syns)))
  
  # Check coordinate columns are numeric
  expect_true(is.numeric(banc.syns$X))
  expect_true(is.numeric(banc.syns$Y))
  expect_true(is.numeric(banc.syns$Z))
  
  # Check prepost is binary
  expect_true(all(banc.syns$prepost %in% c(0, 1)))
})

test_that("bc.neurons.flow contains split neurons", {
  expect_s3_class(bc.neurons.flow, "neuronlist")
  expect_true(length(bc.neurons.flow) > 0)
  
  # Check first neuron has split labels
  first_neuron <- bc.neurons.flow[[1]]
  expect_s3_class(first_neuron, "neuron")
  expect_true("Label" %in% colnames(first_neuron$d))
  
  # Check for different compartment labels
  labels <- unique(first_neuron$d$Label)
  # Should have multiple compartment types
  expect_true(length(labels) > 1)
})

test_that("banc_view is a valid rotation matrix", {
  expect_true(is.matrix(banc_view))
  expect_equal(dim(banc_view), c(4, 4))
  expect_true(is.numeric(banc_view))
  
  # Last row should be [0, 0, 0, 1] for homogeneous coordinates
  expect_equal(banc_view[4, ], c(0, 0, 0, 1))
  
  # Check rotation part is roughly orthonormal (allowing for numerical error)
  rotation_part <- banc_view[1:3, 1:3]
  identity_approx <- rotation_part %*% t(rotation_part)
  expect_true(all(abs(diag(identity_approx) - 1) < 0.01))  # Diagonal close to 1
})

test_that("all data objects have consistent BANC IDs", {
  # Check that neuron names are consistent
  expect_equal(names(banc.skels), names(banc.meshes))
  expect_equal(names(banc.skels), names(bc.neurons.flow))
  
  # All should be BANC IDs (large numbers as strings)
  ids <- names(banc.skels)
  expect_true(all(nchar(ids) > 10))  # BANC IDs are long
})