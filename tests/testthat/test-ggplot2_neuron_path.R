test_that("ggplot2_neuron_path works with neuron objects", {
  # Test with single neuron
  result <- ggplot2_neuron_path(banc.skels[[1]])
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("X", "Y", "Z", "group") %in% colnames(result)))
  expect_true(nrow(result) > 0)
  expect_true(is.numeric(result$X))
  expect_true(is.numeric(result$Y))
  expect_true(is.numeric(result$Z))
  expect_true(is.numeric(result$group))
})

test_that("ggplot2_neuron_path works with neuronlist objects", {
  # Test with neuronlist
  result <- ggplot2_neuron_path(banc.skels)
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("X", "Y", "Z", "group", "id") %in% colnames(result)))
  expect_true(nrow(result) > 0)
  expect_true(length(unique(result$id)) == length(banc.skels))
})

test_that("ggplot2_neuron_path works with mesh3d objects", {
  # Test with mesh3d
  result <- ggplot2_neuron_path(banc.brain_neuropil)
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("X", "Y", "Z", "group") %in% colnames(result)))
  expect_true(nrow(result) > 0)
  # Mesh should create triangles (3 points per group)
  group_counts <- table(result$group)
  expect_true(all(group_counts == 3))
})

test_that("ggplot2_neuron_path handles NULL input", {
  result <- ggplot2_neuron_path(NULL)
  expect_null(result)
})

test_that("ggplot2_neuron_path applies rotation matrix correctly", {
  # Test rotation application
  original <- ggplot2_neuron_path(banc.skels[[1]])
  rotated <- ggplot2_neuron_path(banc.skels[[1]], rotation_matrix = banc_view)
  
  # Coordinates should be different after rotation
  expect_false(isTRUE(all.equal(original$X, rotated$X)))
  expect_false(isTRUE(all.equal(original$Y, rotated$Y)))
  
  # But structure should be the same
  expect_equal(nrow(original), nrow(rotated))
  expect_equal(original$group, rotated$group)
})

test_that("ggplot2_neuron_path handles invalid mesh3d objects", {
  # Test with empty mesh
  empty_mesh <- list(
    vb = matrix(0, nrow = 4, ncol = 0),
    it = matrix(0, nrow = 3, ncol = 0)
  )
  class(empty_mesh) <- "mesh3d"
  
  expect_warning(result <- ggplot2_neuron_path(empty_mesh))
  expect_s3_class(result, "data.frame")
  expect_true(is.na(result$X[1]))
})
