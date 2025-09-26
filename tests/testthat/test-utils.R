test_that("rgl_view returns expected structure", {
  # Test that rgl_view returns a list with userMatrix
  result <- rgl_view()
  
  expect_type(result, "list")
  expect_true("userMatrix" %in% names(result))
  
  if (!is.null(result$userMatrix)) {
    expect_true(is.matrix(result$userMatrix))
    expect_equal(dim(result$userMatrix), c(4, 4))
  }
})

test_that("g.anat is a valid ggplot object", {
  expect_s3_class(g.anat, "ggplot")
  
  # Check that it has expected theme elements
  expect_true("theme" %in% names(g.anat))
  expect_true("coordinates" %in% names(g.anat))
  
  # Check coord_fixed is applied
  expect_s3_class(g.anat$coordinates, "CoordFixed")
  
  # Can add layers to it
  p <- g.anat + ggplot2::geom_point(data = data.frame(x = 1, y = 1), aes(x, y))
  expect_s3_class(p, "ggplot")
})

test_that("is_named_all works correctly", {
  # Test with named vector
  named_vec <- c(a = 1, b = 2, c = 3)
  expect_true(is_named_all(named_vec))
  
  # Test with partially named vector
  partial_named <- c(a = 1, 2, c = 3)
  expect_false(is_named_all(partial_named))
  
  # Test with unnamed vector
  unnamed_vec <- c(1, 2, 3)
  expect_false(is_named_all(unnamed_vec))
  
  # Test with named list
  named_list <- list(a = 1, b = 2)
  expect_true(is_named_all(named_list))
  
  # Test with NULL
  expect_false(is_named_all(NULL))
  
  # Test with empty vector
  expect_false(is_named_all(c()))
})

test_that("ggneuron wrapper function works", {
  # Test basic usage
  p <- ggneuron(banc.skels[[1]])
  expect_s3_class(p, "ggplot")
  
  # Test with volume
  p_with_volume <- ggneuron(
    banc.skels[[1]], 
    volume = banc.brain_neuropil
  )
  expect_s3_class(p_with_volume, "ggplot")
  
  # Test with rotation
  p_rotated <- ggneuron(
    banc.skels, 
    rotation_matrix = banc_view
  )
  expect_s3_class(p_rotated, "ggplot")
  
  # Test with custom colours
  p_colours <- ggneuron(
    banc.skels[[1]],
    cols1 = c("purple", "orange"),
    cols2 = c("grey80", "grey60")
  )
  expect_s3_class(p_colours, "ggplot")
  
  # Test with title
  p_title <- ggneuron(
    banc.skels[[1]],
    info = "Test neuron"
  )
  expect_s3_class(p_title, "ggplot")
  expect_equal(p_title$labels$title, "Test neuron")
})

test_that("ggneuron passes additional arguments", {
  # Test that size parameter is passed through
  p_size <- ggneuron(banc.skels[[1]], size = 2)
  expect_s3_class(p_size, "ggplot")
  
  # Test that alpha is applied
  p_alpha <- ggneuron(banc.skels[[1]], alpha = 0.2)
  expect_s3_class(p_alpha, "ggplot")
})
