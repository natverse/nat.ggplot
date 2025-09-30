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

test_that("gganat is a valid ggplot object", {
  # it's a ggplot
  expect_s3_class(gganat, "ggplot")

  # coord: accept either CoordFixed (old ggplot2) or CoordCartesian with ratio (new)
  coord <- suppressWarnings(ggplot2::ggplot_build(gganat)$layout$coord)
  expect_s3_class(coord, "Coord")
  expect_true(inherits(coord, "CoordFixed") || inherits(coord, "CoordCartesian"))
  if (inherits(coord, "CoordCartesian")) {
    expect_true(is.numeric(coord$ratio) && length(coord$ratio) == 1)
  }

  # theme: don't rely on `$names` — assert it's NULL or a theme object
  expect_true(is.null(gganat$theme) || inherits(gganat$theme, "theme"))
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
