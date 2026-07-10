# Tests aligned with the corresponding plotting source domain.
test_that("plot_enzyme_comparison returns a patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p <- plot_enzyme_comparison(comp)
  expect_s3_class(p, "patchwork")
})

test_that("plot_enzyme_comparison: recommend = FALSE suppresses badge", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  expect_no_error(plot_enzyme_comparison(comp, recommend = FALSE))
})

test_that("plot_enzyme_comparison: subset of scores renders without error", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  expect_no_error(
    plot_enzyme_comparison(comp, scores = c("S_coverage", "S_count"))
  )
})

test_that("plot_enzyme_comparison: custom title is accepted", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p <- plot_enzyme_comparison(comp, title = "My comparison")
  expect_s3_class(p, "patchwork")
})

test_that("plot_enzyme_comparison: two-enzyme comparison works", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison(enzymes = c("trypsin", "lysc"))
  expect_no_error(plot_enzyme_comparison(comp))
})

test_that("plot_enzyme_comparison errors on non-data.frame input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  expect_error(
    plot_enzyme_comparison("not a tibble"),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_enzyme_comparison errors on missing required columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  bad <- data.frame(enzyme = c("a", "b"), foo = 1:2)
  expect_error(
    plot_enzyme_comparison(bad),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_enzyme_comparison errors when fewer than 2 enzymes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  bad <- data.frame(
    enzyme = "trypsin", composite_score = 0.6,
    S_coverage = 0.8, S_length = 0.7,
    S_count = 0.5, S_hydro = 0.6, S_charge = 0.4
  )
  expect_error(
    plot_enzyme_comparison(bad),
    class = "pepvet_error_invalid_comparison"
  )
})
