test_that("digest_protein reports scaffold state", {
  expect_error(
    digest_protein("MKWVTFISLLFLFSSAYSR"),
    class = "pepvet_not_implemented"
  )
})
