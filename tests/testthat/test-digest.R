test_that(
  "digest_protein reports scaffold state with an informative classed error",
  {
    expect_error(
      digest_protein("MKWVTFISLLFLFSSAYSR"),
      regexp = "does not implement",
      class = "pepvet_not_implemented"
    )
  }
)

test_that("digest_protein reports scaffold state for documented input shapes", {
  expect_error(
    digest_protein(
      Biostrings::AAStringSet("MKWVTFISLLFLFSSAYSR"),
      enzyme = "trypsin",
      min_length = 5L,
      max_length = 25L
    ),
    class = "pepvet_not_implemented"
  )
})
