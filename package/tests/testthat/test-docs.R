test_that("every man page parses cleanly", {
  # A `\%` in markdown roxygen becomes `\\%` in Rd, whose `%` then comments
  # out the rest of the line and can close a section early; R CMD INSTALL only
  # warns ("unknown macro '\item'"). Runs from a source checkout.
  man <- testthat::test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "no man/ directory (installed package)")
  for (f in list.files(man, pattern = "\\.Rd$", full.names = TRUE)) {
    problems <- tools::checkRd(f)
    expect_length(problems, 0L)
    if (length(problems)) message(basename(f), ": ", paste(problems, collapse = "; "))
  }
})
