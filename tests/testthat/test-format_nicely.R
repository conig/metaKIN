test_that("transf works in format_nicely", {
  require(metaKIN)
  dat <- metaSEM::Bornmann07

  dat$study_year <- stringr::str_extract(dat$Study, "[0-9]{1,4}")
  dat$study_author <- gsub("\\s.*", "", dat$Study)

  m <- meta3L(logOR, v, cluster = Cluster, data = dat) |>
    moderate(Country = ~ Country - 1)

  testthat::expect_error(
    format_nicely(m, transform = function(y) y * 2)
  )

  testthat::expect_no_error(
    nice_tab <- format_nicely(
      m,
      transf = function(x) 999,
      transf_name = "999"
    )
  )
  testthat::expect_true(grepl("999", nice_tab$"999"[1]))
})
