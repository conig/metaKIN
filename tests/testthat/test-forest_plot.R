test_that("transf works in forest plot", {

  require(metaKIN)
  dat <- metaSEM::Bornmann07

  dat$study_year <- stringr::str_extract(dat$Study, "[0-9]{1,4}")
  dat$study_author <- gsub("\\s.*","",dat$Study)

  m <- meta3(logOR, v, cluster = Cluster ,data = dat)

  p1 <- metaKIN:::forest_plot(m,
                        author = "study_author", year = "study_year")

  p2 <- metaKIN:::forest_plot(m,
                        author = "study_author", year = "study_year",
                        transf = function(x) x * 2)

  testthat::expect_identical(p1$data$est * 2, p2$data$est)

})
