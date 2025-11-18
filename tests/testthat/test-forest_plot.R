test_that("transf works in forest plot", {
  require(metaKIN)
  dat <- metaSEM::Bornmann07

  dat$study_year <- stringr::str_extract(dat$Study, "[0-9]{1,4}")
  dat$study_author <- gsub("\\s.*", "", dat$Study)

  m <- meta3L(logOR, v, cluster = Cluster, data = dat)

  p1 <- metaKIN:::forest_plot(m, author = "study_author", year = "study_year")

  p2 <- metaKIN:::forest_plot(
    m,
    author = "study_author",
    year = "study_year",
    transf = function(x) x * 2
  )

  testthat::expect_identical(p1$data$est * 2, p2$data$est)
})

test_that("Incorrect transf argument throws error", {
  require(metaKIN)
  dat <- metaSEM::Bornmann07[1:5, ]

  dat$study_year <- stringr::str_extract(dat$Study, "[0-9]{1,4}")
  dat$study_author <- gsub("\\s.*", "", dat$Study)

  m <- meta3L(
    logOR,
    v,
    cluster = Cluster,
    data = dat
  )

  testthat::expect_error(
    metaKIN:::forest_plot(
      m,
      author = "study_author",
      year = "study_year",
      transform = function(x) x * 2
    )
  )
})

test_that("Moderators in forest_plot works", {
  require(metaKIN)
  dat <- metaSEM::Bornmann07
  dat$study_year <- stringr::str_extract(dat$Study, "[0-9]{1,4}")
  dat$study_author <- gsub("\\s.*", "", dat$Study)

  m <- meta3L(
    logOR,
    v,
    cluster = Cluster,
    data = dat
  )
  m_mod <- m |>
    metaKIN::moderate(Country = ~ Country - 1)

  testthat::expect_no_error({
    plot_i <- metaKIN:::forest_plot(
      m_mod,
      "Country",
      author = "study_author",
      year = "study_year"
    )
  })
  testthat::expect_true("Country" %in% plot_i$data$moderation)
  dat
})
