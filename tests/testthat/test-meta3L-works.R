test_that("Moderate works with meta3L", {

  dat <- metaSEM::Bornmann07

  meta3L(logOR, v, cluster = Cluster, data = dat)  |>
    metaKIN::moderate(Discipline = ~ Discipline) |>
    testthat::expect_no_error()

})
