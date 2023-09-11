  require(metaKIN)

  dat <- metaSEM::Bornmann07
  dat$n <- seq_len(nrow(dat))

  m0 <- meta3L(logOR, v, cluster = Cluster, data = dat)
  mod <- m0 |>
  moderate(
    Discipline = ~ Discipline - 1,
    store_var = list(
      within = ~ max(n, na.rm = TRUE),
      between = ~ format(sum(., na.rm = TRUE), big.mark = ",")
    )
  ) |> suppressWarnings()

  dat_missing <- dat
  dat_missing$Discipline[dat_missing$Cluster == 1] <- NA

  m0_missing <- meta3L(logOR, v, cluster = Cluster, data = dat_missing)
  mod_missing <- m0_missing |>
  moderate(
    Discipline = ~ Discipline - 1,
    store_var = list(
      within = ~ max(n, na.rm = TRUE),
      between = ~ format(sum(., na.rm = TRUE), big.mark = ",")
    )
  ) |> suppressWarnings()


test_that("Store var works", {

  store_var <- attr(mod$models$Baseline, "store_var")

  manual <- as.character(sum(data.table(dat)[, list(res = max(n)) , by = "Cluster"]$res))

  testthat::expect_equal(store_var, manual)

})

test_that("Store var responsive to missing moderator data", {
  manual <-
    as.character(sum(data.table(dat_missing)[Cluster != 1, list(res = max(n)), by = "Cluster"]$res))

  testthat::expect_equal(attr(mod_missing$models$Discipline, "store_var"),
                         manual)
})

test_that("Store var does not appear uninvited into format_nicely", {

  testthat::expect_false("store_var" %in% names(format_nicely(mod_missing)))

})


test_that("Store var does appear invited into format_nicely", {
  testthat::expect_true("TESTME" %in% names(format_nicely(mod_missing, store_var = "TESTME")))
})
