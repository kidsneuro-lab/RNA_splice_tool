library(testthat)
library(cortar)

# Test is_debug_enabled() -----------------------------------------------------

test_that("is_debug_enabled returns TRUE for valid debug paths", {
  expect_true(is_debug_enabled("/path/to/debug"))
  expect_true(is_debug_enabled("debug"))
  expect_true(is_debug_enabled("./output/debug"))
})

test_that("is_debug_enabled returns FALSE for invalid debug values", {
  expect_false(is_debug_enabled(FALSE))
  expect_false(is_debug_enabled(""))
  expect_false(is_debug_enabled(TRUE))
  expect_false(is_debug_enabled(NULL))
})

# Test get_coverage_threshold() -----------------------------------------------

test_that("get_coverage_threshold returns correct values for standard types", {
  expect_equal(get_coverage_threshold("het"), 60)
  expect_equal(get_coverage_threshold("hom"), 30)
  expect_equal(get_coverage_threshold("hemi"), 30)
  expect_equal(get_coverage_threshold(""), 0)
})

test_that("get_coverage_threshold handles custom numeric values", {
  expect_equal(get_coverage_threshold("50"), 50)
  expect_equal(get_coverage_threshold("100"), 100)
  expect_equal(get_coverage_threshold("15"), 15)
})

# Test %nin% operator ---------------------------------------------------------

test_that("%nin% operator works correctly", {
  expect_true("a" %nin% c("b", "c", "d"))
  expect_false("a" %nin% c("a", "b", "c"))
  expect_true(1 %nin% c(2, 3, 4))
  expect_false(1 %nin% c(1, 2, 3))
})

test_that("%nin% handles edge cases", {
  expect_true("x" %nin% character(0))
  expect_false(NA %nin% c(NA, "a", "b"))
  expect_true("test" %nin% NULL)
})
