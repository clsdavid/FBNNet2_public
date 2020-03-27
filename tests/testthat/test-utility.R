context("utility")

test_that("run utility should succeed", {
  data('yeastTimeSeries')
  expect_false(isBooleanTypeTimeseriesData(yeastTimeSeries))
  ## discrete the result
  timeseries_data <- BoolNet::binarizeTimeSeries(yeastTimeSeries,
                                                 method = "kmean")$binarizedMeasurements
  expect_true(isBooleanTypeTimeseriesData(timeseries_data))
  
  expect_true(is.null(checkProbabilityTypeData(0.6)))
  expect_error(checkProbabilityTypeData(2), "*not a type of probability*")
  expect_error(checkProbabilityTypeData(-0.2), "*not a type of probability*")
  expect_error(checkProbabilityTypeData(-2), "*not a type of probability*")
  expect_error(checkProbabilityTypeData("X"), "*not a type of probability*")
  expect_true(is.null(checkNumeric(1)))
  expect_error(checkNumeric("X"), "*type of numeric*")
  expect_true(is.null(CheckRightTypeTimeseriesData(list(yeastTimeSeries))))
  expect_error(CheckRightTypeTimeseriesData(yeastTimeSeries), 
               "*timeseries_data must be LIST*")
  expect_error(CheckRightTypeTimeseriesData(list(1,2)), 
               "*must be a matrix*")
})
