context("utility")

test_that("run utility should succeed", {
    data("yeastTimeSeries")
    expect_false(isBooleanTypeTimeseriesData(yeastTimeSeries))
    ## discrete the result
    timeseries_data <- BoolNet::binarizeTimeSeries(yeastTimeSeries, method = "kmean")$binarizedMeasurements
    expect_true(isBooleanTypeTimeseriesData(timeseries_data))
    
    expect_true(is.null(checkProbabilityTypeData(0.6)))
    expect_error(checkProbabilityTypeData(2), "*not a type of probability*")
    expect_error(checkProbabilityTypeData(-0.2), "*not a type of probability*")
    expect_error(checkProbabilityTypeData(-2), "*not a type of probability*")
    expect_error(checkProbabilityTypeData("X"), "*not a type of probability*")
    expect_true(is.null(checkNumeric(1)))
    expect_error(checkNumeric("X"), "*type of numeric*")
    expect_true(is.null(CheckRightTypeTimeseriesData(list(yeastTimeSeries))))
    expect_error(CheckRightTypeTimeseriesData(yeastTimeSeries), "*timeseries_data must be LIST*")
    expect_error(CheckRightTypeTimeseriesData(list(1, 2)), "*must be a matrix*")
    
    t1 <- expect_error(dividedVectorIntoSmallgroups(c(1, 2, 3, 4, 5, 6, 7, 8), 2), NA)
    expect_true(length(t1$clusters) == 4)
    
    t1 <- list(yeastTimeSeries, yeastTimeSeries, yeastTimeSeries)
    t2 <- row.names(yeastTimeSeries)
    t3 <- expect_error(dividedVectorIntoSmallgroups(t2, 2), NA)
    t4 <- expect_error(getRelatedGeneTimeseries(t1, t3$clusters[[1]]), NA)
    expect_true(all(row.names(t4[[1]]) == t3$clusters[[1]]))
    t5 <- expect_error(generateAllCombinationBinary(t3$clusters[[1]], 1, 0), NA)
    expect_true(length(t5) == 4)
    t6 <- expect_error(randomGenerateBinary(t3$clusters[[1]], 4), NA)
})
