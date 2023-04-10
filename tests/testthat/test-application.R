context("application")

test_that("run application should succeed", {
    data("yeastTimeSeries")
    expect_error(generateFBMNetwork(yeastTimeSeries), NA)
    expect_error(generateFBMNetwork(yeastTimeSeries, verbose = TRUE), NA)
    expect_error(generateFBMNetwork(yeastTimeSeries, method = "kmeans"), NA)
    expect_error(generateFBMNetwork(yeastTimeSeries, method = "edgeDetector"), NA)
    expect_error(generateFBMNetwork(yeastTimeSeries, method = "scanStatistic"), NA)
    expect_error(generateFBMNetwork(yeastTimeSeries, network_only = FALSE), NA)

    expect_warning(generateFBMNetwork(yeastTimeSeries), NA)
    expect_warning(generateFBMNetwork(yeastTimeSeries, verbose = TRUE), NA)
    expect_warning(generateFBMNetwork(yeastTimeSeries, method = "kmeans"), NA)
    expect_warning(generateFBMNetwork(yeastTimeSeries, method = "edgeDetector"), NA)
    ##expect_warning(generateFBMNetwork(yeastTimeSeries, method = "scanStatistic"), "*show a uniform behaviour and should*")
    expect_warning(generateFBMNetwork(yeastTimeSeries, network_only = FALSE), NA)
})
