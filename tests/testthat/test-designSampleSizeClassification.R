context("designSampleSizeClassification")


test_that("designSampleSizeClassification works", {

    output <- designSampleSizeClassification(simulations = MSstatsSampleSize::simulated_datasets)
    expect_equal(names(output), names(MSstatsSampleSize::classification_results))

})

test_that("designSampleSizeClassification handle missing information", {

    tmp <- simulated_datasets
    tmp[[1]] <- NULL

    expect_error(designSampleSizeClassification(data = tmp))

})


