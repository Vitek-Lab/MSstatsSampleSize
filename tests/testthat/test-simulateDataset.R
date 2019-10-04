context("simulateDataset")


test_that("simulateDataset works", {

    # missing last row in data
    expect_error(simulateDataset(data = MSstatsSampleSize::OV_SRM_train[, -nrow(MSstatsSampleSize::OV_SRM_train)],
                                 annotation = MSstatsSampleSize::OV_SRM_train_annotation))

    # missing 1st row in annotation
    expect_error(simulateDataset(data = MSstatsSampleSize::OV_SRM_train,
                                 annotation = MSstatsSampleSize::OV_SRM_train_annotation[, -1]))

})
