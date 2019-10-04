context("estimateVar")

test_that("estimateVar works", {

    output <- estimateVar(data = MSstatsSampleSize::OV_SRM_train,
                          annotation = MSstatsSampleSize::OV_SRM_train_annotation)
    expect_equal(output, MSstatsSampleSize::variance_estimation, tolerance=1e-5)

})

test_that("estimateVar handles information", {

    # missing last row in data
    expect_error(estimateVar(data = MSstatsSampleSize::OV_SRM_train[, -nrow(MSstatsSampleSize::OV_SRM_train)],
                             annotation = MSstatsSampleSize::OV_SRM_train_annotation))

    # missing 1st row in annotation
    expect_error(estimateVar(data = MSstatsSampleSize::OV_SRM_train,
                             annotation = MSstatsSampleSize::OV_SRM_train_annotation[, -1]))

})


