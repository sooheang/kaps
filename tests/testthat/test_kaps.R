require(kaps)

context("kaps")

data(toy)
f <- Surv(time, status) ~ meta

# Fit kaps algorithm without cross-validation.
# It means the step to finding optimal K is not entered.
fit1 <- kaps(f, data = toy, K = 3)

test_that('Unit test for a main wrapper function, kaps()', {
  fit1 = kaps(f, data = toy, K = 3)  
  preds = predict(fit1)
})
    

