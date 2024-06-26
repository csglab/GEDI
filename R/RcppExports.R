# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

eigenMatProduct <- function(A, B) {
    .Call('_GEDI_eigenMatProduct', PACKAGE = 'GEDI', A, B)
}

eigenMatTcrossprod <- function(A, B) {
    .Call('_GEDI_eigenMatTcrossprod', PACKAGE = 'GEDI', A, B)
}

eigenMatCrossprod <- function(A, B) {
    .Call('_GEDI_eigenMatCrossprod', PACKAGE = 'GEDI', A, B)
}

eigenMatVecProduct <- function(A, b) {
    .Call('_GEDI_eigenMatVecProduct', PACKAGE = 'GEDI', A, b)
}

eigenVecMatProduct <- function(a, B) {
    .Call('_GEDI_eigenVecMatProduct', PACKAGE = 'GEDI', a, B)
}

eigenVecVecProduct <- function(a, b) {
    .Call('_GEDI_eigenVecVecProduct', PACKAGE = 'GEDI', a, b)
}

eigenRowL2 <- function(A) {
    .Call('_GEDI_eigenRowL2', PACKAGE = 'GEDI', A)
}

eigenColL2 <- function(A) {
    .Call('_GEDI_eigenColL2', PACKAGE = 'GEDI', A)
}

eigenColL1 <- function(A) {
    .Call('_GEDI_eigenColL1', PACKAGE = 'GEDI', A)
}

eigenRowMult <- function(A, b) {
    .Call('_GEDI_eigenRowMult', PACKAGE = 'GEDI', A, b)
}

solveBi <- function(Yi, D, Z, Qi, diagK, si, o, oi) {
    .Call('_GEDI_solveBi', PACKAGE = 'GEDI', Yi, D, Z, Qi, diagK, si, o, oi)
}

Yi_resZ <- function(Yi, QiDBi, si, o, oi) {
    .Call('_GEDI_Yi_resZ', PACKAGE = 'GEDI', Yi, QiDBi, si, o, oi)
}

solveZ_noC <- function(Y_res, D, B, diagK, lambda) {
    .Call('_GEDI_solveZ_noC', PACKAGE = 'GEDI', Y_res, D, B, diagK, lambda)
}

solveZ_wC <- function(Y_res, D, B, C, A, diagK, lambda) {
    .Call('_GEDI_solveZ_wC', PACKAGE = 'GEDI', Y_res, D, B, C, A, diagK, lambda)
}

solveU <- function(Y_res, S, B) {
    .Call('_GEDI_solveU', PACKAGE = 'GEDI', Y_res, S, B)
}

solveQi_noH <- function(Yi, ZDBi, si, o, oi, D, Bi, diagK, lambda) {
    .Call('_GEDI_solveQi_noH', PACKAGE = 'GEDI', Yi, ZDBi, si, o, oi, D, Bi, diagK, lambda)
}

solveQi_wH <- function(Yi, ZDBi, si, o, oi, D, Bi, Qi_hat, diagK, lambda) {
    .Call('_GEDI_solveQi_wH', PACKAGE = 'GEDI', Yi, ZDBi, si, o, oi, D, Bi, Qi_hat, diagK, lambda)
}

solveOi_noH <- function(Yi, ZDBi, QiDBi, ni, si, o, Ni, lambda) {
    .Call('_GEDI_solveOi_noH', PACKAGE = 'GEDI', Yi, ZDBi, QiDBi, ni, si, o, Ni, lambda)
}

solveOi_wH <- function(Yi, ZDBi, QiDBi, ni, si, o, oi_hat, Ni, lambda) {
    .Call('_GEDI_solveOi_wH', PACKAGE = 'GEDI', Yi, ZDBi, QiDBi, ni, si, o, oi_hat, Ni, lambda)
}

solveSi <- function(Yi, ZDBi, QiDBi, j, o, oi, si_0, J, lambda) {
    .Call('_GEDI_solveSi', PACKAGE = 'GEDI', Yi, ZDBi, QiDBi, j, o, oi, si_0, J, lambda)
}

Yi_resO_rowSum <- function(Yi, ZDBi, QiDBi, ni, si, oi) {
    .Call('_GEDI_Yi_resO_rowSum', PACKAGE = 'GEDI', Yi, ZDBi, QiDBi, ni, si, oi)
}

solveYi <- function(Mi, Yi, ZDBi, QiDBi, si, o, oi, sigma2) {
    .Call('_GEDI_solveYi', PACKAGE = 'GEDI', Mi, Yi, ZDBi, QiDBi, si, o, oi, sigma2)
}

solveYi_paired <- function(M1i, M2i, Yi, ZDBi, QiDBi, si, o, oi, sigma2) {
    .Call('_GEDI_solveYi_paired', PACKAGE = 'GEDI', M1i, M2i, Yi, ZDBi, QiDBi, si, o, oi, sigma2)
}

vecL2_noPrior <- function(x) {
    .Call('_GEDI_vecL2_noPrior', PACKAGE = 'GEDI', x)
}

vecL2_wPrior <- function(x, prior) {
    .Call('_GEDI_vecL2_wPrior', PACKAGE = 'GEDI', x, prior)
}

matL2_noPrior <- function(X) {
    .Call('_GEDI_matL2_noPrior', PACKAGE = 'GEDI', X)
}

matL2_wPrior <- function(X, A, B) {
    .Call('_GEDI_matL2_wPrior', PACKAGE = 'GEDI', X, A, B)
}

Yi_SSE_fixed <- function(Yi, ZDBi, QiDBi, si, o, oi) {
    .Call('_GEDI_Yi_SSE_fixed', PACKAGE = 'GEDI', Yi, ZDBi, QiDBi, si, o, oi)
}

Yi_SSE_M <- function(Yi, ZDBi, QiDBi, si, o, oi, sigma2) {
    .Call('_GEDI_Yi_SSE_M', PACKAGE = 'GEDI', Yi, ZDBi, QiDBi, si, o, oi, sigma2)
}

Yi_SSE_M_paired <- function(Yi, M1i, M2i, ZDBi, QiDBi, si, o, oi, sigma2) {
    .Call('_GEDI_Yi_SSE_M_paired', PACKAGE = 'GEDI', Yi, M1i, M2i, ZDBi, QiDBi, si, o, oi, sigma2)
}

Yi_var <- function(Yi, sigma2) {
    .Call('_GEDI_Yi_var', PACKAGE = 'GEDI', Yi, sigma2)
}

Yi_var_paired <- function(Yi, M1i, M2i, sigma2) {
    .Call('_GEDI_Yi_var_paired', PACKAGE = 'GEDI', Yi, M1i, M2i, sigma2)
}

predict_Yhat <- function(ZDBi, QiDBi, si, o, oi) {
    .Call('_GEDI_predict_Yhat', PACKAGE = 'GEDI', ZDBi, QiDBi, si, o, oi)
}

matRMSD <- function(A, B) {
    .Call('_GEDI_matRMSD', PACKAGE = 'GEDI', A, B)
}

vecRMSD <- function(A, B) {
    .Call('_GEDI_vecRMSD', PACKAGE = 'GEDI', A, B)
}

