# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

SADMVN <- function(M, C) {
    .Call('lxspline_SADMVN', PACKAGE = 'lxspline', M, C)
}

BivNProb <- function(mean, cv) {
    .Call('lxspline_BivNProb', PACKAGE = 'lxspline', mean, cv)
}

TriNProb <- function(mean, cv) {
    .Call('lxspline_TriNProb', PACKAGE = 'lxspline', mean, cv)
}

Tnorm <- function(mean, sd) {
    .Call('lxspline_Tnorm', PACKAGE = 'lxspline', mean, sd)
}

sampleBetas <- function(ttY, ttX, tbetas, LAM, intLAM, p, tau) {
    .Call('lxspline_sampleBetas', PACKAGE = 'lxspline', ttY, ttX, tbetas, LAM, intLAM, p, tau)
}

sinsertBeta <- function(tY, Xidx, ctau, tp, lam) {
    .Call('lxspline_sinsertBeta', PACKAGE = 'lxspline', tY, Xidx, ctau, tp, lam)
}

sdeleteBeta <- function(tY, Xidx, ctau, tp, lam) {
    .Call('lxspline_sdeleteBeta', PACKAGE = 'lxspline', tY, Xidx, ctau, tp, lam)
}

shapesplineInsert <- function(k, tp, txi, tdeg, tCBX, tpos) {
    .Call('lxspline_shapesplineInsert', PACKAGE = 'lxspline', k, tp, txi, tdeg, tCBX, tpos)
}

shapesplineDelete <- function(k, tp, txi, tdeg, tCBX, tpos) {
    .Call('lxspline_shapesplineDelete', PACKAGE = 'lxspline', k, tp, txi, tdeg, tCBX, tpos)
}

rtmvn <- function(tMean, tVar) {
    .Call('lxspline_rtmvn', PACKAGE = 'lxspline', tMean, tVar)
}

qcopy <- function(a, b, c, d) {
    .Call('lxspline_qcopy', PACKAGE = 'lxspline', a, b, c, d)
}

