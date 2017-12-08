#!/bin/bash
cd ..
rm -rf penaltyLearning-release
cp -r penaltyLearning penaltyLearning-release
## remove long running tests and the packages that they use in
## Suggests.
grep -v PeakSegDP penaltyLearning/DESCRIPTION > penaltyLearning-release/DESCRIPTION
rm penaltyLearning-release/tests/testthat/test-IntervalRegression.R
rm penaltyLearning-release/tests/testthat/test-labelError.R
rm penaltyLearning-release/tests/testthat/test-modelSelection.R
rm penaltyLearning-release/tests/testthat/test-peaks.R
rm penaltyLearning-release/tests/testthat/test-ROChange.R
rm -rf penaltyLearning-release/vignettes
PKG_TGZ=$(R CMD build penaltyLearning-release|grep penaltyLearning_|sed 's/.*‘//'|sed 's/’.*//')
R CMD check --as-cran $PKG_TGZ
