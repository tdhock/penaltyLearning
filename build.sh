#!/bin/bash
cd ..
rm -rf penaltyLearning-release
cp -r penaltyLearning penaltyLearning-release
grep -v PeakSegDP penaltyLearning/DESCRIPTION > penaltyLearning-release/DESCRIPTION
pushd penaltyLearning-release/tests/testthat
rm test-IntervalRegression.R test-labelError.R test-modelSelection.R test-peaks.R test-ROChange.R
popd
PKG_TGZ=$(R CMD build penaltyLearning-release|grep penaltyLearning_|sed 's/.*‘//'|sed 's/’.*//')
R CMD check --as-cran $PKG_TGZ
