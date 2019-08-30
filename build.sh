#!/bin/bash
set -o errexit
PREGEX="^Package: "
PKG=$(grep $PREGEX DESCRIPTION|sed "s/$PREGEX//")
echo Package from DESCRIPTION: $PKG
R -e "if(require(inlinedocs))package.skeleton.dx('.')"
cd ..

RELEASE=$PKG-release
echo Copying $PKG to $RELEASE
rm -rf $RELEASE
cp -r $PKG $RELEASE

echo Editing $RELEASE for CRAN submission
## remove long running tests and the packages that they use in
## Suggests.
grep -v PeakSegDP penaltyLearning/DESCRIPTION > $RELEASE/DESCRIPTION
rm $RELEASE/tests/testthat/test-IntervalRegression.R
rm $RELEASE/tests/testthat/test-labelError.R
rm $RELEASE/tests/testthat/test-modelSelection.R
rm $RELEASE/tests/testthat/test-peaks.R
rm $RELEASE/tests/testthat/test-ROChange.R
rm -rf $RELEASE/vignettes

echo Building $RELEASE
RCMD="R --vanilla CMD"
$RCMD build $RELEASE | tee build.out
PKG_TGZ=$(grep building build.out|sed "s/.*\($PKG.*.tar.gz\).*/\1/")

echo Installing $PKG_TGZ
$RCMD INSTALL $PKG_TGZ

echo Checking $PKG_TGZ
$RCMD check --as-cran $PKG_TGZ
