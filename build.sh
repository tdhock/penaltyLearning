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
rm $RELEASE/tests/testthat/test-peaks.R
rm $RELEASE/tests/testthat/test-ROChange-jointseg.R
rm -rf $RELEASE/vignettes

echo Building $RELEASE
RCMD="R --vanilla CMD"
$RCMD build $RELEASE | tee build.out
PKG_TGZ=$(grep building build.out|sed "s/.*\($PKG.*.tar.gz\).*/\1/")

# https://www.stats.ox.ac.uk/pub/bdr/noRemap/README.txt
# https://cloud.r-project.org/doc/manuals/R-exts.html#The-R-API
export _R_CXX_USE_NO_REMAP_=true

echo Installing $PKG_TGZ
$RCMD INSTALL $PKG_TGZ

echo Checking $PKG_TGZ
$RCMD check --as-cran $PKG_TGZ

echo Checking without any Suggests
R -e "if('check_without_suggests' %in% ls())check_without_suggests('$PKG_TGZ')"
