#!/bin/sh
export SKFlat_BUILDDIR=$SKFlat_WD/build

echo @@@@ Prepare to build SKNanaAnalyzer in $SKNANO_BUILDDIR
rm -rf $SKFlat_BUILDDIR $SKFlat_LIB_PATH
mkdir $SKFlat_BUILDDIR && cd $SKFlat_BUILDDIR

cmake -DCMAKE_INSTALL_PREFIX=$SKFlat_WD $SKFlat_WD

echo @@@@ make -j6
make -j6

echo @@@@ install to $SKFlat_LIB_PATH
make install
