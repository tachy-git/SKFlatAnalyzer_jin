#!/bin/sh
# clear the lhapdf directory
rm -rf $SKFlat_WD/external/lhapdf

# download and install lhapdf
cd $SKFlat_WD/external
wget https://lhapdf.hepforge.org/downloads/\?f\=LHAPDF-6.5.4.tar.gz -O LHAPDF-6.5.4.tar.gz
tar xf LHAPDF-6.5.4.tar.gz
cd LHAPDF-6.5.4/config
wget https://git.savannah.gnu.org/gitweb/\?p\=config.git\;a\=blob_plain\;f\=config.sub\;hb\=HEAD -O config.sub
wget https://git.savannah.gnu.org/gitweb/\?p\=config.git\;a\=blob_plain\;f\=config.guess\;hb\=HEAD -O config.guess
cd ..
ARCH=`arch`
if [ $ARCH == "arm64" ]; then
    ./configure --prefix=$SKFlat_WD/external/lhapdf CXX=clang++
else
    ./configure --prefix=$SKFlat_WD/external/lhapdf
fi
make -j4 && make install
cd $SKFlat_WD/external
rm -rf LHAPDF*

# download CMS central PDF sets
mkdir -p $SKFlat_WD/external/lhapdf/data
cd $SKFlat_WD/external/lhapdf/data
wget https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_hessian_pdfas.tar.gz
wget https://lhapdfsets.web.cern.ch/current/NNPDF31_nlo_hessian_pdfas.tar.gz
tar xf NNPDF31_nnlo_hessian_pdfas.tar.gz
tar xf NNPDF31_nlo_hessian_pdfas.tar.gz
cd $SKFlat_WD
