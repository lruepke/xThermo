cd gsl
./autogen.sh
mkdir build
cd build
../configure --prefix=${PWD}/../../install/gsl_`uname -s`_`arch`
# add -fPIC to cblas compile flag
sed -i '/mode=compile/a\        -fPIC \\' cblas/Makefile
# make
make -j 8
make install