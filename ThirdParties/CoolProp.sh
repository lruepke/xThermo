cd coolprop
git apply ../CoolProp_patches/381c8535.patch
system_name=`uname -s`

for arc in arm64 x86_64
do
    mkdir build_${system_name}_${arc}
    cd build_${system_name}_${arc}
    cmake -DCMAKE_OSX_ARCHITECTURES=${arc} -DCOOLPROP_INSTALL_PREFIX=${PWD}/../../install/CoolProp_${system_name}_${arc} -DCOOLPROP_STATIC_LIBRARY=ON ..
    make -j 8
    make install
    rm -rf *
    cmake -DCMAKE_OSX_ARCHITECTURES=${arc} -DCOOLPROP_INSTALL_PREFIX=${PWD}/../../install/CoolProp_${system_name}_${arc} -DCOOLPROP_SHARED_LIBRARY=ON ..
    make -j 8
    make install
    cd ..
done
