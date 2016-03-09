# Compiling TESSE

* **Setup**: Create a directory for the entire project and put its location to an environment variable for convenience. Also create an installation directory for all the software inside and we will configure them accordingly. We will also add the binary part of the installation directory to our PATH.
```shell
  export TESSE_DIR=*my development directory*
  export PATH=${PATH}:${TESSE_DIR}/install/bin
```

* **Compile PaRSEC**: 

  * Obtain the latest and greatest PaRSEC source from the main bitbucket repo.
```shell
  git clone git@bitbucket.org:icldistcomp/parsec.git
```

  * Let’s use a VPATH compilation to keep the sources clean.
```shell
  cd parsec
  mkdir build
  cd build
```

  * Configure PaRSEC. This requires to have -fPIC in the CFLAGS, install the development headers, and then finally set the PKG_CONFIG_PATH to the right location. If you want a debug build add “-DCMAKE_BUILD_TYPE=Debug“ before the “../“.
```shell
cmake -G 'Unix Makefiles' -DBUILD_SHARED_LIBS=ON —DDAGUE_WITH_DEVEL_HEADERS=ON DHWLOC_DIR=/Users/bosilca/opt -DDAGUE_DIST_WITH_MPI=OFF -DDAGUE_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX:PATH=${TESSE_DIR}/install -DBUILD_DPLASMA=OFF ../
```
PaRSEC will complain about missing BLAS libraries, but for the purpose of TESSE we don’t need them. Take a careful look at the output, to make sure HWLOC has been correctly found.
  * Install PaRSEC and add its pkgconfig location to the PKG_CONFIG_PATH environment variable:
```shell
make install
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${TESSE_DIR}/install/lib/pkgconfig
```

* **Compile MADNESS**.
  * Clone the TESSE repo, and switch to the correct branch.
```shell
git clone -b madness_over_parsec git@github.com:TESSEorg/madness.git
cd madness
```

  * Configure. There are 2 ways to do that: using autotools or CMake.
  
    * Using autotools (building in the source):
```shell
sh ./autogen.sh
./configure --without-tbb  —with-parsec=${TESSE_DIR}
```
To enable debugging add these to the configure command: ```--disable-optimization --enable-debugging```. See ```./configure --help``` for additional options.

    * Using CMake (with a separate build directory):
```shell
mkdir build; cd build
cmake -Wno-dev -D ENABLE_TBB=OFF ..
```
To enable debugging use add the following command-line option to ```cmake```: ```-DCMAKE_BUILD_TYPE=Debug```. It is recommended to collect all CMake variable definitions in a toolchain file, e.g. see ```madness/cmake/toolchains/generic-mkl-tbb.cmake```. You can specify the toolchain file to use by providing a command-line option to ```cmake``` as ```-DCMAKE_TOOLCHAIN_FILE=<toolchain-file-path>```.

  * Build.
```
make
```
