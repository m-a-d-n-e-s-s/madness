# Compiling TESSE

* **Setup**: Create a directory for the entire project and store its location into an environment variable for convenience. Create an instal directory for all the TESSE-related software. We will configure them accordingly, to minimize the impact on the different standard path. First thing, we will add the binary part of the install directory to our PATH, to get first hand access to the software. In case you compile everything as DSO, you might also need to update the LD_LIBRARY_PATH.
  ```shell
  export TESSE_DIR=*my development directory*
  export PATH=${PATH}:${TESSE_DIR}/install/bin
  ```

* **Compile PaRSEC**: 

  * Obtain the latest version of the PaRSEC source directly from the main bitbucket repo.

    ```shell
    git clone git@bitbucket.org:icldistcomp/parsec.git
    ```

  * I suggest to keep the sources clean and use a VPATH compilation, but you are not required to do so. In case you choose to follow my advice:

    ```shell
    cd parsec
    mkdir build
    cd build
    ```

  * Configuring PaRSEC should go smoothly. In the context of TESSE, it is necessary to produce shared libraries and to force the installation of the PaRSEC development headers. As a hint, if you want a debug build add ```-DCMAKE_BUILD_TYPE=Debug`` before the ```../```.

    ```shell
    cmake -G 'Unix Makefiles' -DBUILD_SHARED_LIBS=ON —DDAGUE_WITH_DEVEL_HEADERS=ON -DHWLOC_DIR=<hwloc-install-prefix> -DDAGUE_DIST_WITH_MPI=ON -DDAGUE_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX:PATH=${TESSE_DIR}/install -DBUILD_DPLASMA=OFF ../
    ```
  
    PaRSEC will complain about missing BLAS libraries, but for the purpose of TESSE we don’t need them. Take a careful look at the output, to make sure the important pieces are indeed correctly found. Look for HWLOC and MPI at a minimum, CUDA if you need GPU support. It is also important to have support for atomic operations (especially 128 bits, such a double CAS to prevent ABA issues).
  
  * Install PaRSEC and update your environment variables. Add PaRSEC pkgconfig location to the PKG_CONFIG_PATH environment variable so that MADNESS will find it automatically:

    ```shell
    make install
    export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${TESSE_DIR}/install/lib/pkgconfig
    ```
  * Congratulation! You just got yourself a clean installation of PaRSEC.
  

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
    
      To enable debugging add the following to the configure command: ```--disable-optimization --enable-debugging```. See ```./configure --help``` for additional options.

    * Using CMake (with a separate build directory):
    
      ```shell
      mkdir build; cd build
      cmake -Wno-dev -D ENABLE_TBB=OFF ..
      ```
    
      To enable debugging  add the following command-line option to ```cmake```: ```-DCMAKE_BUILD_TYPE=Debug```. It is recommended to collect all CMake variable definitions in a toolchain file, e.g. see ```madness/cmake/toolchains/generic-mkl-tbb.cmake```. You can specify the toolchain file to use by providing a command-line option to ```cmake``` as ```-DCMAKE_TOOLCHAIN_FILE=<toolchain-file-path>```.

  * Build.

    ```
    make
    ```
 * You are now ready to hack into the TESSE source code!
  
