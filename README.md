
SMART-SED is a program designed to model hydrological and erosion processes. 
It was developed using C++11 and R.

For the R part you need to install the following libraries,
raster, gstat, compositions, dissever, fields, soiltexture, viridis, psych, latex2exp

You need to install the gdal binaries.

The C++11 part is self contained, i.e. the github repository is self-contained, it relies on the Eigen library 
to solve the linear suystem arising in the discretization of the superficial run-off dynamics.


From vim type this to format to a given file formatting style
:!clang-format -i %


- Release build (optimized, default)
cmake -S . -B build
cmake --build build -j

- Debug build
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug -j

- Enable CUDA
cmake -S . -B build-cuda -DENABLE_CUDA=ON
cmake --build build-cuda -j

- Debug + CUDA
cmake -S . -B build-cuda-debug \
      -DENABLE_CUDA=ON \
      -DCMAKE_BUILD_TYPE=Debug
cmake --build build-cuda-debug -j


