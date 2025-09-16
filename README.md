## libNeuralFoil

libNeuralFoil is a C++ library implementation of Peter Sharpe's [NeuralFoil](https://github.com/peterdsharpe/NeuralFoil).
NeuralFoil is a neural network model of airfoil performance trained on a large set of
[XFoil](https://web.mit.edu/drela/Public/web/xfoil/) runs.

This library version of NeuralFoil is intended for tight integration into C++ programs without introducing a Python
dependency on the resulting binary.

libNeuralFoil uses the NeuralFoil model as stored in compressed NumPy data files distributed with the NeuralFoil project.
[libnpy](https://github.com/matajoh/libnpy) is used to read these files.  libpny uses a bundled version of
[miniz](https://github.com/richgel999/miniz) to decompress the `*.npz` files.

libNeuralFoil is currently based on NeuralFoil version 0.3.2.  It may be updated in the future as NeuralFoil continues
to be developed.

## Compiling libNeuralFoil

libNeuralFoil uses CMake as a build system.  An example project using CMake is included in the `example` folder.
libNeuralFoil can be added to a CMake project via `FetchContent`.

```
include( FetchContent )

FetchContent_Declare( libneuralfoil
        GIT_REPOSITORY https://github.com/ramcdona/libNeuralFoil
        GIT_TAG origin/main
)

FetchContent_MakeAvailable( libneuralfoil )

target_link_libraries( mytarget PUBLIC
        neuralfoil
)
```

The libNeuralFoil CMake project automatically takes care of building a bundled copy of libnpy and its bundled miniz.
If you do not want to use CMake for some reason, libNeuralFoil is just two files (one source and one header) that can
easily be included in your project.  However, the libpny / miniz dependency is slightly more complex and will still
need to be built and linked into your project.

The `neuralfoil` CMake target includes a property `RESOURCE_DIR` that identifies a directory that contains
the NeuralFoil data files.  You can add a CMake command to copy the data files to your build directory, or into
your program's installation package.

```
get_target_property( neuralfoil_RESOURCE_DIR neuralfoil RESOURCE_DIR )
```

## Contact

Questions about libNeuralFoil are best directed to the [OpenVSP Google Group](https://groups.google.com/g/openvsp).

## License

Like NeuralFoil, libNeuralFoil is licensed under [the MIT license](LICENSE.txt).
