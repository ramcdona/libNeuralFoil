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

##### All project variables:

   - `LNF_BUILD_TESTS` -- Set this variable to enable building some low-level tests.

## Usage

libNeuralFoil provides a single class `nf::neuralfoil` in the `nf` namespace.  `nf::neuralfoil` provides a simple
public API with methods to `load` a NeuralFoil model from a file and to `evaluate` that model.

Loading a model from file is substantially more expensive than a single evaluation of the model.  Consequently,
if your program is going to evaluate the NeuralFoil model repeatedly (the expected use case), you should take care
to load the model once and keep the `nf::neuralfoil` class instance in scope for the duration of your program.

The example use of libNeuralFoil from the `example` folder is included below.  In this situation, the majority of
code is setting up the sample inputs for the model and variables to hold the model outputs.  Instancing the class,
loading the model, and evaluating the model each take a single line of code.

```
#include "neuralfoil.h"
#include <vector>

int main()
{
    // CST coefficients for NACA 4412 airfoil.
    std::vector < double > CST_up = {
        0.18109497123192536, 0.2126841855378795, 0.28098503347171055,
        0.24864887442714956, 0.2402814006497832, 0.2726284343621469,
        0.2577647399772107, 0.27817638080391327
    };

    std::vector < double > CST_low = {
        -0.16965145937914555, -0.09364138259246524, -0.06345895875250755,
        -0.006796600807436377, -0.09024470493931744, 0.020818449122134496,
        -0.03575216051938629, -0.0022362344439678192
    };

    double CST_le = 0.10647339061374347;
    double CST_te = 0.0025720113171501307;

    double alpha = 5.159747295611583;
    double Re = 1000302.1671445925;

    double n_crit = 9.0;
    double xtr_upper = 1.0;
    double xtr_lower = 1.0;

    // Outputs
    double analysis_confidence;
    double CL;
    double CD;
    double CM;
    double Top_Xtr;
    double Bot_Xtr;
    std::vector < double > upper_bl_ue_over_vinf;
    std::vector < double > lower_bl_ue_over_vinf;
    std::vector < double > upper_theta;
    std::vector < double > lower_theta;
    std::vector < double > upper_H;
    std::vector < double > lower_H;

    nf::neuralfoil nf;
    nf.load( "nn-xlarge.npz", "scaled_input_distribution.npz" );

    nf.evaluate( analysis_confidence,
                 CL,
                 CD,
                 CM,
                 Top_Xtr,
                 Bot_Xtr,
                 upper_bl_ue_over_vinf,
                 lower_bl_ue_over_vinf,
                 upper_theta,
                 lower_theta,
                 upper_H,
                 lower_H,
                 CST_up,
                 CST_low,
                 CST_le,
                 CST_te,
                 alpha * M_PI / 180.0,
                 Re,
                 n_crit,
                 xtr_upper,
                 xtr_lower );

    return 0;
}
```

## Contact

Questions about libNeuralFoil are best directed to the [OpenVSP Google Group](https://groups.google.com/g/openvsp).

## License

Like NeuralFoil, libNeuralFoil is licensed under [the MIT license](LICENSE.txt).
