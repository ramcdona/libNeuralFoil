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

## Contact

Questions about libNeuralFoil are best directed to the [OpenVSP Google Group](https://groups.google.com/g/openvsp).

## License

Like NeuralFoil, libNeuralFoil is licensed under [the MIT license](LICENSE.txt).
