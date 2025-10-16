//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#ifndef NEURALFOIL_TEST_H
#define NEURALFOIL_TEST_H

#include "../include/neuralfoil.h"

namespace nf
{

class neuralfoil_public : public neuralfoil
{
public:
    using neuralfoil::sigmoid;
    using neuralfoil::swish;
    using neuralfoil::multiply;
    using neuralfoil::squared_mahalanobis_distance;
    using neuralfoil::net;
    using neuralfoil::evaluate;
    using neuralfoil::inputs;
    using neuralfoil::flipx;
    using neuralfoil::unflipy;
    using neuralfoil::fusey;
    using neuralfoil::unpacky;
};

} // namespace nf

#endif //NEURALFOIL_TEST_H
