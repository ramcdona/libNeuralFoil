//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#ifndef NEURALFOIL_TEST_H
#define NEURALFOIL_TEST_H

#include "../include/neuralfoil.h"

namespace nf
{

class neuralfoil_test : public neuralfoil
{
public:

    bool test() const;
};

} // namespace nf

#endif //NEURALFOIL_TEST_H