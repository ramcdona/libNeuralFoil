//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include <cstdio>

#include "neuralfoil_test.h"
#include "test_cst.h"


int main()
{
    nf::neuralfoil_test nf;
    nf.load( "nn-xlarge.npz", "scaled_input_distribution.npz" );

    bool pass = true;

    pass = nf.test();

    if ( !pass )
    {
        printf( "Neuralfoil test failed!\n" );
    }

    pass = test_cst();

    if ( !pass )
    {
        printf( "CST test failed!\n" );
    }

    return pass ? 0 : 1;
}
