//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include <cstdio>

#include "test_neuralfoil.h"
#include "test_cst.h"


int main()
{
    bool pass = true;

    pass = test_neuralfoil();

    if ( !pass )
    {
        printf( "Neuralfoil test failed!\n" );
    }

    pass = test_cst();

    if ( !pass )
    {
        printf( "CST test failed!\n" );
    }

    // test_cst_converge();

    return pass ? 0 : 1;
}
