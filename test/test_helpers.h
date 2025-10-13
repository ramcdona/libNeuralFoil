//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include <vector>

bool compare( double x, double xref, bool print = false );

bool compare( const std::vector < double > & x, const std::vector < double > & xref, bool print = false );

bool compare( const std::vector < std::vector < double > > &X, const std::vector < std::vector < double > > &Xref, bool print = false );

#endif //TEST_HELPERS_H
