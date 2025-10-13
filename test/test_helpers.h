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

std::vector < std::vector < double > > generate_NACA4( double m, double p, double t, int n_points );

void write_matlab_plot( const std::vector < std::vector < double > > &orig, const std::vector < std::vector < double > > &recon, const std::string &fname );

#endif //TEST_HELPERS_H
