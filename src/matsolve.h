//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#ifndef MATSOLVE_H
#define MATSOLVE_H

#include <vector>

namespace nf
{

void form_normal_equations( std::vector < std::vector < double > > &ATA, std::vector < double > &ATb, const std::vector < std::vector < double > > &A, const std::vector < double > &b );
bool solve_linear_system( std::vector < double > &x, std::vector < std::vector < double > > &A, std::vector < double > &b );
void solve_least_squares( std::vector < double > &x, const std::vector < std::vector < double > > &A, const std::vector < double > &b );

void quadratic_critical_point( double & xcrit, double & ycrit, double x0, double y0, double y1, double x2, double y2 );
double quadratic_eval( double x, double x0, double y0, double y1, double x2, double y2 );

} // namespace nf
#endif //MATSOLVE_H
