//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include "matsolve.h"

#include <cmath>
#include <vector>

#include <algorithm>  // for std::swap

namespace nf
{

// Form normal equations for least squares
void form_normal_equations( std::vector < std::vector < double > > &ATA, std::vector < double > &ATb, const std::vector < std::vector < double > > &A, const std::vector < double > &b )
{
    int m = ( int )A.size();
    if ( m == 0 )
    {
        ATA.clear();
        ATb.clear();
        return;
    }
    int n = ( int )A[ 0 ].size();

    ATA.assign( n, std::vector<double>( n, 0.0 ) );
    ATb.assign( n, 0.0 );

    for ( int i = 0; i < m; ++i )
        {
        for ( int j = 0; j < n; ++j )
        {
            ATb[ j ] += A[ i ][ j ] * b[ i ];
            for ( int k = 0; k < n; ++k )
            {
                ATA[ j ][ k ] += A[ i ][ j ] * A[ i ][ k ];
            }
        }
    }
}

// Solve generic linear system A x = b
bool solve_linear_system( std::vector < double > &x, std::vector < std::vector < double > > &A, std::vector < double > &b )
{
    int n = ( int )A.size();
    if ( n == 0 )
    {
        x.clear();
        return false;
    }

    const double EPS = 1e-15;

    // Forward elimination with partial pivoting
    for ( int col = 0; col < n; ++col )
    {
        // pivot search
        int pivot = col;
        double maxv = std::fabs( A[ col ][ col ] );
        for ( int r = col + 1; r < n; ++r )
        {
            double v = std::fabs( A[ r ][ col ] );
            if ( v > maxv )
            {
                maxv = v;
                pivot = r;
            }
        }

        // swap pivot row
        if ( pivot != col )
        {
            std::swap( A[ pivot ], A[ col ] );
            std::swap( b[ pivot ], b[ col ] );
        }

        // avoid singularities
        if ( std::fabs( A[ col ][ col ] ) < EPS )
        {
            A[ col ][ col ] = ( A[ col ][ col ] >= 0.0 ) ? EPS : -EPS;
        }

        // normalize pivot row
        double piv = A[ col ][ col ];
        for ( int j = col; j < n; ++j )
        {
            A[ col ][ j ] /= piv;
        }
        b[ col ] /= piv;

        // eliminate below
        for ( int r = col + 1; r < n; ++r )
        {
            double f = A[ r ][ col ];
            if ( f == 0.0 )
            {
                continue;
            }
            for ( int j = col; j < n; ++j )
            {
                A[ r ][ j ] -= f * A[ col ][ j ];
            }
            b[ r ] -= f * b[ col ];
        }
    }

    // back substitution
    x.assign( n, 0.0 );
    for ( int i = n - 1; i >= 0; --i )
    {
        double s = b[ i ];
        for ( int j = i + 1; j < n; ++j )
        {
            s -= A[ i ][ j ] * x[ j ];
        }
        x[ i ] = s; // pivot row normalized to 1
    }

    return true;
}

// Solve min ||A x - b|| with A ( m x n ), b ( m )
// returns x ( n ). Uses normal equations ATA x = ATb with pivoting Gaussian elimination.
void solve_least_squares( std::vector < double > &x, const std::vector < std::vector < double > > &A, const std::vector < double > &b )
{
    std::vector < std::vector < double > > ATA;
    std::vector < double > ATb;

    form_normal_equations( ATA, ATb, A, b );
    solve_linear_system( x, ATA, ATb );
}

} // namespace nf
