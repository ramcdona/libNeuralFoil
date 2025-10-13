//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include "test_helpers.h"

#include <cstdio>
#include <cmath>
#include <vector>

bool compare( const double x, const double xref, const bool print )
{
    bool success = true;

    const double err = std::abs( xref - x );

    if ( print )
    {
        printf( "%g %g %g\n", xref, x, err );
    }

    if ( err > 1e-12 )
    {
        success = false;
    }

    printf( "abserr: %g\n", err );

    return success;
}

bool compare( const std::vector < double > & x, const std::vector < double > & xref, bool print )
{
    bool success = true;

    if ( x.size() != xref.size() )
    {
        printf( "size mismatch\n" );
        return false;
    }

    double maxabserr = 0;
    for( size_t i = 0; i < xref.size(); i++ )
    {
        const double err = std::abs( xref[i] - x[i] );
        if( err > maxabserr )
        {
            maxabserr = err;
        }

        if ( print )
        {
            printf( "%ld %g %g %g\n", i, xref[i], x[i], err );
        }
    }

    if ( maxabserr > 1e-12 )
    {
        success = false;
    }

    printf( "maxabserr: %g\n", maxabserr );

    return success;
}


bool compare( const std::vector<std::vector<double>> &X, const std::vector<std::vector<double>> &Xref, bool print )
{
    bool success = true;

    if ( X.size() != Xref.size() )
    {
        printf( "row size mismatch: %zu vs %zu\n", X.size(), Xref.size() );
        return false;
    }
    if ( !X.empty() && X[0].size() != Xref[0].size() )
    {
        printf( "column size mismatch: %zu vs %zu\n", X[0].size(), Xref[0].size() );
        return false;
    }

    const size_t nrows = Xref.size();
    const size_t ncols = nrows > 0 ? Xref[0].size() : 0;

    double maxabserr = 0.0;

    // --- Elementwise comparison ---
    for ( size_t i = 0; i < nrows; ++i )
    {
        for ( size_t j = 0; j < ncols; ++j )
        {
            const double err = std::abs( Xref[i][j] - X[i][j] );
            if ( err > maxabserr )
            {
                maxabserr = err;
            }

            if ( print && err > 1e-12 )
            {
                printf("(%zu,%zu) %g %g %g\n", i, j, Xref[i][j], X[i][j], err);
            }
        }
    }

    if ( maxabserr > 1e-12 )
    {
        success = false;
    }

    printf( "maxabserr: %g\n", maxabserr );
    return success;
}
