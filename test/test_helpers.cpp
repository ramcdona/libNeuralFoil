//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include "test_helpers.h"

#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

bool compare( const double x, const double xref, double tol, const bool print )
{
    bool success = true;

    const double err = std::abs( xref - x );

    if ( print )
    {
        printf( "%g %g %g\n", xref, x, err );
    }

    if ( err > tol )
    {
        success = false;
    }

    printf( "abserr: %g\n", err );

    return success;
}

bool compare( const std::vector < double > & x, const std::vector < double > & xref, double tol, bool print )
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

    if ( maxabserr > tol )
    {
        success = false;
    }

    printf( "maxabserr: %g\n", maxabserr );

    return success;
}


bool compare( const std::vector<std::vector<double>> &X, const std::vector<std::vector<double>> &Xref, double tol, bool print )
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

            if ( print && err > tol )
            {
                printf("(%zu,%zu) %g %g %g\n", i, j, Xref[i][j], X[i][j], err);
            }
        }
    }

    if ( maxabserr > tol )
    {
        success = false;
    }

    printf( "maxabserr: %g\n", maxabserr );
    return success;
}

std::vector < std::vector < double > > generate_NACA4( double m, double p, double t, int n_points, bool cluster )
{
    // returns 2*n_points x 2 coordinates: upper reversed ( TE->LE ) then lower ( LE->TE )
    if ( n_points < 2 ) throw std::runtime_error( "n_points must be >=2" );

    std::vector < double > x( n_points );
    if ( cluster )
    {
        std::vector < double > beta( n_points );
        for ( int i = 0; i < n_points; ++i )
        {
            beta[ i ] = double( i ) * M_PI / double( n_points - 1 );
        }

        for ( int i = 0; i < n_points; ++i )
        {
            x[ i ] = 0.5 * ( 1.0 - std::cos( beta[ i ] ) ); // cosine spacing 0..1
        }
    }
    else
    {
        for ( int i = 0; i < n_points; ++i )
        {
            x[ i ] = double( i ) / double( n_points - 1 );
        }
    }

    std::vector < double > y_c( n_points ), dyc_dx( n_points ), y_t( n_points ), theta( n_points );
    for ( int i = 0; i < n_points; ++i )
    {
        double xi = x[ i ];
        if ( xi < p && p > 0.0 )
        {
            y_c[ i ] = ( m / ( p * p ) ) * ( 2.0 * p * xi - xi * xi );
            dyc_dx[ i ] = ( 2.0 * m / ( p * p ) ) * ( p - xi );
        }
        else if ( p > 0.0 )
        {
            y_c[ i ] = ( m / ( ( 1 - p ) * ( 1 - p ) ) ) * ( ( 1 - 2 * p ) + 2 * p * xi - xi * xi );
            dyc_dx[ i ] = ( 2.0 * m / ( ( 1 - p ) * ( 1 - p ) ) ) * ( p - xi );
        }
        else
        {
            y_c[ i ] = 0.0;
            dyc_dx[ i ] = 0.0;
        }
        y_t[ i ] = 5.0 * t * ( 0.2969 * std::sqrt( xi )
                              - 0.1260 * xi
                              - 0.3516 * xi * xi
                              + 0.2843 * xi * xi * xi
                              - 0.1015 * xi * xi * xi * xi );
        theta[ i ] = std::atan( dyc_dx[ i ] );
    }

    std::vector < std::vector < double > > upper( n_points, std::vector < double >( 2 ) );
    std::vector < std::vector < double > > lower( n_points, std::vector < double >( 2 ) );
    for ( int i = 0; i < n_points; ++i )
    {
        upper[ i ][ 0 ] = x[ i ] - y_t[ i ] * std::sin( theta[ i ] );
        upper[ i ][ 1 ] = y_c[ i ] + y_t[ i ] * std::cos( theta[ i ] );
        lower[ i ][ 0 ] = x[ i ] + y_t[ i ] * std::sin( theta[ i ] );
        lower[ i ][ 1 ] = y_c[ i ] - y_t[ i ] * std::cos( theta[ i ] );
    }

    std::vector < std::vector < double > > coords( 2 * n_points - 1, std::vector < double >( 2 ) );
    // upper reversed: TE -> LE
    for ( int i = 0; i < n_points; ++i )
    {
        coords[ i ][ 0 ] = upper[ n_points - 1 - i ][ 0 ];
        coords[ i ][ 1 ] = upper[ n_points - 1 - i ][ 1 ];
    }
    // lower LE -> TE
    for ( int i = 0; i < n_points - 1; ++i )
    {
        coords[ n_points + i ][ 0 ] = lower[ i + 1 ][ 0 ];
        coords[ n_points + i ][ 1 ] = lower[ i + 1 ][ 1 ];
    }
    return coords;
}

void write_matlab_plot( const std::vector < std::vector < double > > &orig, const std::vector < std::vector < double > > &recon, const std::string &fname )
{
    std::ofstream f( fname );
    f << "% Auto-generated MATLAB script\n";
    f << "orig = [\n";
    for ( auto &p: orig ) f << p[ 0 ] << " " << p[ 1 ] << ";\n";
    f << "];\n";
    f << "recon = [\n";
    for ( auto &p: recon ) f << p[ 0 ] << " " << p[ 1 ] << ";\n";
    f << "];\n";

    f << "figure( 1 ); hold on; grid on;\n";
    f << "plot( orig( :,1 ), orig( :,2 ), 'b.-', 'DisplayName','Original' );\n";
    f << "plot( recon( :,1 ), recon( :,2 ), 'r.-', 'DisplayName','Reconstructed' );\n";
    f << "axis equal; xlabel( 'x' ); ylabel( 'y' ); legend show;\n";
    f << "title( 'Original vs Reconstructed Airfoil' );\n";

    f.close();
}

void write_matlab( const std::vector < int > &m, const std::string & name )
{
    printf( "%s = [", name.c_str() );

    for ( int i = 0; i < m.size(); i++ )
    {
        printf( "%d ", m[ i ] );
    }
    printf( "];\n" );
}

void write_matlab( const std::vector < double > &m, const std::string & name )
{
    printf( "%s = [", name.c_str() );

    for ( int i = 0; i < m.size(); i++ )
    {
        printf( "%.16g ", m[ i ] );
    }
    printf( "];\n" );
}

void write_matlab( const std::vector < std::vector < double > > &m, const std::string & name )
{
    printf( "%s = [", name.c_str() );

    int ni = m.size();
    for ( int i = 0; i < ni; i++ )
    {
        int n = m[ i ].size();
        for ( int j = 0; j < n - 1; j++ )
        {
            printf( "%.16g ", m[ i ][ j ] );
        }

        if ( i < ni - 1 )
        {
            printf( "%.16g;", m[ i ][ n - 1 ] );
        }
        else
        {
            printf( "%.16g", m[ i ][ n - 1 ] );
        }
    }
    printf( "];\n" );
}

