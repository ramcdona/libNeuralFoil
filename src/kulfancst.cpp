//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#define _USE_MATH_DEFINES
#include <cmath>

#include "../include/neuralfoil.h"

#include "matsolve.h"

#include <vector>

#include <algorithm>  // for std::swap

namespace nf
{

KulfanCST::KulfanCST()
{
    m_N1 = 0.5;
    m_N2 = 1.0;
    m_LE_Weight = 0.0;
    m_TE_Thickness = 0.0;
}

void KulfanCST::fit_kulfan( const std::vector < std::vector < double > > &coords, int n_weights_per_side )
{
    int n_coordinates = ( int ) coords.size();

    // extract x and y arrays
    std::vector < double > x( n_coordinates ), y( n_coordinates );
    for ( int i = 0; i < n_coordinates; ++i )
    {
        x[ i ] = coords[ i ][ 0 ];
        y[ i ] = coords[ i ][ 1 ];
    }

    int LE_index = 0;
    double xmin = x[ 0 ];
    for ( int i = 1; i < n_coordinates; ++i )
    {
        if ( x[ i ] < xmin )
        {
            xmin = x[ i ];
            LE_index = i;
        }
    }


    // Class function C = x^N1 * ( 1-x )^N2
    std::vector < double > Cvec( n_coordinates );
    for ( int i = 0; i < n_coordinates; ++i )
    {
        Cvec[ i ] = std::pow( x[ i ], m_N1 ) * std::pow( 1.0 - x[ i ], m_N2 );
    }

    // Bernstein polynomial order
    int N = n_weights_per_side - 1;

    // Precompute S_matrix ( n_weights_per_side x n_coordinates ): S( k,i ) = K_k * x_i^k * ( 1-x_i )^( N-k )
    std::vector < std::vector < double > > S_matrix( n_weights_per_side, std::vector < double >( n_coordinates, 0.0 ) );
    for ( int k = 0; k <= N; ++k )
    {
        double K = binomial_coefficient( N, k );
        for ( int i = 0; i < n_coordinates; ++i )
        {
            double val = K * std::pow( x[ i ], k ) * std::pow( 1.0 - x[ i ], N - k );
            S_matrix[ k ][ i ] = val;
        }
    }

    // leading edge weight row
    std::vector < double > leading_edge_weight_row( n_coordinates, 0.0 );
    double le_exp = double( n_weights_per_side ) + 0.5;
    for ( int i = 0; i < n_coordinates; ++i )
    {
        leading_edge_weight_row[ i ] = x[ i ] * std::pow( std::max( 0.0, 1.0 - x[ i ] ), le_exp );
    }

    // trailing edge thickness row: x/2 for upper, -x/2 for lower
    std::vector < double > trailing_edge_row( n_coordinates, 0.0 );
    for ( int i = 0; i < n_coordinates; ++i )
    {
        double sign;
        if ( i <= LE_index )
        {
            sign = 1.0;
        }
        else
        {
            sign = -1.0;
        }
        trailing_edge_row[ i ] = sign * x[ i ] / 2.0;
    }

    // Build A matrix: columns = [ lower_weights ( n_w ), upper_weights ( n_w ), leading_edge, trailing_edge ]
    int n_unknowns = 2 * n_weights_per_side + 2;
    std::vector < std::vector < double > > A( n_coordinates, std::vector < double >( n_unknowns, 0.0 ) );
    for ( int row = 0; row < n_coordinates; ++row )
    {
        if ( row > LE_index ) // lower weights columns ( 0..n_w-1 ) only if not is_upper
        {
            for ( int k = 0; k < n_weights_per_side; ++k )
            {
                A[ row ][ k ] = Cvec[ row ] * S_matrix[ k ][ row ];
            }
        }
        else // upper weights columns ( n_w .. 2*n_w -1 ) only if is_upper
        {
            for ( int k = 0; k < n_weights_per_side; ++k )
            {
                A[ row ][ n_weights_per_side + k ] = Cvec[ row ] * S_matrix[ k ][ row ];
            }
        }
        // leading edge column
        A[ row ][ 2 * n_weights_per_side ] = leading_edge_weight_row[ row ];
        // trailing edge
        A[ row ][ 2 * n_weights_per_side + 1 ] = trailing_edge_row[ row ];
    }

    // b vector is y
    std::vector < double > b = y;

    // Solve least squares
    std::vector < double > xsol;
    solve_least_squares( xsol, A, b );
    // Extract results
    m_Lower_Weights.resize( n_weights_per_side );
    m_Upper_Weights.resize( n_weights_per_side );
    for ( int i = 0; i < n_weights_per_side; ++i )
    {
        m_Lower_Weights[ i ] = xsol[ i ];
    }
    for ( int i = 0; i < n_weights_per_side; ++i )
    {
        m_Upper_Weights[ i ] = xsol[ n_weights_per_side + i ];
    }
    m_LE_Weight = xsol[ 2 * n_weights_per_side ];
    m_TE_Thickness = xsol[ 2 * n_weights_per_side + 1 ];

    // If negative TE_thickness, refit without last column
    if ( m_TE_Thickness < 0.0 )
    {
        // build A2 without last column
        int n_unknowns2 = n_unknowns - 1;
        std::vector < std::vector < double > > A2( n_coordinates, std::vector < double >( n_unknowns2, 0.0 ) );
        for ( int row = 0; row < n_coordinates; ++row )
        {
            for ( int col = 0; col < n_unknowns2; ++col )
            {
                A2[ row ][ col ] = A[ row ][ col ];
            }
        }
        std::vector < double > xsol2;
        solve_least_squares( xsol2, A2, b );
        for ( int i = 0; i < n_weights_per_side; ++i )
        {
            m_Lower_Weights[ i ] = xsol2[ i ];
        }
        for ( int i = 0; i < n_weights_per_side; ++i )
        {
            m_Upper_Weights[ i ] = xsol2[ n_weights_per_side + i ];
        }
        m_LE_Weight = xsol2[ 2 * n_weights_per_side ];
        m_TE_Thickness = 0.0;
    }
}

void KulfanCST::evaluate_kulfan( std::vector < std::vector < double > > & coords_out, int n_points_per_side ) const
{
    // produce x linspace 0..1 length = n_points_per_side
    std::vector < double > x( n_points_per_side );
    for ( int i = 0; i < n_points_per_side; ++i )
    {
        x[ i ] = double( i ) / double( n_points_per_side - 1 );
    }

    int n_w = ( int ) m_Upper_Weights.size();

    // compute y_upper, y_lower in normalized space
    std::vector < double > yU( n_points_per_side, 0.0 ), yL( n_points_per_side, 0.0 );
    for ( int i = 0; i < n_points_per_side; ++i )
    {
        double xi = x[ i ];
        double C = std::pow( xi, m_N1 ) * std::pow( 1.0 - xi, m_N2 );
        // shape function sum
        double SU = 0.0, SL = 0.0;
        for ( int k = 0; k < n_w; ++k )
        {
            double K = binomial_coefficient( n_w - 1, k );
            double basis = K * std::pow( xi, k ) * std::pow( 1.0 - xi, n_w - 1 - k );
            SU += m_Upper_Weights[ k ] * basis;
            SL += m_Lower_Weights[ k ] * basis;
        }
        yU[ i ] = C * SU;
        yL[ i ] = C * SL;
    }

    // Add TE thickness
    for ( int i = 0; i < n_points_per_side; ++i )
    {
        yU[ i ] += x[ i ] * m_TE_Thickness / 2.0;
        yL[ i ] -= x[ i ] * m_TE_Thickness / 2.0;
    }

    // Add LEM
    double exp_LE = double( n_w ) + 0.5;
    for ( int i = 0; i < n_points_per_side; ++i )
    {
        double xi = x[ i ];
        double le_term = xi * std::pow( 1.0 - xi, exp_LE );
        yU[ i ] += m_LE_Weight * le_term;
        yL[ i ] += m_LE_Weight * le_term;
    }

    // Build coords in same ordering as generator: upper reversed ( TE->LE ) then lower ( LE->TE )
    coords_out.assign( 2 * n_points_per_side, std::vector < double > ( 2, 0.0 ) );
    for ( int i = 0; i < n_points_per_side; ++i )
    {
        coords_out[ i ][ 0 ] = x[ n_points_per_side - 1 - i ];
        coords_out[ i ][ 1 ] = yU[ n_points_per_side - 1 - i ];
    }
    for ( int i = 0; i < n_points_per_side; ++i )
    {
        coords_out[ n_points_per_side + i ][ 0 ] = x[ i ];
        coords_out[ n_points_per_side + i ][ 1 ] = yL[ i ];
    }
}

double KulfanCST::binomial_coefficient( int n, int k )
{
    if ( k < 0 || k > n )
    {
        return 0.0;
    }
    k = std::min( k, n - k );
    double res = 1.0;
    for ( int i = 1; i <= k; ++i )
    {
        res *= double( n - ( k - i ) );
        res /= double( i );
    }
    return res;
}

} // namespace nf
