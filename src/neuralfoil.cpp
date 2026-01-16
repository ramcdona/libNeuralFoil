//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include "../include/neuralfoil.h"

#include <map>
#include <set>
#include <filesystem>
#include <cmath>

#include "npy/tensor.h"
#include "npy/npy.h"
#include "npy/npz.h"

namespace nf
{

neuralfoil::neuralfoil()
{
    m_LoadSuccess = false;
}

bool neuralfoil::load( const std::string & datafile, const std::string & scalefile )
{
    // Assume success, set to false for any detected problem.
    m_LoadSuccess = true;

    if ( std::filesystem::exists( scalefile ) )
    {
        try
        {
            npy::inpzstream input( scalefile );

            const npy::tensor < float > & mean_inputs_scaled = input.read < float >( "mean_inputs_scaled.npy" );

            std::vector < size_t > shape = mean_inputs_scaled.shape();

            if ( shape.size() == 1 )
            {
                const int n = shape[ 0 ];
                m_InputScale.resize( n );

                for ( int i = 0; i < n; ++i )
                {
                    m_InputScale[ i ] = mean_inputs_scaled( i );
                }
            }
            else
            {
                m_LoadSuccess = false;
            }


            const npy::tensor < double > & inv_cov_inputs_scaled = input.read < double >( "inv_cov_inputs_scaled.npy" );

            shape = inv_cov_inputs_scaled.shape();
            int rows = static_cast < int >( shape[ 0 ] );
            int cols = static_cast < int >( shape[ 1 ] );

            m_InputCovScale.resize( rows, std::vector < double > ( cols, 0.0 ) );

            for ( int i = 0; i < rows; ++i )
            {
                for ( int k = 0; k < cols; ++k )
                {
                    m_InputCovScale[ i ][ k ] = inv_cov_inputs_scaled( i, k );
                }
            }

        }
        catch ( const std::exception & e )
        {
            printf( "neuralfoil::load() error reading %s.\n%s\n", scalefile.c_str(), e.what() );
            m_LoadSuccess = false;
        }
    }
    else
    {
        printf( "neuralfoil::load() %s not found.\n", scalefile.c_str() );
        m_LoadSuccess = false;
    }

    if ( std::filesystem::exists( datafile ) )
    {
        try
        {
            npy::inpzstream input( datafile );

            std::vector < std::string > keys = input.keys();

            std::map < std::string, npy::tensor < float >  > tensors;

            std::set < int > layer_indices;

            for( size_t i = 0; i < keys.size(); i++ )
            {
                std::string key = keys[ i ];
                size_t firstDot = key.find( '.' );
                size_t secondDot = key.find( '.', firstDot + 1 );
                std::string numberStr = key.substr( firstDot + 1, secondDot - ( firstDot + 1 ) );
                layer_indices.insert( std::stoi( numberStr ) );

                tensors.emplace( keys[i], input.read < float >( keys[i] ) );
            }

            std::vector < int > layer_indx_vec;
            layer_indx_vec.reserve( layer_indices.size() );
            for ( auto it = layer_indices.begin(); it != layer_indices.end(); it++ )
            {
                layer_indx_vec.push_back( *it );
            }


            m_Weights.resize( layer_indx_vec.size() );
            m_Biases.resize( layer_indx_vec.size() );
            for ( size_t ilayer = 0; ilayer < layer_indx_vec.size(); ilayer++ )
            {
                std::string wkey = "net." + std::to_string( layer_indx_vec[ilayer] ) + ".weight.npy";
                const npy::tensor < float > &w = tensors.at( wkey );
                std::string bkey = "net." + std::to_string( layer_indx_vec[ilayer] ) + ".bias.npy";
                const npy::tensor < float > &b = tensors.at( bkey );


                std::vector < size_t > wshape = w.shape();
                int rows = static_cast < int >( wshape[ 0 ] );
                int cols = static_cast < int >( wshape[ 1 ] );

                m_Weights[ ilayer ].resize( rows, std::vector < double > ( cols, 0.0 ) );

                for ( int i = 0; i < rows; ++i )
                {
                    for ( int k = 0; k < cols; ++k )
                    {
                        m_Weights[ ilayer ][ i ][ k ] = w( i, k );
                    }
                }

                std::vector < size_t > bshape = b.shape();
                size_t brows = bshape[ 0 ];
                m_Biases[ ilayer ].resize( brows, 0.0 );

                for ( int i = 0; i < brows; ++i )
                {
                    m_Biases[ ilayer ][ i ] = b( i );
                }
            }
        }
        catch ( std::exception & e )
        {
            printf( "neuralfoil::load() error reading %s.\n%s\n", datafile.c_str(), e.what() );
            m_LoadSuccess = false;
        }
    }
    else
    {
        printf( "neuralfoil::load() %s not found.\n", datafile.c_str() );
        m_LoadSuccess = false;
    }

    return m_LoadSuccess;
}

void neuralfoil::evaluate( double & analysis_confidence,
               double & CL,
               double & CD,
               double & CM,
               double & Top_Xtr,
               double & Bot_Xtr,
               std::vector < double > & upper_bl_ue_over_vinf,
               std::vector < double > & lower_bl_ue_over_vinf,
               std::vector < double > & upper_theta,
               std::vector < double > & lower_theta,
               std::vector < double > & upper_H,
               std::vector < double > & lower_H,
               const std::vector < double > & CST_up,
               const std::vector < double > & CST_low,
               const double & CST_le,
               const double & CST_te,
               const double & alpharad,
               const double & Re,
               const double & n_crit,
               const double & xtr_upper,
               const double & xtr_lower ) const
{
    if ( !m_LoadSuccess )
    {
        printf( "NeuralFoil model not loaded.\n" );
        return;
    }

    std::vector < double > x;

    inputs( x, CST_up, CST_low, CST_le,
            CST_te, alpharad,
            Re, n_crit, xtr_upper,
            xtr_lower );

    std::vector < double > y;
    evaluate( x, y );

    unpacky( y, Re,
           analysis_confidence,
           CL,
           CD,
           CM,
           Top_Xtr,
           Bot_Xtr,
           upper_bl_ue_over_vinf,
           lower_bl_ue_over_vinf,
           upper_theta,
           lower_theta,
           upper_H,
           lower_H );
}

void neuralfoil::evaluate( double & analysis_confidence,
               double & CL,
               double & CD,
               double & CM,
               double & Top_Xtr,
               double & Bot_Xtr,
               const std::vector < double > & CST_up,
               const std::vector < double > & CST_low,
               const double & CST_le,
               const double & CST_te,
               const double & alpharad,
               const double & Re,
               const double & n_crit,
               const double & xtr_upper,
               const double & xtr_lower ) const
{
    if ( !m_LoadSuccess )
    {
        printf( "NeuralFoil model not loaded.\n" );
        return;
    }

    std::vector < double > x;

    inputs( x, CST_up, CST_low, CST_le,
            CST_te, alpharad,
            Re, n_crit, xtr_upper,
            xtr_lower );

    std::vector < double > y;
    evaluate( x, y );

    unpacky( y,
           analysis_confidence,
           CL,
           CD,
           CM,
           Top_Xtr,
           Bot_Xtr );
}

double neuralfoil::sigmoid( const double x )
{
    return 1.0 / ( 1.0 + std::exp( -x ) );
}

double neuralfoil::sigmoid_with_derivative( double x, double & dsigmoid_dx )
{
    double s = sigmoid( x );
    dsigmoid_dx = s * ( 1.0 - s );
    return s;
}

double neuralfoil::swish( const double x, const double beta )
{
    return x * sigmoid( beta * x );
}

double neuralfoil::swish_with_derivative( double x, double & dsw_dx, double beta )
{
    double s = sigmoid( beta * x );
    double sw = s * x;
    // swish'(x) = beta * swish(x) + sigmoid(beta * x) * (1 - beta * swish(x))
    dsw_dx = beta * sw + s * ( 1.0 - beta * sw );
    return sw;
}

void neuralfoil::multiply( const std::vector < std::vector < double > > & A, std::vector < double > & B, std::vector < double > & C )
{
    const size_t aRows = A.size();
    const size_t aCols = A[0].size();

    const size_t bRows = B.size();

    if ( aCols != bRows )
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication");

    C.clear();
    C.resize( aRows, 0.0 );

    for ( int i = 0; i < aRows; ++i )
    {
        for ( int k = 0; k < aCols; ++k )
        {
            C[ i ] += A[ i ][ k ] * B[ k ];
        }
    }
}

void neuralfoil::multiply( const std::vector < std::vector < double > > & A, const std::vector < std::vector < double > > & B, std::vector < std::vector < double > > & C )
{
    const size_t aRows = A.size();
    const size_t aCols = A[ 0 ].size();
    const size_t bCols = B[ 0 ].size();

    C.assign( aRows, std::vector < double > ( bCols, 0.0 ) );

    for ( size_t i = 0; i < aRows; ++i )
    {
        for ( size_t k = 0; k < aCols; ++k )
        {
            for ( size_t j = 0; j < bCols; ++j )
            {
                C[ i ][ j ] += A[ i ][ k ] * B[ k ][ j ];
            }
        }
    }
}

void neuralfoil::add( const std::vector < std::vector < double > > & A,
                      const std::vector < std::vector < double > > & B,
                      std::vector < std::vector < double > > & C )
{
    const size_t rows = A.size();
    const size_t cols = A[ 0 ].size();
    C.assign( rows, std::vector < double > ( cols, 0.0 ) );
    for ( size_t i = 0; i < rows; i++ )
    {
        for ( size_t j = 0; j < cols; j++ )
        {
            C[ i ][ j ] = A[ i ][ j ] + B[ i ][ j ];
        }
    }
}

// v += u * s   where v, u are vectors and s is scalar
void neuralfoil::plus_equals( std::vector < double > & v, const std::vector < double > & u, const double s )
{
    const size_t n = v.size();
    for ( size_t i = 0; i < n; i++ )
    {
        v[ i ] += s * u[ i ];
    }
}

double neuralfoil::squared_mahalanobis_distance( const std::vector < double > & x ) const
{
    const size_t N_inputs = x.size();

    std::vector < double > diff( N_inputs );
    for ( size_t j = 0; j < N_inputs; ++j )
    {
        diff[ j ] = x[ j ] - m_InputScale[ j ];
    }

    std::vector < double > temp;
    multiply( m_InputCovScale, diff, temp );

    double dist = 0.0;
    for ( size_t j = 0; j < N_inputs; ++j )
    {
        dist += diff[ j ] * temp[ j ];
    }

    return dist;
}

void neuralfoil::mahalanobis_distance_derivatives( std::vector < double > & ddist_dx,
                                                   const std::vector < double > & x ) const
{
    const size_t n = x.size();
    std::vector < double > diff( n );
    for ( size_t j = 0; j < n; ++j )
    {
        diff[ j ] = x[ j ] - m_InputScale[ j ];
    }

    // Distance is (x-mu)^T * S^-1 * (x-mu).
    // Since S^-1 is symmetric (covariance inverse), the derivative is 2 * S^-1 * (x-mu).
    std::vector < double > tmp;
    multiply( m_InputCovScale, diff, tmp );

    ddist_dx.resize( n );
    for ( size_t j = 0; j < n; ++j )
    {
        ddist_dx[ j ] = 2.0 * tmp[ j ];
    }
}

void neuralfoil::net( std::vector < double > x, std::vector < double > & y ) const
{

    for ( int ilayer = 0; ilayer < m_Biases.size(); ilayer++ )
    {
        multiply( m_Weights[ ilayer ], x, y );

        for ( int i = 0; i < y.size(); i++ )
        {
            y[ i ] += m_Biases[ ilayer ][ i ];
        }

        if ( ilayer < m_Biases.size() - 1 )
        {
            x.resize( y.size() );
            for ( size_t i = 0; i < y.size(); i++ )
            {
                x[ i ] = swish( y[ i ] );
            }
        }
    }

}

void neuralfoil::net_with_derivatives( std::vector < double > x,
                                       std::vector < double > & y,
                                       std::vector < std::vector < double > > & dy_dx ) const
{

    for ( int ilayer = 0; ilayer < m_Weights.size(); ilayer++ )
    {
        multiply( m_Weights[ ilayer ], x, y );

        for ( int i = 0; i < y.size(); i++ )
        {
            y[ i ] += m_Biases[ ilayer ][ i ];
        }

        std::vector < std::vector < double > > dlayer_dprev = m_Weights[ ilayer ];

        if ( ilayer < m_Weights.size() - 1 )
        {
            x.resize( y.size() );
            for ( size_t i = 0; i < y.size(); i++ )
            {
                double dsw_dx;
                x[ i ] = swish_with_derivative( y[ i ], dsw_dx );
                for ( size_t j = 0; j < dlayer_dprev[ i ].size(); j++ )
                {
                    dlayer_dprev[ i ][ j ] *= dsw_dx;
                }
            }
        }

        if ( ilayer == 0 )
        {
            dy_dx = dlayer_dprev;
        }
        else
        {
            std::vector < std::vector < double > > next_J;
            multiply( dlayer_dprev, dy_dx, next_J );
            dy_dx = next_J;
        }
    }
}

void neuralfoil::evaluate( const std::vector < double > & x, std::vector < double > & y ) const
{
    int N_Input = 25;

    std::vector < double > ynormal;
    std::vector < double > xflip;
    std::vector < double > yflip;
    std::vector < double > yunflip;

    net( x, ynormal );

    // This was baked into training in order to ensure the network asymptotes to zero analysis confidence far away from the training data.
    ynormal[ 0 ] = ynormal[ 0 ] - squared_mahalanobis_distance( x ) / ( 2.0 * N_Input );

    flipx( x, xflip );
    net( xflip, yflip );

    // This was baked into training in order to ensure the network asymptotes to zero analysis confidence far away from the training data.
    yflip[ 0 ] = yflip[ 0 ] - squared_mahalanobis_distance( xflip ) / ( 2.0 * N_Input );

    unflipy( yflip, yunflip );

    fusey( ynormal, yunflip, y );
}

void neuralfoil::inputs( std::vector < double > & x,
                         const std::vector < double > & CST_up,
                         const std::vector < double > & CST_low,
                         const double & CST_le,
                         const double & CST_te,
                         const double & alpharad,
                         const double & Re,
                         const double & n_crit,
                         const double & xtr_upper,
                         const double & xtr_lower )
{
    const double ca = cos ( alpharad );
    const std::vector < double > xother = { CST_le,
                                            CST_te * 50.0,
                                            sin( 2.0 * alpharad ),
                                            ca,
                                            1.0 - ca * ca,
                                            ( log( Re ) - 12.5 ) / 3.5,
                                            ( n_crit - 9.0 ) / 4.5,
                                            xtr_upper,
                                            xtr_lower };
    x = CST_up;
    x.insert( x.end(), CST_low.begin(), CST_low.end() );
    x.insert( x.end(), xother.begin(), xother.end() );
}

void neuralfoil::inputs_derivatives( std::vector < std::vector < double > > & dx_du,
                                     const double & alpharad,
                                     const double & Re )
{
    dx_du.assign( 25, std::vector < double > ( 23, 0.0 ) );
    for ( int i = 0; i < 16; ++i )
    {
        dx_du[ i ][ i ] = 1.0;                                   // CST_up, CST_low
    }
    dx_du[ 16 ][ 16 ] = 1.0;                                     // CST_le
    dx_du[ 17 ][ 17 ] = 50.0;                                    // CST_te
    dx_du[ 18 ][ 18 ] = 2.0 * cos( 2.0 * alpharad );             // sin( 2 * alpha )
    dx_du[ 19 ][ 18 ] = -sin( alpharad );                        // cos( alpha )
    dx_du[ 20 ][ 18 ] = 2.0 * sin( alpharad ) * cos( alpharad ); // 1 - cos^2( alpha )
    dx_du[ 21 ][ 19 ] = 1.0 / ( 3.5 * Re );                      // ( log( Re ) - 12.5 ) / 3.5
    dx_du[ 22 ][ 20 ] = 1.0 / 4.5;                               // ( n_crit - 9.0 ) / 4.5
    dx_du[ 23 ][ 21 ] = 1.0;                                     // xtr_upper
    dx_du[ 24 ][ 22 ] = 1.0;                                     // xtr_lower
}

void neuralfoil::flipx( const std::vector < double > & x, std::vector < double > & xflip )
{
    // Initialize xflip as copy of x
    xflip = x;

    // Flip Kulfan CST coefficients
    const int ncst = 8;
    for ( int i = 0; i < ncst; i++ )
    {
        xflip[ i ] = -1.0 * x[ i + ncst ];
        xflip[ i + ncst ] = - 1.0 * x[ i ];
    }

    xflip[ 16 ] = -1.0 * x [ 16 ];  // CST_le
    xflip[ 18 ] = -1.0 * x [ 18 ];  // sin ( 2.0 * alpharad )

    xflip[ 23 ] = x [ 24 ];  // xtr_upper / xtr_lower
    xflip[ 24 ] = x [ 23 ];  // xtr_lower / xtr_upper
}

void neuralfoil::flipx_derivatives( std::vector < std::vector < double > > & df_dx )
{
    df_dx.assign( 25, std::vector < double > ( 25, 0.0 ) );
    for ( int i = 0; i < 8; i++ )
    {
        df_dx[ i ][ i + 8 ] = -1.0;
        df_dx[ i + 8 ][ i ] = -1.0;
    }
    df_dx[ 16 ][ 16 ] = -1.0;
    df_dx[ 17 ][ 17 ] = 1.0;
    df_dx[ 18 ][ 18 ] = -1.0;
    for ( int i = 19; i <= 22; i++ )
    {
        df_dx[ i ][ i ] = 1.0;
    }
    df_dx[ 23 ][ 24 ] = 1.0;
    df_dx[ 24 ][ 23 ] = 1.0;
}

void neuralfoil::unflipy( const std::vector < double > & y, std::vector < double > & yflip )
{
    yflip = y;

    yflip[ 1 ] = - 1.0 * y[ 1 ];  // CL
    yflip[ 3 ] = - 1.0 * y[ 3 ];  // CM
    yflip[ 4 ] = y[ 5 ];  // switch Top_Xtr with Bot_Xtr
    yflip[ 5 ] = y[ 4 ];  // switch Bot_Xtr with Top_Xtr

    constexpr int isurf = 6; // Index for start of surface data.
    constexpr int nsurf = 32; // Number of surface data points.

    for ( int i = 0; i < nsurf; i++ )
    {
        int j = i + isurf;

        // switch upper and lower Ret, H
        yflip[ j ] = y[ j + nsurf * 3 ];
        yflip[ j + nsurf ] = y[ j + nsurf * 4 ];

        yflip[ j + nsurf * 3 ] = y[ j ];
        yflip[ j + nsurf * 4 ] = y[ j + nsurf ];

        // switch upper_bl_ue/vinf with lower_bl_ue/vinf
        yflip[ j + nsurf * 2 ] = -1.0 * y[ j + nsurf * 5 ];
        yflip[ j + nsurf * 5 ] = -1.0 * y[ j + nsurf * 2 ];
    }
}

void neuralfoil::unflipy_derivatives( std::vector < std::vector < double > > & du_dy )
{
    constexpr int nsurf = 32;
    constexpr int n_out = 6 + nsurf * 6;
    du_dy.assign( n_out, std::vector < double > ( n_out, 0.0 ) );
    du_dy[ 0 ][ 0 ] = 1.0;
    du_dy[ 1 ][ 1 ] = -1.0;
    du_dy[ 2 ][ 2 ] = 1.0;
    du_dy[ 3 ][ 3 ] = -1.0;
    du_dy[ 4 ][ 5 ] = 1.0;
    du_dy[ 5 ][ 4 ] = 1.0;

    constexpr int isurf = 6;
    for ( int i = 0; i < nsurf; i++ )
    {
        int j = i + isurf;
        du_dy[ j ][ j + nsurf * 3 ] = 1.0;
        du_dy[ j + nsurf ][ j + nsurf * 4 ] = 1.0;
        du_dy[ j + nsurf * 2 ][ j + nsurf * 5 ] = -1.0;
        du_dy[ j + nsurf * 3 ][ j ] = 1.0;
        du_dy[ j + nsurf * 4 ][ j + nsurf ] = 1.0;
        du_dy[ j + nsurf * 5 ][ j + nsurf * 2 ] = -1.0;
    }
}

void neuralfoil::fusey( const std::vector < double > & y, const std::vector < double > & yflip, std::vector < double > & yfuse )
{
    yfuse.resize( y.size() );

    for ( int i = 0; i < y.size(); i++ )
    {
        yfuse[ i ] = 0.5 * ( y[ i ] + yflip[ i ] );
    }

    yfuse[ 0 ] = sigmoid( yfuse[ 0 ] );  // Analysis confidence, a binary variable
    yfuse[ 4 ] = std::clamp( yfuse[ 4 ], 0.0, 1.0 );  // Top_Xtr
    yfuse[ 5 ] = std::clamp( yfuse[ 5 ], 0.0, 1.0 );  // Bot_Xtr
}

void neuralfoil::fusey_derivatives( std::vector < std::vector < double > > & dfuse_dy,
                                    std::vector < std::vector < double > > & dfuse_dyunflip,
                                    const std::vector < double > & y,
                                    const std::vector < double > & yunflip )
{
    const size_t n = y.size();
    dfuse_dy.assign( n, std::vector < double > ( n, 0.0 ) );
    dfuse_dyunflip.assign( n, std::vector < double > ( n, 0.0 ) );

    for ( size_t i = 0; i < n; i++ )
    {
        double val = 0.5;
        if ( i == 0 )
        {
            double dsig;
            sigmoid_with_derivative( 0.5 * ( y[ 0 ] + yunflip[ 0 ] ), dsig );
            val = 0.5 * dsig;
        }
        else if ( i == 4 || i == 5 )
        {
            double fused = 0.5 * ( y[ i ] + yunflip[ i ] );
            if ( fused < 0.0 || fused > 1.0 )
            {
                val = 0.0;
            }
        }
        dfuse_dy[ i ][ i ] = val;
        dfuse_dyunflip[ i ][ i ] = val;
    }
}

void neuralfoil::unpacky( const std::vector < double > & y, const double Re,
                          double & analysis_confidence,
                          double & CL,
                          double & CD,
                          double & CM,
                          double & Top_Xtr,
                          double & Bot_Xtr,
                          std::vector < double > & upper_bl_ue_over_vinf,
                          std::vector < double > & lower_bl_ue_over_vinf,
                          std::vector < double > & upper_theta,
                          std::vector < double > & lower_theta,
                          std::vector < double > & upper_H,
                          std::vector < double > & lower_H )
{
    unpacky( y, analysis_confidence, CL, CD, CM, Top_Xtr, Bot_Xtr );

    constexpr int isurf = 6; // Index for start of surface data.
    constexpr int nsurf = 32; // Number of surface data points.

    upper_bl_ue_over_vinf.resize( nsurf );
    lower_bl_ue_over_vinf.resize( nsurf );
    upper_theta.resize( nsurf );
    lower_theta.resize( nsurf );
    upper_H.resize( nsurf );
    lower_H.resize( nsurf );

    for ( int i = 0; i < nsurf; i++ )
    {
        const int j = i + isurf;

        upper_bl_ue_over_vinf[ i ] = y[ j + nsurf * 2 ];
        lower_bl_ue_over_vinf[ i ] = y[ j + nsurf * 5 ];

        upper_theta[ i ] = ( pow( 10.0, y[ j ] ) - 0.1 ) / ( abs( upper_bl_ue_over_vinf[ i ] ) * Re );
        upper_H[ i ] = 2.6 * exp( y[ j + nsurf ] );

        lower_theta[ i ] = ( pow( 10.0, y[ j + nsurf * 3 ] ) - 0.1 ) / ( abs( lower_bl_ue_over_vinf[ i ] ) * Re );
        lower_H[ i ] = 2.6 * exp( y[ j + nsurf * 4 ] );
    }
}

void neuralfoil::unpacky( const std::vector < double > & y,
                          double & analysis_confidence,
                          double & CL,
                          double & CD,
                          double & CM,
                          double & Top_Xtr,
                          double & Bot_Xtr )
{
    analysis_confidence = y[ 0 ];
    CL = 0.5 * y[ 1 ];
    CD = exp( ( y[ 2 ] - 2.0 ) * 2.0 );
    CM = y[ 3 ] / 20.0;
    Top_Xtr = y[ 4 ];
    Bot_Xtr = y[ 5 ];
}

void neuralfoil::unpacky_derivatives( std::vector < std::vector < double > > & dyf_dy,
                                      const std::vector < double > & y_net )
{
    dyf_dy.assign( 6, std::vector < double > ( y_net.size(), 0.0 ) );
    dyf_dy[ 0 ][ 0 ] = 1.0;
    dyf_dy[ 1 ][ 1 ] = 0.5;
    double CD_val = exp( ( y_net[ 2 ] - 2.0 ) * 2.0 );
    dyf_dy[ 2 ][ 2 ] = 2.0 * CD_val;
    dyf_dy[ 3 ][ 3 ] = 1.0 / 20.0;
    dyf_dy[ 4 ][ 4 ] = 1.0;
    dyf_dy[ 5 ][ 5 ] = 1.0;
}

void neuralfoil::unpacky_derivatives( std::vector < std::vector < double > > & dyf_dy,
                                      const std::vector < double > & y_net,
                                      const double Re )
{
    const size_t n_out = y_net.size();
    // Number of scalar outputs (6) + detailed BL outputs (nsurf*6)
    constexpr int nsurf = 32;
    dyf_dy.assign( 6 + nsurf * 6, std::vector < double > ( n_out, 0.0 ) );

    // Scalar outputs
    dyf_dy[ 0 ][ 0 ] = 1.0;
    dyf_dy[ 1 ][ 1 ] = 0.5;
    double CD = exp( ( y_net[ 2 ] - 2.0 ) * 2.0 );
    dyf_dy[ 2 ][ 2 ] = 2.0 * CD;
    dyf_dy[ 3 ][ 3 ] = 1.0 / 20.0;
    dyf_dy[ 4 ][ 4 ] = 1.0;
    dyf_dy[ 5 ][ 5 ] = 1.0;

    // Vector/Surface outputs
    constexpr int isurf = 6;
    for ( int i = 0; i < nsurf; i++ )
    {
        const int j = i + isurf;
        const int row_ue_u = 6 + i;
        const int row_ue_l = 6 + nsurf + i;
        const int row_theta_u = 6 + nsurf * 2 + i;
        const int row_theta_l = 6 + nsurf * 3 + i;
        const int row_H_u = 6 + nsurf * 4 + i;
        const int row_H_l = 6 + nsurf * 5 + i;

        // upper_bl_ue_over_vinf[i] = y[j + nsurf * 2]
        dyf_dy[ row_ue_u ][ j + nsurf * 2 ] = 1.0;
        // lower_bl_ue_over_vinf[i] = y[j + nsurf * 5]
        dyf_dy[ row_ue_l ][ j + nsurf * 5 ] = 1.0;

        double ue_u = y_net[ j + nsurf * 2 ];
        double ue_l = y_net[ j + nsurf * 5 ];

        // upper_theta[i] = (pow(10, y[j]) - 0.1) / (abs(ue_u) * Re)
        double term_u = pow( 10.0, y_net[ j ] );
        dyf_dy[ row_theta_u ][ j ] = ( log( 10.0 ) * term_u ) / ( abs( ue_u ) * Re );
        dyf_dy[ row_theta_u ][ j + nsurf * 2 ] = - ( term_u - 0.1 ) / ( ue_u * abs( ue_u ) * Re );

        // lower_theta[i] = (pow(10, y[j + nsurf * 3]) - 0.1) / (abs(ue_l) * Re)
        double term_l = pow( 10.0, y_net[ j + nsurf * 3 ] );
        dyf_dy[ row_theta_l ][ j + nsurf * 3 ] = ( log( 10.0 ) * term_l ) / ( abs( ue_l ) * Re );
        dyf_dy[ row_theta_l ][ j + nsurf * 5 ] = - ( term_l - 0.1 ) / ( ue_l * abs( ue_l ) * Re );

        // upper_H[i] = 2.6 * exp(y[j + nsurf])
        dyf_dy[ row_H_u ][ j + nsurf ] = 2.6 * exp( y_net[ j + nsurf ] );
        // lower_H[i] = 2.6 * exp(y[j + nsurf * 4])
        dyf_dy[ row_H_l ][ j + nsurf * 4 ] = 2.6 * exp( y_net[ j + nsurf * 4 ] );
    }
}

} // namespace nf
