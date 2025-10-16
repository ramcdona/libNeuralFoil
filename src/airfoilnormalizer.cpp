//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include "../include/neuralfoil.h"

#include "matsolve.h"

#include <cmath>
#include <vector>

namespace nf
{

AirfoilNormalizer::AirfoilNormalizer()
{
    m_X_Trans = 0.0;
    m_Y_Trans = 0.0;
    m_Scale = 1.0;
    m_Rotate_Deg = 0.0;
}

void AirfoilNormalizer::normalize_airfoil( std::vector < std::vector < double > > & coords_out, const std::vector < std::vector < double > > &coords_in )
{
    int n = ( int ) coords_in.size();

    // TE midpoint ( first and last coordinate points )
    double x_te = 0.5 * ( coords_in[ 0 ][ 0 ] + coords_in[ n - 1 ][ 0 ] );
    double y_te = 0.5 * ( coords_in[ 0 ][ 1 ] + coords_in[ n - 1 ][ 1 ] );

    // distances to TE, find LE index ( max distance )
    double maxd = -1.0;
    int le_index = 0;
    for ( int i = 0; i < n; ++i )
    {
        double dx = coords_in[ i ][ 0 ] - x_te;
        double dy = coords_in[ i ][ 1 ] - y_te;
        double d = std::sqrt( dx * dx + dy * dy );
        if ( d > maxd )
        {
            maxd = d;
            le_index = i;
        }
    }
    // translations to move LE to origin
    m_X_Trans = -coords_in[ le_index ][ 0 ];
    m_Y_Trans = -coords_in[ le_index ][ 1 ];

    // translate
    std::vector < std::vector < double > > translated = std::vector < std::vector < double > >( n, std::vector < double >( 2, 0 ) );
    for ( int i = 0; i < n; ++i )
    {
        translated[ i ][ 0 ] = coords_in[ i ][ 0 ] + m_X_Trans;
        translated[ i ][ 1 ] = coords_in[ i ][ 1 ] + m_Y_Trans;
    }

    // scale so chord = 1 ( chord was maxd )
    m_Scale = 1.0 / maxd;
    std::vector < std::vector < double > > scaled = std::vector < std::vector < double > >( n, std::vector < double >( 2, 0 ) );
    for ( int i = 0; i < n; ++i )
    {
        scaled[ i ][ 0 ] = translated[ i ][ 0 ] * m_Scale;
        scaled[ i ][ 1 ] = translated[ i ][ 1 ] * m_Scale;
    }

    // recompute TE after scaling
    x_te = 0.5 * ( scaled[ 0 ][ 0 ] + scaled[ n - 1 ][ 0 ] );
    y_te = 0.5 * ( scaled[ 0 ][ 1 ] + scaled[ n - 1 ][ 1 ] );

    // rotation to move TE to ( 1,0 )
    double rotation_angle = -std::atan2( y_te, x_te ); // radians
    double cos_a = std::cos( rotation_angle );
    double sin_a = std::sin( rotation_angle );

    coords_out.assign( n, std::vector < double > ( 2, 0 ) );
    for ( int i = 0; i < n; ++i )
    {
        double x = scaled[ i ][ 0 ];
        double y = scaled[ i ][ 1 ];
        coords_out[ i ][ 0 ] = x * cos_a - y * sin_a;
        coords_out[ i ][ 1 ] = x * sin_a + y * cos_a;
    }

    m_Rotate_Deg = rotation_angle * 180.0 / M_PI;
}

void AirfoilNormalizer::denormalize_airfoil( std::vector < std::vector < double > > & coords_out, const std::vector < std::vector < double > > &coords_in ) const
{
    int n = ( int ) coords_in.size();
    coords_out.assign( n, std::vector < double > ( 2, 0 ) );

    double theta = m_Rotate_Deg * M_PI / 180.0;
    double cos_a = std::cos( theta );
    double sin_a = std::sin( theta );

    // inverse rotation: rotate by +theta ( since normalize rotated by -theta )
    for ( int i = 0; i < n; ++i )
    {
        double xn = coords_in[ i ][ 0 ];
        double yn = coords_in[ i ][ 1 ];
        double xr = xn * cos_a + yn * sin_a;
        double yr = -xn * sin_a + yn * cos_a;
        // inverse scale
        xr /= m_Scale;
        yr /= m_Scale;
        // inverse translation
        xr -= m_X_Trans;
        yr -= m_Y_Trans;
        coords_out[ i ][ 0 ] = xr;
        coords_out[ i ][ 1 ] = yr;
    }
}

} // namespace nf
