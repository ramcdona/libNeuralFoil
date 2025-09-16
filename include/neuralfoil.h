//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#ifndef NEURALFOIL_H
#define NEURALFOIL_H

#include <vector>
#include <string>

namespace nf
{

class neuralfoil
{
public:
    neuralfoil();

    bool load( const std::string & datafile, const std::string & scalefile );

    void evaluate( double & analysis_confidence,
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
                   const double & n_crit = 9.0,
                   const double & xtr_upper = 1.0,
                   const double & xtr_lower = 1.0 ) const;

    void evaluate( double & analysis_confidence,
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
                   const double & n_crit = 9.0,
                   const double & xtr_upper = 1.0,
                   const double & xtr_lower = 1.0 ) const;

protected:

    static double sigmoid( double x );
    static double swish( double x, double beta = 1.0 );

    static void multiply( const std::vector < std::vector < double > > & A, std::vector < double > & B, std::vector < double > & C );

    double squared_mahalanobis_distance( const std::vector < double > & x ) const;

    void net( std::vector < double > x, std::vector < double > & y ) const;

    void evaluate( const std::vector < double > & x, std::vector < double > & y ) const;

    static void inputs( std::vector < double > & x,
                        const std::vector < double > & CST_up,
                        const std::vector < double > & CST_low,
                        const double & CST_le,
                        const double & CST_te,
                        const double & alpharad,
                        const double & Re,
                        const double & n_crit,
                        const double & xtr_upper,
                        const double & xtr_lower );

    static void flipx( const std::vector < double > & x, std::vector < double > & xflip );
    static void unflipy( const std::vector < double > & y, std::vector < double > & yflip );
    static void fusey( const std::vector < double > & y, const std::vector < double > & yflip, std::vector < double > & yfuse );

    static void unpacky( const std::vector < double > & y, double Re,
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
                         std::vector < double > & lower_H );

    static void unpacky( const std::vector < double > & y,
                         double & analysis_confidence,
                         double & CL,
                         double & CD,
                         double & CM,
                         double & Top_Xtr,
                         double & Bot_Xtr );

    std::vector < std::vector < std::vector < double > > > m_Weights;
    std::vector < std::vector < double > > m_Biases;

    std::vector < double > m_InputScale;
    std::vector < std::vector < double > > m_InputCovScale;

    bool m_LoadSuccess;
};

} // namespace nf

#endif //NEURALFOIL_H
