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

/**
 * @class AirfoilNormalizer
 * @brief Class for normalizing and denormalizing airfoil coordinates.
 *
 * This class handles the transformation of airfoil coordinates to a normalized
 * form and back. Normalization includes translation, scaling, and rotation operations.
 */
class AirfoilNormalizer
{
public:
    /**
     * @brief Default constructor for AirfoilNormalizer.
     */
    AirfoilNormalizer();

    /**
     * @brief Normalize airfoil coordinates.
     *
     * Transforms the input airfoil coordinates to a normalized form by applying
     * translation, scaling, and rotation. The transformation parameters are stored
     * internally for later denormalization.
     *
     * @param[out] coords_out Normalized airfoil coordinates (2D vector: [[x1,y1], [x2,y2], ...])
     * @param[in] coords_in Input airfoil coordinates to be normalized
     * @param[in] fit_le Flag to use a quadratic fit to find the LE point
     */
    void normalize_airfoil( std::vector < std::vector < double > > & coords_out, std::vector < std::vector < double > > coords_in, bool fit_le = true );

    /**
     * @brief Denormalize airfoil coordinates.
     *
     * Transforms normalized airfoil coordinates back to their original form using
     * the transformation parameters stored during normalization.
     *
     * @param[out] coords_out Denormalized airfoil coordinates (2D vector: [[x1,y1], [x2,y2], ...])
     * @param[in] coords_in Normalized airfoil coordinates to be denormalized
     */
    void denormalize_airfoil( std::vector < std::vector < double > > & coords_out, const std::vector < std::vector < double > > &coords_in ) const;

    double m_X_Trans;       ///< X-axis translation applied during normalization
    double m_Y_Trans;       ///< Y-axis translation applied during normalization
    double m_Scale;         ///< Scale factor applied during normalization (1/chord)
    double m_Rotate_Deg;    ///< Rotation angle in degrees (counter-clockwise positive)
};

/**
 * @class KulfanCST
 * @brief Class for Kulfan CST (Class Shape Transformation) airfoil parameterization.
 *
 * This class implements the Kulfan CST method for representing airfoil shapes using
 * Bernstein polynomials. It can fit CST parameters to airfoil coordinates and
 * evaluate coordinates from CST parameters.
 */
class KulfanCST
{
public:
    /**
     * @brief Default constructor for KulfanCST.
     */
    KulfanCST();

    /**
     * @brief Fit Kulfan CST parameters to airfoil coordinates.
     *
     * Computes the CST weights that best represent the given airfoil coordinates
     * using a least-squares fitting approach.
     *
     * @param[in] coords Airfoil coordinates to fit (2D vector: [[x1,y1], [x2,y2], ...])
     * @param[in] n_weights_per_side Number of CST weights for each surface (upper/lower), default is 8
     */
    void fit_kulfan( const std::vector < std::vector < double > > &coords, int n_weights_per_side = 8 );

    /**
     * @brief Evaluate airfoil coordinates from Kulfan CST parameters.
     *
     * Generates airfoil coordinates by evaluating the CST representation using
     * the stored weight parameters.
     *
     * @param[out] coords_out Generated airfoil coordinates (2D vector: [[x1,y1], [x2,y2], ...])
     * @param[in] n_points_per_side Number of points to generate for each surface (upper/lower)
     */
    void evaluate_kulfan( std::vector < std::vector < double > > & coords_out, int n_points_per_side ) const;

    std::vector < double > m_Lower_Weights;  ///< CST weights for the lower surface
    std::vector < double > m_Upper_Weights;  ///< CST weights for the upper surface
    double m_LE_Weight;                      ///< Leading edge weight parameter
    double m_TE_Thickness;                   ///< Trailing edge thickness
    double m_N1;                             ///< CST class function exponent N1
    double m_N2;                             ///< CST class function exponent N2

protected:
    static double binomial_coefficient( int n, int k );
};

/**
 * @class neuralfoil
 * @brief Neural network-based airfoil aerodynamic analysis.
 *
 * This class provides a machine learning-based approach to predict airfoil
 * aerodynamic characteristics including lift coefficient (CL), drag coefficient (CD),
 * moment coefficient (CM), transition locations, and boundary layer properties.
 * The neural network is trained on airfoil analysis data and uses CST parameterization
 * as input along with flow conditions.
 */
class neuralfoil
{
public:
    /**
     * @brief Default constructor for neuralfoil.
     */
    neuralfoil();

    /**
     * @brief Load neural network model and scaling parameters.
     *
     * Loads the trained neural network weights and biases from the data file,
     * and the input scaling parameters from the scale file. Both files must be
     * successfully loaded for the model to work.
     *
     * @param[in] datafile Path to the file containing neural network weights and biases
     * @param[in] scalefile Path to the file containing input scaling parameters
     * @return true if both files were loaded successfully, false otherwise
     */
    bool load( const std::string & datafile, const std::string & scalefile );

    /**
     * @brief Evaluate airfoil aerodynamics with detailed boundary layer output.
     *
     * Predicts aerodynamic characteristics and boundary layer properties for an airfoil
     * defined by CST parameters at specified flow conditions. This version returns
     * detailed boundary layer information including edge velocity, momentum thickness,
     * and shape factor distributions.
     *
     * @param[out] analysis_confidence Confidence level of the analysis (0-1, higher is better)
     * @param[out] CL Lift coefficient
     * @param[out] CD Drag coefficient
     * @param[out] CM Moment coefficient
     * @param[out] Top_Xtr Transition location on upper surface (x/c)
     * @param[out] Bot_Xtr Transition location on lower surface (x/c)
     * @param[out] upper_bl_ue_over_vinf Upper surface boundary layer edge velocity ratio (ue/Vinf)
     * @param[out] lower_bl_ue_over_vinf Lower surface boundary layer edge velocity ratio (ue/Vinf)
     * @param[out] upper_theta Upper surface boundary layer momentum thickness
     * @param[out] lower_theta Lower surface boundary layer momentum thickness
     * @param[out] upper_H Upper surface boundary layer shape factor
     * @param[out] lower_H Lower surface boundary layer shape factor
     * @param[in] CST_up CST weights for upper surface
     * @param[in] CST_low CST weights for lower surface
     * @param[in] CST_le Leading edge CST weight
     * @param[in] CST_te Trailing edge thickness parameter
     * @param[in] alpharad Angle of attack in radians
     * @param[in] Re Reynolds number
     * @param[in] n_crit Transition criterion (Ncrit), default is 9.0
     * @param[in] xtr_upper Forced transition location on upper surface (x/c), default is 1.0 (free transition)
     * @param[in] xtr_lower Forced transition location on lower surface (x/c), default is 1.0 (free transition)
     */
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

    /**
     * @brief Evaluate airfoil aerodynamics without boundary layer details.
     *
     * Predicts basic aerodynamic characteristics for an airfoil defined by CST parameters
     * at specified flow conditions. This is a simplified version that returns only the
     * integral quantities (CL, CD, CM) and transition locations.
     *
     * @param[out] analysis_confidence Confidence level of the analysis (0-1, higher is better)
     * @param[out] CL Lift coefficient
     * @param[out] CD Drag coefficient
     * @param[out] CM Moment coefficient
     * @param[out] Top_Xtr Transition location on upper surface (x/c)
     * @param[out] Bot_Xtr Transition location on lower surface (x/c)
     * @param[in] CST_up CST weights for upper surface
     * @param[in] CST_low CST weights for lower surface
     * @param[in] CST_le Leading edge CST weight
     * @param[in] CST_te Trailing edge thickness parameter
     * @param[in] alpharad Angle of attack in radians
     * @param[in] Re Reynolds number
     * @param[in] n_crit Transition criterion (Ncrit), default is 9.0
     * @param[in] xtr_upper Forced transition location on upper surface (x/c), default is 1.0 (free transition)
     * @param[in] xtr_lower Forced transition location on lower surface (x/c), default is 1.0 (free transition)
     */
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

    static void multiply( const std::vector < std::vector < double > > & A, const std::vector < std::vector < double > > & B, std::vector < std::vector < double > > & C );

    static void add( const std::vector < std::vector < double > > & A, const std::vector < std::vector < double > > & B, std::vector < std::vector < double > > & C );

    static void plus_equals( std::vector < double > & v, const std::vector < double > & u, double s );
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
