//
// This file is released under the terms of the MIT License as detailed in the
// LICENSE.txt file which accompanies this software.
//

#include "neuralfoil.h"
#include <vector>

int main()
{
    // CST coefficients for NACA 4412 airfoil.
    std::vector < double > CST_up = {
        0.18109497123192536, 0.2126841855378795, 0.28098503347171055,
        0.24864887442714956, 0.2402814006497832, 0.2726284343621469,
        0.2577647399772107, 0.27817638080391327
    };

    std::vector < double > CST_low = {
        -0.16965145937914555, -0.09364138259246524, -0.06345895875250755,
        -0.006796600807436377, -0.09024470493931744, 0.020818449122134496,
        -0.03575216051938629, -0.0022362344439678192
    };

    double CST_le = 0.10647339061374347;
    double CST_te = 0.0025720113171501307;

    double alpha = 5.159747295611583;
    double Re = 1000302.1671445925;

    double n_crit = 9.0;
    double xtr_upper = 1.0;
    double xtr_lower = 1.0;

    // Outputs
    double analysis_confidence;
    double CL;
    double CD;
    double CM;
    double Top_Xtr;
    double Bot_Xtr;
    std::vector < double > upper_bl_ue_over_vinf;
    std::vector < double > lower_bl_ue_over_vinf;
    std::vector < double > upper_theta;
    std::vector < double > lower_theta;
    std::vector < double > upper_H;
    std::vector < double > lower_H;

    nf::neuralfoil nf;
    nf.load( "nn-xlarge.npz", "scaled_input_distribution.npz" );

    nf.evaluate( analysis_confidence,
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
                 lower_H,
                 CST_up,
                 CST_low,
                 CST_le,
                 CST_te,
                 alpha * M_PI / 180.0,
                 Re,
                 n_crit,
                 xtr_upper,
                 xtr_lower );

    return 0;
}
