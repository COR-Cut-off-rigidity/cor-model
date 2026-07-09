#include <math.h>
#include <string.h>
#include "crustal.h"

#include <chaos_crustal_7_16.h>
#include <chaos_crustal_7_17.h>
#include <chaos_crustal_7_18.h>
#include <chaos_crustal_8_1.h>
#include <chaos_crustal_8_2.h>
#include <chaos_crustal_8_3.h>

void chaos_geo(double XGEO, double YGEO, double ZGEO, const CrustalModelCoefs *coefs, double *BXGEO, double *BYGEO, double *BZGEO) {
    double RHO2 = XGEO * XGEO + YGEO * YGEO;
    double R = sqrt(RHO2 + ZGEO * ZGEO);
    const double C = ZGEO / R;
    double r = R * DEF_RE_KM_D, theta = 0., phi = 0., S = 0., CF = 1., SF = 0.;
    
    if (RHO2 != 0.) {
        const double RHO = sqrt(RHO2);
        S = RHO / R;
        CF = XGEO/RHO;
        SF = YGEO/RHO;
        phi = atan2(YGEO, XGEO);
        theta = atan2(RHO, ZGEO);
        
        if (phi < 0.) {
            phi += DEF_2PI_D;
        }
    } else if (ZGEO < 0.) {
        theta = DEF_PI_D;
    }
    
    double B_r = 0.;
    double B_theta = 0.;
    double B_phi = 0.;
    int nm_max = coefs->nm_max;

    double tem = DEF_RE_KM_D / r;
    double tem_r[nm_max];
    tem_r[0] = tem;

    double dtem_Rr[nm_max];
    dtem_Rr[0] = -tem/r;

    for (int i = 1; i < nm_max; i++) {
        double tem_r_new = tem_r[i - 1] * tem;
        tem_r[i] = tem_r_new;
        dtem_Rr[i] = -(i + 1) * tem_r_new / r;
    }

    double pmn2[nm_max][nm_max];
    double dpmn2[nm_max][nm_max];
    memset(pmn2, 0, nm_max*nm_max*sizeof(double));
    memset(dpmn2, 0, nm_max*nm_max*sizeof(double));

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    pmn2[0][0] = 1.;
    pmn2[1][0] = cos_theta;
    pmn2[1][1] = sin_theta;

    dpmn2[0][0] = 0.;
    dpmn2[1][0] = -sin_theta;
    dpmn2[1][1] = cos_theta;

    for (int n = 2; n < nm_max; n++) {
        double n2 = n * 2;
        double n2_1 = n2 - 1;
        int n_1 = n - 1;

        double temp_mn_2 = sqrt(n2_1 / n2);
        double temp_mn_2_pmn2_dec = temp_mn_2 * pmn2[n_1][n_1];
        pmn2[n][n] = sin_theta * temp_mn_2_pmn2_dec;
        dpmn2[n][n] = cos_theta * temp_mn_2_pmn2_dec + sin_theta * temp_mn_2 * dpmn2[n_1][n_1];

        int n_2 = n - 2;
        double nn = (double) n * n;
        for (int m = 0; m < n; m++) {
            double mm = m * m;
            double nn_mm = nn - mm;
            double temp_mn_2 = n2_1 / sqrt(nn_mm);
            double temp_mn_22 = sqrt((n_1 * n_1 - mm) / nn_mm);
            double temp_pmn2 = pmn2[n_1][m];
            pmn2[n][m] = temp_mn_2 * cos_theta * temp_pmn2 - temp_mn_22 * pmn2[n_2][m];
            dpmn2[n][m] = temp_mn_2 * (-sin_theta * temp_pmn2 + cos_theta * dpmn2[n_1][m]) - temp_mn_22 * dpmn2[n_2][m];
        }
    }

    const double (*G)[nm_max] = coefs->G;
    const double (*H)[nm_max] = coefs->H;

    for (int n = 0; n < nm_max; n++) {
        for (int m = 0; m <= n; m++) {
            double m_phi = m*phi;
            double sin_m_phi = sin(m_phi);
            double cos_m_phi = cos(m_phi);
            double temp_G = G[n][m];
            double temp_H = H[n][m];
            double tem = temp_G * cos_m_phi + temp_H * sin_m_phi;
            double temp_pmn2 = pmn2[n][m];
            double temp_tem_r = tem_r[n];

            B_r += tem * temp_pmn2 * dtem_Rr[n];
            B_theta += tem * dpmn2[n][m] * temp_tem_r;
            B_phi += (-temp_G * m * sin_m_phi + temp_H * m * cos_m_phi) * temp_pmn2 * temp_tem_r;
        }
    }

    B_theta = -B_theta * DEF_RE_KM_D / r;
    B_phi = -B_phi * DEF_RE_KM_D / (r * sin_theta);
    B_r = -B_r * DEF_RE_KM_D;
    const double HE = B_r*S + B_theta*C;

    *BXGEO = HE*CF - B_phi*SF;
    *BYGEO = HE*SF + B_phi*CF;
    *BZGEO = B_r * C - B_theta*S;
}
