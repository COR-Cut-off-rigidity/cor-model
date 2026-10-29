#include <math.h>
#include <stdlib.h>
#include <macro_defs.h>
#include "internal.h"

// IGRF coefficients
#include <igrf_9.h>
#include <igrf_10.h>
#include <igrf_11.h>
#include <igrf_12.h>
#include <igrf_13.h>
#include <igrf_14.h>

// Historical coefficients
#include <historical.h>
#include <cals10k_2.h>

// CHAOS coefficients
#include <chaos_core_7_16.h>
#include <chaos_core_7_17.h>
#include <chaos_core_7_18.h>
#include <chaos_core_8_1.h>
#include <chaos_core_8_2.h>
#include <chaos_core_8_3.h>

int bin_search_index(double year, const GHTime **coefs, int coefs_size) {
    int low_index = 0;
    int high_index = coefs_size - 1;

    if (year < coefs[0]->year || year >= coefs[high_index]->year) {
        return -1;
    }

    while (low_index <= high_index) {
        int mid_index = low_index + (high_index - low_index) / 2;
        double mid_year = coefs[mid_index]->year;
        double mid_year2 = coefs[mid_index + 1]->year;

        if (mid_year <= year && mid_year2 > year) {
            return mid_index;
        }

        if (mid_year < year) {
            low_index = mid_index + 1;
            continue;
        }

        high_index = mid_index - 1;
    }

    return -1;
}

int internal_recalc(const InternalModelCoefs *model_coefs, const CalcTime *time, InternalModelContext *internal_model_context) {
    double year = (double) time->year;
    double year_decimal = year + (time->year_day - 1) / 365.25;

    int coef_index = bin_search_index(year_decimal, model_coefs->coefs, model_coefs->coefs_size);
    if (coef_index == -1) {
        return -1;
    }

    int nm_max = model_coefs->nm_max;
    int nm_max_1 = nm_max + 1;
    int coefs_size = (nm_max_1 * nm_max) / 2 + nm_max_1;

    double *coefs = malloc(coefs_size*3 * sizeof(double));
    if (coefs == NULL) {
        return -1;
    }

    const GHTime *older_coefs = model_coefs->coefs[coef_index];
    const GHTime *newer_coefs = model_coefs->coefs[coef_index + 1];
    const double *Gx = older_coefs->G;
    const double *Hx = older_coefs->H;
    const double *Gy = newer_coefs->G;
    const double *Hy = newer_coefs->H;

    if (model_coefs->extrapolate && model_coefs->coefs_size - coef_index <= 2) {
        // extrapolate
        double DT = year_decimal - older_coefs->year;
        for (int i = 0, j = 0; i < coefs_size; ++i, j += 3) {
            coefs[j] = Gx[i] + Gy[i] * DT;
            coefs[j + 1] = Hx[i] + Hy[i] * DT;
        }
    } else {
        // interpolate
        double years_diff = newer_coefs->year - older_coefs->year;
        double F2 = (year_decimal - older_coefs->year) / years_diff;
        double F1 = 1. - F2;
        for (int i = 0, j = 0; i < coefs_size; ++i, j += 3) {
            coefs[j] = Gx[i] * F1 + Gy[i] * F2;
            coefs[j + 1] = Hx[i] * F1 + Hy[i] * F2;
        }
    }

    int index = 3;
    double S = 1., P;
    for (int i = 1; i < 14; ++i) {
        S = S * (2*i - 1) / i;
        coefs[index] *= S;
        coefs[index + 1] *= S;
        index += 3;
        P = S;
        
        for (int j = 1; j <= i; ++j) {
            double AA = 1.;
            if (j == 1) {
                AA = 2.;
            }
            P = P * sqrt(AA * (i - j + 1) / (i + j));
            coefs[index] *= P;
            coefs[index + 1] *= P;
            index += 3;
        }
    }

    index = 2;
    for (int i = 0; i < 14; ++i) {
        int ii = i*i;
        double N2 = (double) (4*ii - 1);

        for (int j = 0; j <= i; ++j) {
            int jj = j*j;
            coefs[index] = (ii - jj) / N2;
            index += 3;
        }
    }

    internal_model_context->coefs = coefs;
    internal_model_context->nm_max = nm_max;
    return 0;
}

void internal_clear(InternalModelContext *internal_model_context) {
    free(internal_model_context->coefs);
    internal_model_context->coefs = NULL;
}

void internal_geo(const InternalModelContext *internal_model_context, double XGEO, double YGEO, double ZGEO, double *BXGEO, double *BYGEO, double *BZGEO) {
    const double RHO2 = XGEO*XGEO + YGEO*YGEO;
    const double R = sqrt(RHO2 + ZGEO*ZGEO);
    const double C = ZGEO/R;
    const double RHO = sqrt(RHO2);
    const double S = RHO/R;

    double CF, SF;
    if(S < 1.E-5) {
        CF = 1.;
        SF = 0.;
    } else {
        CF = XGEO/RHO;
        SF = YGEO/RHO;
    }

    const int IRP3 = (int) R + 2;
    int NM = 3 + 30/IRP3;
    if (NM > internal_model_context->nm_max) {
        NM = internal_model_context->nm_max;
    }
    
    const int K = NM+1;

    double A[2*K];
    const double PP = 1./R;
    double P = PP;
    for (int i = 1, j = 0; i <= K; ++i, j += 2) {
        P = P * PP;
        A[j] = P;
        A[j + 1] = P * i;
    }

    P = 1.;

    double D = 0., BBR = 0., BBT = 0., BBF = 0., X = 0., Y = 1., Q = 1., Z = 0.;
    for (int i = 0; i < K; i++) {
        double BI = 0.;
        double P2 = 0.;
        double D2 = 0.;

        for (int j = i; j < K; j++) {
            int MN = (j*(j + 1)/2 + i)*3;
            double E = internal_model_context->coefs[MN];
            double HH = internal_model_context->coefs[MN+1];
            double XK = internal_model_context->coefs[MN+2];
            double W = E*Y + HH*X;

            int a_j = 2*j;
            double AN = A[a_j];
            BBR = BBR + A[a_j + 1]*W*Q;
            BBT = BBT - AN*W*Z;
            BI = BI + AN*(E*X - HH*Y)*Q;

            double DP = C*Z - S*Q - XK*D2;
            double PM = C*Q - XK*P2;
            D2 = Z;
            P2 = Q;
            Z = DP;
            Q = PM;
        }

        D = S*D + C*P;
        P = S*P;
        Q = P;
        Z = D;
        
        BBF = BBF + BI*i;

        double W = X;
        X = W*CF + Y*SF;
        Y = Y*CF - W*SF;
    }

    const double BR = BBR;
    const double BT = BBT;
    double BF;
    if (S < 1.E-5) {
        if (C < 0.) {
            BBF = -BBF;
        }
        BF = BBF;
    } else {
        BF = BBF/S;
    }

    const double HE = BR*S + BT*C;

    *BXGEO = HE*CF - BF*SF;
    *BYGEO = HE*SF + BF*CF;
    *BZGEO = BR*C - BT*S;
}
