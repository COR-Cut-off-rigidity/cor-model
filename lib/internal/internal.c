#include <math.h>
#include "internal.h"
#include "coefs/igrf_9_coefs.h"
#include "coefs/igrf_10_coefs.h"
#include "coefs/igrf_11_coefs.h"
#include "coefs/igrf_12_coefs.h"
#include "coefs/igrf_13_coefs.h"
#include "coefs/hist_coefs.h"
#include "coefs/cals10k_2_coefs.h"

void interpolate(double *G, double *H, const double *Gx, const double *Hx, const double *Gy, const double *Hy, int edge_year, int year, int iday, double years_step) {
    // TODO: interpolacia na urovni dni v roku?
    // const double F2 = (year - edge_year) / (double) years_step;
    const double F2 = (year + (iday-1) / 365.25 - edge_year) / years_step;
    const double F1 = 1. - F2;

    for (int i = 0; i < INTERNAL_COEFS_SIZE; ++i) {
        G[i] = Gx[i] * F1 + Gy[i] * F2;
        H[i] = Hx[i] * F1 + Hy[i] * F2;
    }
}

void extrapolate(double *G, double *H, const double *Gx, const double *Hx, const double *Gy, const double *Hy, int max_year, int year, int iday, double years_step) {
    const double DT = year + (iday-1) / 365.25 - (max_year-years_step);

    for (int i = 0; i < INTERNAL_COEFS_SIZE; ++i) {
        G[i] = Gx[i];
        H[i] = Hx[i];
        if(i > 45) {
            continue;
        }
        G[i] += Gy[i] * DT;
        H[i] += Hy[i] * DT;
    }
}

void sun(int iyear, int iday, int ihour, int min, int isec, double *GST, double *SLONG, double *SRASN, double *SDEC) {
    const double FDAY = (ihour*3600 + min*60 + isec) / 86400.;
    const double DJ = 365*(iyear-1900) + (int)((iyear-1901)/4) + iday - 0.5 + FDAY;
    const double T = DJ/36525.;
    *GST = fmod(279.690983 + 0.9856473354 * DJ + 360. * FDAY + 180., 360.) / DEF_TO_DEG;

    const double VL = fmod(279.696678 + 0.9856473354 * DJ, 360.);
    const double G = fmod(358.475845 + 0.985600267 * DJ, 360.) / DEF_TO_DEG;
    *SLONG = (VL + (1.91946 - 0.004789*T)*sin(G) + 0.020094*sin(2.*G)) / DEF_TO_DEG;
    if(*SLONG > 6.2831853){
        *SLONG = *SLONG - 6.2831853;
    } else if(*SLONG < 0){
        *SLONG = *SLONG + 6.2831853;
    }

    const double OBLIQ = (23.45229 - 0.0130125*T) / DEF_TO_DEG;
    const double SLP = *SLONG-9.924E-5;
    const double SOB = sin(OBLIQ);

    const double SIND = SOB*sin(SLP);
    const double COSD = sqrt(1. - SIND*SIND);
    const double SC = SIND/COSD;
    *SDEC = atan(SC);
    *SRASN = 3.141592654 - atan2(cos(OBLIQ)/SOB*SC,-cos(SLP)/COSD);
}

int interpolate_or_extrapolate(const InternalModelCoefs *coefs, int year, int iday, double *G, double *H) {
    const double **GHs = coefs->coefs;

    if(coefs->extrapolate && coefs->max_year-year <= coefs->years_step) {
        extrapolate(G, H, GHs[coefs->coefs_size-4], GHs[coefs->coefs_size-3], GHs[coefs->coefs_size-2], GHs[coefs->coefs_size-1], coefs->max_year, year, iday, coefs->years_step);
    } else {
        if(year >= coefs->max_year) {
            return -1;
        }

        int GH_year;
        if(year <= coefs->min_year) {
            GH_year = coefs->max_year / coefs->years_step - (coefs->min_year + 1) / coefs->years_step;
        } else {
            GH_year = coefs->max_year / coefs->years_step - year / coefs->years_step;
        }

        const int GH_position = coefs->coefs_size - 2*GH_year - 1;
        const int edge_year = coefs->min_year + ((GH_position - 1) / 2) * coefs->years_step;
        interpolate(G, H, GHs[GH_position-1], GHs[GH_position], GHs[GH_position+1], GHs[GH_position+2], edge_year, year, iday, coefs->years_step);
    }

    return 0;
}

int internal_recalc(const InternalModelCoefs *coefs, const RecalcTime *time, double *G, double *H, double *REC, double *A1, double *A2, double *A3) {
    if(time->iy > coefs->max_year || time->iy < coefs->min_year) {
        return -1;
    }

    double N2;
    for (int i = 1; i < 15; ++i) {
        N2 = 2*i - 1;
        N2 = N2 * (N2-2);
        int MN = 0, j;
        for (j = 1; j < i+1; ++j) {
            MN = i*(i-1) / 2+j-1;
            REC[MN] = ((i-j)*(i+j-2)) / N2;
        }
    }

    if(interpolate_or_extrapolate(coefs, time->iy, time->id, G, H) == -1) {
        return -1;
    }

    double S = 1., P, AA;
    int MN = 0, MNN;
    for (int i = 2; i <= 14; ++i) {
        MN = (i*(i-1)/2+1)-1;
        S = S * (2*i-3) / (i-1);
        G[MN] *= S;
        H[MN] *= S;
        P = S;
        for (int j = 2; j <= i; j++){
            AA = 1.;
            if (j == 2){
                AA = 2.;
            }
            P = P * sqrt(AA * (i-j+1) / (i+j-2));
            MNN = MN+j-1;
            G[MNN] *= P;
            H[MNN] *= P;
        }
    }

    const double G10 = -G[1];
    const double G11 = G[2];
    const double H11 = H[2];

    const double SQ = G11*G11 + H11*H11;
    const double SQQ = sqrt(SQ);
    const double SQR = sqrt(G10*G10 + SQ);

    double GST = 0., SLONG = 0., SRASN = 0., SDEC = 0.;
    if(time->iy > 1900 && time->iy < 2100) {
        sun(time->iy, time->id, time->ih, time->min, time->is, &GST, &SLONG, &SRASN, &SDEC);
    }

    const double CT0 = G10/SQR;
    const double CL0 = -G11/SQQ;
    const double SL0 = -H11/SQQ;
    const double ST0 = SQQ/SQR;
    const double STCL = ST0*CL0;
    const double STSL = ST0*SL0;
    const double CGST = cos(GST);
    const double SGST = sin(GST);
    const double DIP1 = STCL*CGST - STSL*SGST;
    const double DIP2 = STCL*SGST + STSL*CGST;
    const double DIP3 = CT0;

    const double CSDEC = cos(SDEC);

    const double S1 = cos(SRASN) * CSDEC;
    const double S2 = sin(SRASN) * CSDEC;
    const double S3 = sin(SDEC);
    double Y1 = DIP2*S3 - DIP3*S2;
    double Y2 = DIP3*S1 - DIP1*S3;
    double Y3 = DIP1*S2 - DIP2*S1;
    const double Y = sqrt(Y1*Y1 + Y2*Y2 + Y3*Y3);
    Y1 = Y1/Y;
    Y2 = Y2/Y;
    Y3 = Y3/Y;
    const double Z1 = S2*Y3 - S3*Y2;
    const double Z2 = S3*Y1 - S1*Y3;
    const double Z3 = S1*Y2 - S2*Y1;

    A1[0] = S1*CGST + S2*SGST;
    A1[1] = -S1*SGST + S2*CGST;
    A1[2] = S3;
    A2[0] = Y1*CGST + Y2*SGST;
    A2[1] = -Y1*SGST + Y2*CGST;
    A2[2] = Y3;
    A3[0] = Z1*CGST + Z2*SGST;
    A3[1] = -Z1*SGST + Z2*CGST;
    A3[2] = Z3;
    return 0;
}

void internal_geo(double XGEO, double YGEO, double ZGEO, const double *G, const double *H, const double *REC, double *BXGEO, double *BYGEO, double *BZGEO) {
    double A[14], B[14];
    const double RHO2 = XGEO*XGEO + YGEO*YGEO;
    const double R = sqrt(RHO2 + ZGEO*ZGEO);
    const double C = ZGEO/R;
    const double RHO = sqrt(RHO2);
    const double S = RHO/R;
    double CF, SF;
    if(S < 1.E-5){
        CF = 1.;
        SF = 0.;
    } else {
        CF = XGEO/RHO;
        SF = YGEO/RHO;
    }

    const double PP = 1./R;
    double P = PP;

    const int IRP3 = (int) R + 2;
    int NM = 3 + 30/IRP3;
    if(NM > 13){
        NM = 13;
    }
    const int K = NM+1;
    for(int i = 1; i <= K; i++) {
        P = P * PP;
        A[i-1] = P;
        B[i-1] = P * i;
    }
    P = 1.;

    int MM = 0, MN;
    double D = 0., BBR = 0., BBT = 0., BBF = 0., X = 0., Y = 1., W, Q, Z, BI, P2, D2, AN, E, HH, QQ, XK, DP, PM;
    for(int i = 1; i <= K; i++){
        if(i != 1){
            MM = i - 1;
            W = X;
            X = W*CF + Y*SF;
            Y = Y*CF - W*SF;
        }
        Q = P;
        Z = D;
        BI = 0.;
        P2 = 0.;
        D2 = 0.;

        for(int j = i; j <= K; j++){
            AN = A[j-1];
            MN = j*(j - 1)/2 + i;
            E = G[MN-1];
            HH = H[MN-1];
            W = E*Y + HH*X;
            BBR = BBR + B[j-1]*W*Q;
            BBT = BBT - AN*W*Z;
            if(i != 1){
                QQ = Q;
                if(S < 1.E-5) {
                    QQ = Z;
                }
                BI = BI + AN*(E*X - HH*Y)*QQ;
            }
            XK = REC[MN-1];
            DP = C*Z - S*Q - XK*D2;
            PM = C*Q - XK*P2;
            D2 = Z;
            P2 = Q;
            Z = DP;
            Q = PM;
        }
        D = S*D + C*P;
        P = S*P;
        if(i != 1){
            BI = BI*MM;
            BBF = BBF + BI;
        }
    }

    const double BR = BBR;
    const double BT = BBT;
    double BF;
    if(S < 1.E-5){
        if(C < 0.){
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

// recalc -> lines 110 - 128 must be commented out
/*void internal_geo2(double XGEO, double YGEO, double ZGEO, const double *G, const double *H, double *BXGEO, double *BYGEO, double *BZGEO) {
    double RHO2 = XGEO * XGEO + YGEO * YGEO;
    double R = sqrt(RHO2 + ZGEO * ZGEO);
    const double C = ZGEO / R;

    double r = R * DEF_RE_KM_DF, theta = 0., phi = 0., S = 0., CF = 1., SF = 0.;
    
    if (RHO2 != 0.) {
        const double RHO = sqrt(RHO2);
        S = RHO / R;
        CF = XGEO/RHO;
        SF = YGEO/RHO;
        phi = atan2(YGEO, XGEO);
        theta = atan2(RHO, ZGEO);
        
        if (phi < 0.) {
            phi += DEF_2PI_DF;
        }
    } else if (ZGEO < 0.) {
        theta = DEF_PI_DF;
    }
    
    double B_r = 0.;
    double B_theta = 0.;
    double B_phi = 0.;

    double tem = DEF_RE_KM_DF / r;
    double tem_r[NM_MAX];
    tem_r[0] = tem;

    double dtem_Rr[NM_MAX];
    dtem_Rr[0] = -tem/r;

    for (int i = 1; i < NM_MAX; i++) {
        double tem_r_new = tem_r[i - 1] * tem;
        tem_r[i] = tem_r_new;
        dtem_Rr[i] = -(i + 1) * tem_r_new / r;
    }

    double pmn2[NM_MAX][NM_MAX];
    double dpmn2[NM_MAX][NM_MAX];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    pmn2[0][0] = 1.;
    pmn2[1][0] = cos_theta;
    pmn2[1][1] = sin_theta;

    dpmn2[0][0] = 0.;
    dpmn2[1][0] = -sin_theta;
    dpmn2[1][1] = cos_theta;

    for (int n = 2; n < NM_MAX; n++) {
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

    int n = 0;
    int m = 0;
    int n_max = 1;

    for (int mm = 0; mm < INTERNAL_COEFS_SIZE; mm++) {
        if(m == n_max) {
            n++;
            n_max++;
            m = 0;
        }

        double m_phi = m*phi;
        double sin_m_phi = sin(m_phi);
        double cos_m_phi = cos(m_phi);
        double temp_G = G[mm];
        double temp_H = H[mm];
        double tem = temp_G * cos_m_phi + temp_H * sin_m_phi;
        double temp_pmn2 = pmn2[n][m];
        double temp_tem_r = tem_r[n];

        B_r += tem * temp_pmn2 * dtem_Rr[n];
        B_theta += tem * dpmn2[n][m] * temp_tem_r;
        B_phi += (-temp_G * m * sin_m_phi + temp_H * m * cos_m_phi) * temp_pmn2 * temp_tem_r;
        m++;
    }

    B_theta = -B_theta * DEF_RE_KM_DF / r;
    B_phi = -B_phi * DEF_RE_KM_DF / (r * sin_theta);
    B_r = -B_r * DEF_RE_KM_DF;

    const double HE = B_r*S + B_theta*C;
    *BXGEO = HE*CF - B_phi*SF;
    *BYGEO = HE*SF + B_phi*CF;
    *BZGEO = B_r * C - B_theta*S;
}*/