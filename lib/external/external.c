#include <math.h>
#include <macro_defs.h>

#include "external.h"

void sun(const CalcTime *calc_time, double *GST, double *SLONG, double *SRASN, double *SDEC) {
    const double FDAY = (calc_time->hour*3600 + calc_time->min*60 + calc_time->sec) / 86400.;
    const double DJ = 365*(calc_time->year-1900) + (int)((calc_time->year-1901)/4) + calc_time->year_day - 0.5 + FDAY;
    const double T = DJ/36525.;
    *GST = fmod(279.690983 + 0.9856473354 * DJ + 360. * FDAY + 180., 360.) / DEF_TO_DEG;

    const double VL = fmod(279.696678 + 0.9856473354 * DJ, 360.);
    const double G = fmod(358.475845 + 0.985600267 * DJ, 360.) / DEF_TO_DEG;
    *SLONG = (VL + (1.91946 - 0.004789*T)*sin(G) + 0.020094*sin(2.*G)) / DEF_TO_DEG;
    if (*SLONG > 6.2831853) {
        *SLONG = *SLONG - 6.2831853;
    } else if (*SLONG < 0) {
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

void external_init_transformation_matrix(const CalcTime *calc_time, double G1_0, double G1_1, double H1_1, ExternalGeoGsmMatrix *mtx) {
    const double SQ = G1_1*G1_1 + H1_1*H1_1;
    const double SQQ = sqrt(SQ);
    const double SQR = sqrt(G1_0*G1_0 + SQ);

    double GST = 0., SLONG = 0., SRASN = 0., SDEC = 0.;
    if (calc_time->year > 1900 && calc_time->year < 2100) {
        sun(calc_time, &GST, &SLONG, &SRASN, &SDEC);
    }

    const double CT0 = -G1_0/SQR;
    const double CL0 = -G1_1/SQQ;
    const double SL0 = -H1_1/SQQ;
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

    mtx->A1[0] = S1*CGST + S2*SGST;
    mtx->A1[1] = -S1*SGST + S2*CGST;
    mtx->A1[2] = S3;
    mtx->A2[0] = Y1*CGST + Y2*SGST;
    mtx->A2[1] = -Y1*SGST + Y2*CGST;
    mtx->A2[2] = Y3;
    mtx->A3[0] = Z1*CGST + Z2*SGST;
    mtx->A3[1] = -Z1*SGST + Z2*CGST;
    mtx->A3[2] = Z3;
}

void external_geo_to_gsm(const ExternalGeoGsmMatrix *mtx, double geo_x, double geo_y, double geo_z, double *gsm_x, double *gsm_y, double *gsm_z) {
    const double *A1 = mtx->A1;
    const double *A2 = mtx->A2;
    const double *A3 = mtx->A3;

    const double nx = A1[0]*geo_x + A1[1]*geo_y + A1[2]*geo_z;
    const double ny = A2[0]*geo_x + A2[1]*geo_y + A2[2]*geo_z;
    const double nz = A3[0]*geo_x + A3[1]*geo_y + A3[2]*geo_z;

    *gsm_x = nx;
    *gsm_y = ny;
    *gsm_z = nz;
}

void external_gsm_to_geo(const ExternalGeoGsmMatrix *mtx, double gsm_x, double gsm_y, double gsm_z, double *geo_x, double *geo_y, double *geo_z) {
    const double *A1 = mtx->A1;
    const double *A2 = mtx->A2;
    const double *A3 = mtx->A3;
    
    const double nx = A1[0]*gsm_x + A2[0]*gsm_y + A3[0]*gsm_z;
    const double ny = A1[1]*gsm_x + A2[1]*gsm_y + A3[1]*gsm_z;
    const double nz = A1[2]*gsm_x + A2[2]*gsm_y + A3[2]*gsm_z;

    *geo_x = nx;
    *geo_y = ny;
    *geo_z = nz;
}
