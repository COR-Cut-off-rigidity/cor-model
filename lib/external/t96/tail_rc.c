#include <math.h>
#include "tail_rc.h"

void ringcurr(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    double F[2] = { 569.895366, -1603.386993 };
    double BETA[2] = { 2.722188, 3.766875 };

    double DZETASS = sqrt(Z * Z + 4.0);
    double RHOS = sqrt(X * X + Y * Y);
    double DDZETADZZ = Z / DZETASS;
    double DRHOSDX = 0.0;
    double DRHOSDY = (Y >= 0) ? 1.0 : -1.0;
    
    if (RHOS >= 1.0e-5) {
        DRHOSDX = X / RHOS;
        DRHOSDY = Y / RHOS;
    }

    double BXX = 0.0;
    double BYY = 0.0;
    double BZZ = 0.0;

    for (int I = 0; I < 2; ++I) {
        double BI = BETA[I];
        const double DZETASS_BI = DZETASS + BI;
        const double RHOS_PBI = RHOS + BI;
        const double RHOS_MBI = RHOS - BI;

        const double DZETASS_BISQ = DZETASS_BI * DZETASS_BI;
        const double BI2 = 2.0 * BI;

        const double S1 = sqrt(DZETASS_BISQ + RHOS_PBI * RHOS_PBI);
        const double S2 = sqrt(DZETASS_BISQ + RHOS_MBI * RHOS_MBI);
        const double DS1DDZ = DZETASS_BI / S1;
        const double DS2DDZ = DZETASS_BI / S2;
        const double DS1DRHOS = RHOS_PBI / S1;
        const double DS2DRHOS = RHOS_MBI / S2;

        const double S1TS2 = S1 * S2;
        const double S1PS2 = S1 + S2;
        const double S1PS2SQ = S1PS2 * S1PS2;
        const double FAC1 = sqrt(S1PS2SQ - BI2 * BI2);
        const double AS = FAC1 / (S1TS2 * S1PS2SQ);
        const double TERM1 = 1.0 / (S1TS2 * S1PS2 * FAC1);
        const double FAC2 = AS / S1PS2SQ;

        const double DASDS1 = TERM1 - FAC2 / S1 * (S2 * S2 + S1 * (3.0 * S1 + 4.0 * S2));
        const double DASDS2 = TERM1 - FAC2 / S2 * (S1 * S1 + S2 * (3.0 * S2 + 4.0 * S1));

        const double DASDX = DASDS1 * (DS1DRHOS * DRHOSDX) + DASDS2 * (DS2DRHOS * DRHOSDX);
        const double DASDY = DASDS1 * (DS1DRHOS * DRHOSDY) + DASDS2 * (DS2DRHOS * DRHOSDY);
        const double DASDZ = DASDS1 * (DS1DDZ * DDZETADZZ) + DASDS2 * (DS2DDZ * DDZETADZZ);

        BXX += -F[I] * X * DASDZ;
        BYY -= F[I] * Y * DASDZ;
        BZZ += F[I] * ((2.0 * AS + Y * DASDY) + X * DASDX);
    }

    *BX = BXX;
    *BY = BYY;
    *BZ = BZZ;
}

static void shlcar3x3(const double *A, double X, double Y, double Z, double *HX, double *HY, double *HZ) {
    double HXX = 0.0;
    double HYY = 0.0;
    double HZZ = 0.0;

    double P = A[18];
    double PP = 1.0 / (P * P);
    double YP = Y / P;
    double SYPI = sin(YP);
    double CYPI = cos(YP);

    double R = A[21];
    double ZR = Z / R;
    double SZRK = sin(ZR);
    double CZRK = cos(ZR);
    double SQPR = sqrt(PP + 1.0 / (R * R));
    double EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[0] + A[1]);
    HYY += (EPR / P * SYPI * SZRK) * (A[0] + A[1]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[0] + A[1]);

    R = A[22];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[2] + A[3]);
    HYY += (EPR / P * SYPI * SZRK) * (A[2] + A[3]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[2] + A[3]);

    R = A[23];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[4] + A[5]);
    HYY += (EPR / P * SYPI * SZRK) * (A[4] + A[5]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[4] + A[5]);

    P = A[19];
    PP = 1.0 / (P * P);
    YP = Y / P;
    SYPI = sin(YP);
    CYPI = cos(YP);

    R = A[21];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[6] + A[7]);
    HYY += (EPR / P * SYPI * SZRK) * (A[6] + A[7]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[6] + A[7]);

    R = A[22];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[8] + A[9]);
    HYY += (EPR / P * SYPI * SZRK) * (A[8] + A[9]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[8] + A[9]);

    R = A[23];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[10] + A[11]);
    HYY += (EPR / P * SYPI * SZRK) * (A[10] + A[11]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[10] + A[11]);

    P = A[20];
    PP = 1.0 / (P * P);
    YP = Y / P;
    SYPI = sin(YP);
    CYPI = cos(YP);

    R = A[21];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[12] + A[13]);
    HYY += (EPR / P * SYPI * SZRK) * (A[12] + A[13]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[12] + A[13]);

    R = A[22];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[14] + A[15]);
    HYY += (EPR / P * SYPI * SZRK) * (A[14] + A[15]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[14] + A[15]);

    R = A[23];
    ZR = Z / R;
    SZRK = sin(ZR);
    CZRK = cos(ZR);
    SQPR = sqrt(PP + 1.0 / (R * R));
    EPR = exp(X * SQPR);
    HXX += (-SQPR * EPR * CYPI * SZRK) * (A[16] + A[17]);
    HYY += (EPR / P * SYPI * SZRK) * (A[16] + A[17]);
    HZZ += (-EPR / R * CYPI * CZRK) * (A[16] + A[17]);

    *HX = HXX;
    *HY = HYY;
    *HZ = HZZ;
}

static void taildisk(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    static const double F[4] = { -745796.7338, 1176470.141, -444610.529, -57508.01028 };
    static const double BETA[4] = { 7.9250000, 8.0850000, 8.4712500, 27.89500 };

    const double Y_20 = Y / 20.0;
    const double D = 2.0 + 10.0 * (Y_20 * Y_20);
    const double DZETAS = sqrt(Z * Z + D * D);
    const double DDZETADY = (D * 10.0 * Y * 0.005) / DZETAS;
    const double DDZETADZ = Z / DZETAS;

    const double XSHIFT = X - 4.5;
    const double RHOS = sqrt(XSHIFT * XSHIFT + Y * Y);
    
    double DRHOSDX = 0.0;
    double DRHOSDY = (Y >= 0) ? 1.0 : -1.0;

    if (RHOS >= 1.0e-5) {
        DRHOSDX = XSHIFT / RHOS;
        DRHOSDY = Y / RHOS;
    }

    double BXX = 0.0;
    double BYY = 0.0;
    double BZZ = 0.0;

    for (int I = 0; I < 4; ++I) {
        const double BI = BETA[I];
        const double BI1 = 2.0 * BI;
        const double DZETAS_PBI = DZETAS + BI;
        const double RHOS_PBI = RHOS + BI;
        const double RHOS_MBI = RHOS - BI;
        const double S1 = sqrt(DZETAS_PBI * DZETAS_PBI + RHOS_PBI * RHOS_PBI);
        const double S2 = sqrt(DZETAS_PBI * DZETAS_PBI + RHOS_MBI * RHOS_MBI);
        const double DS1DDZ = DZETAS_PBI / S1;
        const double DS2DDZ = DZETAS_PBI / S2;
        const double DS1DRHOS = RHOS_PBI / S1;
        const double DS2DRHOS = RHOS_MBI / S2;

        const double S1TS2 = S1 * S2;
        const double S1PS2 = S1 + S2;
        const double S1PS2SQ = S1PS2 * S1PS2;
        const double FAC1 = sqrt(S1PS2SQ - BI1 * BI1);
        const double AS = FAC1 / (S1TS2 * S1PS2SQ);
        const double TERM1 = 1.0 / (S1TS2 * S1PS2 * FAC1);
        const double FAC2 = AS / S1PS2SQ;
        const double DASDS1 = TERM1 - FAC2 / S1 * (S2 * S2 + S1 * (3.0 * S1 + 4.0 * S2));
        const double DASDS2 = TERM1 - FAC2 / S2 * (S1 * S1 + S2 * (3.0 * S2 + 4.0 * S1));

        const double DASDX = DASDS1 * (DS1DRHOS * DRHOSDX) + DASDS2 * (DS2DRHOS * DRHOSDX);
        const double DASDY = DASDS1 * (DS1DDZ * DDZETADY + DS1DRHOS * DRHOSDY) + DASDS2 * (DS2DDZ * DDZETADY + DS2DRHOS * DRHOSDY);
        const double DASDZ = DASDS1 * (DS1DDZ * DDZETADZ) + DASDS2 * (DS2DDZ * DDZETADZ);

        BXX += F[I] * ((4.5 - X) * DASDZ);
        BYY -= F[I] * Y * DASDZ;
        BZZ += F[I] * ((2.0 * AS + Y * DASDY) + XSHIFT * DASDX);
    }

    *BX = BXX;
    *BY = BYY;
    *BZ = BZZ;
}

void tail87(double X, double Z, double *BX, double *BZ) {
    const double ZP = Z - 40.0;
    const double ZM = Z + 40.0;
    const double XNX = -10.0 - X;
    const double XNX2 = XNX * XNX;
    const double XC1 = X + 1.261;
    const double XC2 = X + 0.663;
    const double XC22 = XC2 * XC2;
    const double XR2 = XC2 * -0.1071;
    const double XC12 = XC1 * XC1;
    const double D2 = 3.0 * 3.0;
    const double B20 = Z * Z + D2;
    const double B2P = ZP * ZP + D2;
    const double B2M = ZM * ZM + D2;
    const double B = sqrt(B20);
    const double BP = sqrt(B2P);
    const double BM = sqrt(B2M);
    const double XA1 = XC12 + B20;
    const double XAP1 = XC12 + B2P;
    const double XAM1 = XC12 + B2M;
    const double XA2 = 1. / (XC22 + B20);
    const double XAP2 = 1. / (XC22 + B2P);
    const double XAM2 = 1. / (XC22 + B2M);
    const double XNA = XNX2 + B20;
    const double XNAP = XNX2 + B2P;
    const double XNAM = XNX2 + B2M;
    const double F = B20 - XC22;
    const double FP = B2P - XC22;
    const double FM = B2M - XC22;
    const double XLN1 = log(76.37 / XNA);
    const double XLNP1 = log(76.37 / XNAP);
    const double XLNM1 = log(76.37 / XNAM);
    const double XLN2 = XLN1 + 0.13238005;
    const double XLNP2 = XLNP1 + 0.13238005;
    const double XLNM2 = XLNM1 + 0.13238005;
    const double ALN = 0.25 * (XLNP1 + XLNM1 - 2. * XLN1);
    const double S0 = (atan(XNX / B) + 1.5707963) / B;
    const double S0P = (atan(XNX / BP) + 1.5707963) / BP;
    const double S0M = (atan(XNX / BM) + 1.5707963) / BM;
    const double S1 = (XLN1 * .5 + XC1 * S0) / XA1;
    const double S1P = (XLNP1 * .5 + XC1 * S0P) / XAP1;
    const double S1M = (XLNM1 * .5 + XC1 * S0M) / XAM1;
    const double S2 = (XC2 * XA2 * XLN2 + 0.1071 - F * XA2 * S0) * XA2;
    const double S2P = (XC2 * XAP2 * XLNP2 + 0.1071 - FP * XAP2 * S0P) * XAP2;
    const double S2M = (XC2 * XAM2 * XLNM2 + 0.1071 - FM * XAM2 * S0M) * XAM2;
    const double G1 = (B20 * S0 - 0.5 * XC1 * XLN1) / XA1;
    const double G1P = (B2P * S0P - 0.5 * XC1 * XLNP1) / XAP1;
    const double G1M = (B2M * S0M - 0.5 * XC1 * XLNM1) / XAM1;
    const double G2 = ((0.5 * F * XLN2 + 2.* S0 * B20 * XC2) * XA2 + XR2) * XA2;
    const double G2P = ((0.5 * FP * XLNP2 + 2. * S0P * B2P * XC2) * XAP2 + XR2) * XAP2;
    const double G2M = ((0.5 * FM * XLNM2 + 2. * S0M * B2M * XC2) * XAM2 + XR2) * XAM2;
    
    *BX = 0.391734 * (Z * S0 - 0.5 * (ZP * S0P + ZM * S0M)) + 5.89715 * (Z * S1 - 0.5 * (ZP * S1P + ZM * S1M)) + 24.6833 * (Z * S2 - 0.5 * (ZP * S2P + ZM * S2M));
    *BZ = 0.391734 * ALN + 5.89715 * (G1 - 0.5 * (G1P + G1M)) + 24.6833 * (G2 - 0.5 * (G2P + G2M));
}

void tailrc(double X, double Y, double Z, double *BXRC, double *BYRC, double *BZRC, double *BXT2, double *BYT2, double *BZT2, double *BXT3, double *BYT3, double *BZT3) {
    static const double ARC[] = {
       -3.087699646, 3.516259114, 18.81380577, -13.95772338, -5.497076303,
       0.1712890838, 2.392629189, -2.728020808, -14.79349936, 11.08738083,
       4.388174084, 0.2492163197E-01, 0.7030375685, -.7966023165, -3.835041334,
       2.642228681, -0.2405352424, -0.7297705678, 12.75434981, 11.37659788,
       636.4346279, 1.752483754, 3.604231143, 12.83078674
   };

   static const double ATAIL2[] = {
       0.8747515218, -0.9116821411, 2.209365387, -2.159059518, -7.059828867, 
       5.924671028, -1.916935691, 1.996707344, -3.877101873, 3.947666061, 
       11.38715899, -8.343210833, 1.194109867, -1.244316975, 3.73895491, 
       -4.406522465, -20.66884863, 3.020952989, 5.187666221, 6.802472048, 
       39.13543412, 2.784722096, 6.979576616, 25.71716760
   };

   static const double ATAIL3[] = {
       -19091.95061, -3011.613928, 20582.16203, 4242.918430, -2377.091102, 
       -1504.820043, 19884.04650, 2725.150544, -21389.04845, -3990.475093, 
       2401.610097, 1548.171792, -946.5493963, 490.1528941, 986.9156625, 
       -489.3265930, -67.99278499, 8.711175710, 19.69877970, 20.30095680, 
       86.45407420, 22.50403727, 23.41617329, 48.48140573
   };

   double HX, HY, HZ;
   ringcurr(X, Y, Z, &HX, &HY, &HZ);

   double WX, WY, WZ;
   shlcar3x3(ARC, X, Y, Z, &WX, &WY, &WZ);
   
   *BXRC = WX + HX;
   *BYRC = WY + HY;
   *BZRC = WZ + HZ;

   shlcar3x3(ATAIL2, X, Y, Z, &WX, &WY, &WZ);
   taildisk(X, Y, Z, &HX, &HY, &HZ);
   
   *BXT2 = WX + HX;
   *BYT2 = WY + HY;
   *BZT2 = WZ + HZ;

   shlcar3x3(ATAIL3, X, Y, Z, &WX, &WY, &WZ);
   tail87(X, Z, &HX, &HZ);
   
   *BXT3 = WX + HX;
   *BYT3 = WY;
   *BZT3 = WZ + HZ;
}
