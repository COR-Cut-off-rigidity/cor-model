#include <math.h>
#include "t96.h"

void t96_gsm(const double *PARMOD, double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double PDYN = PARMOD[0];
    const double DEPR = 0.8 * PARMOD[1] - 13.0 * sqrt(PDYN);
    const double BYIMF = PARMOD[2];
    const double BZIMF = PARMOD[3];
    const double Bt = sqrt(BYIMF * BYIMF + BZIMF * BZIMF);

    double THETA = 0.0;
    if (BYIMF != 0.0 || BZIMF != 0.0) {
        THETA = atan2(BYIMF, BZIMF);
        if (THETA <= 0.0) {
            THETA += 6.2831853;
        }
    }

    const double CT = cos(THETA);
    const double ST = sin(THETA);
    const double EPS = 718.5 * sqrt(PDYN) * Bt * sin(THETA / 2.);
    const double YS = Y * CT - Z * ST;
    const double ZS = Z * CT + Y * ST;
    const double YS_DELIMFY = YS / 10.0;
    const double FACTIMF = exp(X / 20.0 - YS_DELIMFY * YS_DELIMFY);
    const double OIMFY = 0.7850 * BYIMF * FACTIMF;
    const double OIMFZ = 0.7850 * BZIMF * FACTIMF;
    const double XAPPA = pow(PDYN / 2.0, 0.14);
    const double XAPPA3 = XAPPA * XAPPA * XAPPA;
    const double AM = 70.0 / XAPPA;
    const double RHO2 = Y * Y + Z * Z;
    const double ASQ = AM * AM;

    double XMXM = AM + X - (5.48 / XAPPA);
    if (XMXM < 0.0) {
        XMXM = 0.0;
    }

    const double AXX0 = XMXM * XMXM;
    const double ARO = ASQ + RHO2;
    const double ARO_AXX0 = ARO + AXX0;
    const double SIGMA = sqrt((ARO + AXX0 + sqrt(ARO_AXX0 * ARO_AXX0 - 4.0 * ASQ * AXX0)) / (2.0 * ASQ));

    if (SIGMA >= 1.08 + 0.005) {
        double QX, QY, QZ;
        dipole_t96(X, Y, Z, &QX, &QY, &QZ);

        *BX = -QX;
        *BY = OIMFY - QY;
        *BZ = OIMFZ - QZ;
        return;
    }

    const double FACTEPS = EPS / 3630.7 - 1.0;
    const double FACTPD = sqrt(PDYN / 2.0) - 1.0;
    const double RCAMPL = -1.162 * DEPR;
    const double TAMPL2 = 22.344 + 18.50 * FACTPD + 2.602 * FACTEPS;
    const double TAMPL3 = 6.903 + 5.287 * FACTPD;
    const double B1AMPL = 0.5790 + 0.4462 * FACTEPS;
    const double B2AMPL = 20.0 * B1AMPL;
    const double XX = X * XAPPA;
    const double YY = Y * XAPPA;
    const double ZZ = Z * XAPPA;

    double CFX, CFY, CFZ;
    cylharm(XX, YY, ZZ, &CFX, &CFY, &CFZ);

    double BXRC, BYRC, BZRC, BXT2, BYT2, BZT2, BXT3, BYT3, BZT3;
    tailrc_96(XX, YY, ZZ, &BXRC, &BYRC, &BZRC, &BXT2, &BYT2, &BZT2, &BXT3, &BYT3, &BZT3);

    double R1X, R1Y, R1Z;
    birk1tot_02(XX, YY, ZZ, &R1X, &R1Y, &R1Z);

    double R2X, R2Y, R2Z;
    birk2shl(XX, YY, ZZ, &R2X, &R2Y, &R2Z);

    double R2X_HX, R2Y_HY, R2Z_HZ;
    r2_birk(XX, YY, ZZ, &R2X_HX, &R2Y_HY, &R2Z_HZ);

    R2X += R2X_HX;
    R2Y += R2Y_HY;
    R2Z += R2Z_HZ;

    double RIMFX, RIMFYS, RIMFZS;
    intercon(XX, YS * XAPPA, ZS * XAPPA, &RIMFX, &RIMFYS, &RIMFZS);

    const double RIMFY = RIMFYS * CT + RIMFZS * ST;
    const double RIMFZ = RIMFZS * CT - RIMFYS * ST;
    const double RIMFAMPL = 0.7850 * Bt;

    const double FX = CFX * XAPPA3 + RCAMPL * BXRC + TAMPL2 * BXT2 + TAMPL3 * BXT3 + B1AMPL * R1X + B2AMPL * R2X + RIMFAMPL * RIMFX;
    const double FY = CFY * XAPPA3 + RCAMPL * BYRC + TAMPL2 * BYT2 + TAMPL3 * BYT3 + B1AMPL * R1Y + B2AMPL * R2Y + RIMFAMPL * RIMFY;
    const double FZ = CFZ * XAPPA3 + RCAMPL * BZRC + TAMPL2 * BZT2 + TAMPL3 * BZT3 + B1AMPL * R1Z + B2AMPL * R2Z + RIMFAMPL * RIMFZ;

    if (SIGMA < 1.08 - 0.005) {
        *BX = FX;
        *BY = FY;
        *BZ = FZ;
        return;
    }

    const double FINT = 0.5 * (1.0 - (SIGMA - 1.08) / 0.005);
    const double FEXT = 0.5 * (1.0 + (SIGMA - 1.08) / 0.005);
    
    double QX, QY, QZ;
    dipole_t96(X, Y, Z, &QX, &QY, &QZ);
    
    *BX = (FX + QX) * FINT - QX;
    *BY = (FY + QY) * FINT + OIMFY * FEXT - QY;
    *BZ = (FZ + QZ) * FINT + OIMFZ * FEXT - QZ;
}

void cylharm(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    static const double A[12] = { 0.24777, -27.003, -0.46815, 7.0637, -1.5918, -0.90317E-01, 57.522, 13.757, 2.0100, 10.458, 4.5798, 2.1695 };
    double RHO = sqrt(Y * Y + Z * Z);
    double SINFI, COSFI;

    if (RHO < 1.0e-8) {
        SINFI = 1.0;
        COSFI = 0.0;
        RHO = 1.0e-8;
    } else {
        SINFI = Z / RHO;
        COSFI = Y / RHO;
    }
    
    double SINFI2 = SINFI * SINFI;
    double SI2CO2 = SINFI2 - COSFI * COSFI;
    
    double BXX = 0.0;
    double BYY = 0.0;
    double BZZ = 0.0;
    
    for (int I = 0; I < 3; ++I) {
        double DZETA = RHO / A[I + 6];
        double XJ0 = bes0(DZETA);
        double XJ1 = bes1(DZETA);
        double XEXP = exp(X / A[I + 6]);
        
        BXX -= A[I] * XJ1 * XEXP * SINFI;
        BYY += A[I] * (2.0 * XJ1 / DZETA - XJ0) * XEXP * SINFI * COSFI;
        BZZ += A[I] * (XJ1 / DZETA * SI2CO2 - XJ0 * SINFI2) * XEXP;
    }
    
    for (int I = 3; I < 6; ++I) {
        double DZETA = RHO / A[I + 6];
        double XKSI = X / A[I + 6];
        double XJ0 = bes0(DZETA);
        double XJ1 = bes1(DZETA);
        double XEXP = exp(XKSI);
        
        double BRHO = (XKSI * XJ0 - (DZETA * DZETA + XKSI - 1.0) * XJ1 / DZETA) * XEXP * SINFI;
        double BPHI = (XJ0 + XJ1 / DZETA * (XKSI - 1.0)) * XEXP * COSFI;
        
        BXX += A[I] * (DZETA * XJ0 + XKSI * XJ1) * XEXP * SINFI;
        BYY += A[I] * (BRHO * COSFI - BPHI * SINFI);
        BZZ += A[I] * (BRHO * SINFI + BPHI * COSFI);
    }

    *BX = BXX;
    *BY = BYY;
    *BZ = BZZ;
}

void cylhar1(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    static const double A[12] = { -0.65385, -18.061, -0.40457, -5.0995, 1.2846, 0.78231E-01, 39.592, 13.291, 1.9970, 10.062, 4.5140, 2.1558 };
    double RHO = sqrt(Y * Y + Z * Z);
    double SINFI = 1.0;
    double COSFI = 0.0;
    
    if (RHO >= 1.0e-10) {
        SINFI = Z / RHO;
        COSFI = Y / RHO;
    }
    
    double BXX = 0.0;
    double BYY = 0.0;
    double BZZ = 0.0;
    
    for (int I = 0; I < 3; ++I) {
        double DZETA = RHO / A[I + 6];
        double XKSI = X / A[I + 6];
        double XJ0 = bes0(DZETA);
        double XJ1 = bes1(DZETA);
        double XEXP = exp(XKSI);
        double BRHO = XJ1 * XEXP;
        
        BXX -= A[I] * XJ0 * XEXP;
        BYY += A[I] * BRHO * COSFI;
        BZZ += A[I] * BRHO * SINFI;
    }
    
    for (int I = 3; I < 6; ++I) {
        double DZETA = RHO / A[I + 6];
        double XKSI = X / A[I + 6];
        double XJ0 = bes0(DZETA);
        double XJ1 = bes1(DZETA);
        double XEXP = exp(XKSI);
        double BRHO = (DZETA * XJ0 + XKSI * XJ1) * XEXP;
        
        BXX += A[I] * (DZETA * XJ1 - XJ0 * (XKSI + 1.0)) * XEXP;
        BYY += A[I] * BRHO * COSFI;
        BZZ += A[I] * BRHO * SINFI;
    }

    *BX = BXX;
    *BY = BYY;
    *BZ = BZZ;
}

double bes0(double X) {
    if (fabs(X) < 3.0) {
        double X32 = (X / 3.0) * (X / 3.0);
        return 1.0 - X32 * (2.2499997 - X32 * (1.2656208 - X32 * (0.3163866 - X32 * (0.0444479 - X32 * (0.0039444 - X32 * 0.00021)))));
    }

    double XD3 = 3.0 / X;
    double F0 = 0.79788456 - XD3 * (0.00000077 + XD3 * (0.00552740 + XD3 * (0.00009512 - XD3 * (0.00137237 - XD3 * (0.00072805 - XD3 * 0.00014476)))));
    double T0 = X - 0.78539816 - XD3 * (0.04166397 + XD3 * (0.00003954 - XD3 * (0.00262573 - XD3 * (0.00054125 + XD3 * (0.00029333 - XD3 * 0.00013558)))));
    return F0 / sqrt(X) * cos(T0);
}

double bes1(double X) {
    if (fabs(X) < 3.0) {
        double X32 = (X / 3.0) * (X / 3.0);
        double BES1XM1 = 0.5 - X32 * (0.56249985 - X32 * (0.21093573 - X32 * (0.03954289 - X32 * (0.00443319 - X32 * (0.00031761 - X32 * 0.00001109)))));
        return BES1XM1 * X;
    }

    double XD3 = 3.0 / X;
    double F1 = 0.79788456 + XD3 * (0.00000156 + XD3 * (0.01659667 + XD3 * (0.00017105 - XD3 * (0.00249511 - XD3 * (0.00113653 - XD3 * 0.00020033)))));
    double T1 = X - 2.35619449 + XD3 * (0.12499612 + XD3 * (0.0000565 - XD3 * (0.00637879 - XD3 * (0.00074348 + XD3 * (0.00079824 - XD3 * 0.00029166)))));
    return F1 / sqrt(X) * cos(T1);
}

void intercon(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    static const double A[15] = {
        -8.411078731, 5932254.951, -9073284.93, -11.68794634,
        6027598.824, -9218378.368, -6.508798398, -11824.42793, 18015.66212,
        7.99754043, 13.9669886, 90.24475036, 16.75728834, 1015.645781,
        1553.493216
    };

    static const double RP[3] = { 0.12503844191332988, 0.071597395204476152, 0.011080976889832743 };
    static const double RR[3] = { 0.059675527522452816, 0.00098459526668325296, 0.00064371058922777988 };

    double BXX = 0.0;
    double BYY = 0.0;
    double BZZ = 0.0;
    int L = 0;

    for (int i = 0; i < 3; i++) {
        double CYPI = cos(Y * RP[i]);
        double SYPI = sin(Y * RP[i]);

        for (int k = 0; k < 3; k++) {
            double SZRK = sin(Z * RR[k]);
            double CZRK = cos(Z * RR[k]);
            double SQPR = sqrt(RP[i] * RP[i] + RR[k] * RR[k]);
            double EPR = exp(X * SQPR);
            BXX += A[L] * (-SQPR * EPR * CYPI * SZRK);
            BYY += A[L] * (RP[i] * EPR * SYPI * SZRK);
            BZZ += A[L] * (-RR[k] * EPR * CYPI * CZRK);
            L++;
        }
    }

    *BX = BXX;
    *BY = BYY;
    *BZ = BZZ;
}

void tailrc_96(double X, double Y, double Z, double *BXRC, double *BYRC, double *BZRC, double *BXT2, double *BYT2, double *BZT2, double *BXT3, double *BYT3, double *BZT3) {
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
    ringcurr_96(X, Y, Z, &HX, &HY, &HZ);

    double WX, WY, WZ;
    shlcar3x3_t96(ARC, X, Y, Z, &WX, &WY, &WZ);
    
    *BXRC = WX + HX;
    *BYRC = WY + HY;
    *BZRC = WZ + HZ;

    shlcar3x3_t96(ATAIL2, X, Y, Z, &WX, &WY, &WZ);
    taildisk_t96(X, Y, Z, &HX, &HY, &HZ);
    
    *BXT2 = WX + HX;
    *BYT2 = WY + HY;
    *BZT2 = WZ + HZ;

    shlcar3x3_t96(ATAIL3, X, Y, Z, &WX, &WY, &WZ);
    tail87(X, Z, &HX, &HZ);
    
    *BXT3 = WX + HX;
    *BYT3 = WY;
    *BZT3 = WZ + HZ;
}

void ringcurr_96(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
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

void taildisk_t96(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
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

void shlcar3x3_t96(const double *A, double X, double Y, double Z, double *HX, double *HY, double *HZ) {
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

void birk1tot_02(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    static const double C1[14] = {
        -0.911582E-03, -0.376654E-02, -0.727423E-02, -0.270084E-02,
        -0.123899E-02, -0.154387E-02, -0.340040E-02, -0.191858E-01,
        -0.518979E-01, 0.635061E-01, 0.440680, -0.396570, -89.5587,
        23.2806
    };

    static const double C2[] = {
        6.04133, 0.305415, 0.606066E-02, 0.128379E-03, -0.179406E-04,
        1.41714, -27.2586, -4.28833, -1.30675, 35.5607,
        8.95792, 0.961617E-03, -0.801477E-03, -0.782795E-03, -1.65242,
        -16.5242, -5.33798, 0.424878E-03, 0.331787E-03, -0.704305E-03,
        0.844342E-03, 0.953682E-04, 0.886271E-03, 25.1120, 20.9299,
        5.14569, -44.1670, -51.0672, -1.87725, 20.2998,
        48.7505, -2.97415, 0.174615E-03, -0.789777E-06, 0.686047E-03,
        0.460104E-04, -0.345216E-02, 0.221871E-02, 0.110078E-01, -0.661373E-02,
        0.249201E-02, 0.343978E-01
    };

    const double TNOONN = 12.0 * 0.01745329;
    const double TNOONS = 3.141592654 - TNOONN;
    const double R2 = X * X + Y * Y + Z * Z;
    const double R = sqrt(R2);
    const double R3 = R * R2;

    double PAS = 0.0;
    if (X != 0.0 || Y != 0.0) {
        PAS = atan2(Y, X);
    }

    const double TAS = atan2(sqrt(X * X + Y * Y), Z);
    const double STAS = sin(TAS);

    double TET0 = asin(STAS / pow((STAS * STAS * STAS * STAS * STAS * STAS) * (1.0 - R3) + R3, 0.1666666667));
    if (TAS > 1.5707963) {
        TET0 = 3.141592654 - TET0;
    }

    const double DTET = (8.0 * 0.01745329) * (sin(PAS * 0.5) * sin(PAS * 0.5));
    const double TETR1N = TNOONN + DTET;
    const double TETR1S = TNOONS - DTET;

    double BXX = 0.0;
    double BYY = 0.0;
    double BZZ = 0.0;
    double D1[3][14], D2[3][42];

    const double T01N = TETR1N - 0.034906;
    const double T01S = TETR1S - 0.034906;
    const double T02S = TETR1S + 0.034906;
    const double T02N = TETR1N + 0.034906;

    if (TET0 < T01N || TET0 > T02S) {
        diploop1(X, Y, Z, D1);
        for (int I = 0; I < 14; I++) {
            BXX += C1[I] * D1[0][I];
            BYY += C1[I] * D1[1][I];
            BZZ += C1[I] * D1[2][I];
        }
    }

    if (TET0 > T02N && TET0 < T01S) {
        condip1(X, Y, Z, D2);
        for (int I = 0; I < 42; I++) {
            BXX += C2[I] * D2[0][I];
            BYY += C2[I] * D2[1][I];
            BZZ += C2[I] * D2[2][I];
        }
    }

    if (TET0 >= T01N && TET0 <= T02N) {
        double SQR = sqrt(R);
        
        double SIN_T01 = sin(T01N);
        double ST01AS = SQR / pow(R3 + 1.0 / (SIN_T01 * SIN_T01 * SIN_T01 * SIN_T01 * SIN_T01 * SIN_T01) - 1.0, 0.1666666667);
        double SIN_T02 = sin(T02N);
        double ST02AS = SQR / pow(R3 + 1.0 / (SIN_T02 * SIN_T02 * SIN_T02 * SIN_T02 * SIN_T02 * SIN_T02) - 1.0, 0.1666666667);

        double CT01AS = sqrt(1.0 - ST01AS*ST01AS);
        double CT02AS = sqrt(1.0 - ST02AS*ST02AS);
        double XAS1 = R * ST01AS * cos(PAS);
        double Y1 = R * ST01AS * sin(PAS);
        double ZAS1 = R * CT01AS;
        double X1 = XAS1;
        double Z1 = ZAS1;
        diploop1(X1, Y1, Z1, D1);

        double BX1 = 0.0;
        double BY1 = 0.0;
        double BZ1 = 0.0;
        for (int I = 0; I < 14; I++) {
            BX1 += C1[I] * D1[0][I];
            BY1 += C1[I] * D1[1][I];
            BZ1 += C1[I] * D1[2][I];
        }

        double XAS2 = R * ST02AS * cos(PAS);
        double Y2 = R * ST02AS * sin(PAS);
        double ZAS2 = R * CT02AS;
        double X2 = XAS2;
        double Z2 = ZAS2;

        condip1(X2, Y2, Z2, D2);
        double BX2 = 0.0;
        double BY2 = 0.0;
        double BZ2 = 0.0;
        for (int I = 0; I < 42; I++) {
            BX2 += C2[I] * D2[0][I];
            BY2 += C2[I] * D2[1][I];
            BZ2 += C2[I] * D2[2][I];
        }

        double X_DIFF2 = X2 - X1;
        double Y_DIFF2 = Y2 - Y1;
        double Z_DIFF2 = Z2 - Z1;
        double SS = sqrt(X_DIFF2 * X_DIFF2 + Y_DIFF2 * Y_DIFF2 + Z_DIFF2 * Z_DIFF2);
        double X_DIFF = X - X1;
        double Y_DIFF = Y - Y1;
        double Z_DIFF = Z - Z1;
        double DS = sqrt(X_DIFF * X_DIFF + Y_DIFF * Y_DIFF + Z_DIFF * Z_DIFF);
        double FRAC = DS / SS;
        BXX = BX1 * (1.0 - FRAC) + BX2 * FRAC;
        BYY = BY1 * (1.0 - FRAC) + BY2 * FRAC;
        BZZ = BZ1 * (1.0 - FRAC) + BZ2 * FRAC;
    }

    if (TET0 >= T01S && TET0 <= T02S) {
        double SQR = sqrt(R);

        double SIN_T01 = sin(T01S);
        double ST01AS = SQR / pow(R3 + 1.0 / (SIN_T01 * SIN_T01 * SIN_T01 * SIN_T01 * SIN_T01 * SIN_T01) - 1.0, 0.1666666667);
        double SIN_T02 = sin(T02S);
        double ST02AS = SQR / pow(R3 + 1.0 / (SIN_T02 * SIN_T02 * SIN_T02 * SIN_T02 * SIN_T02 * SIN_T02) - 1.0, 0.1666666667);

        double CT01AS = -sqrt(1.0 - ST01AS*ST01AS);
        double CT02AS = -sqrt(1.0 - ST02AS*ST02AS);
        double XAS1 = R * ST01AS * cos(PAS);
        double Y1 = R * ST01AS * sin(PAS);
        double ZAS1 = R * CT01AS;
        double X1 = XAS1;
        double Z1 = ZAS1;

        condip1(X1, Y1, Z1, D2);
        double BX1 = 0.0;
        double BY1 = 0.0;
        double BZ1 = 0.0;
        for (int I = 0; I < 42; I++) {
            BX1 += C2[I] * D2[0][I];
            BY1 += C2[I] * D2[1][I];
            BZ1 += C2[I] * D2[2][I];
        }

        double XAS2 = R * ST02AS * cos(PAS);
        double Y2 = R * ST02AS * sin(PAS);
        double ZAS2 = R * CT02AS;
        double X2 = XAS2;
        double Z2 = ZAS2;

        diploop1(X2, Y2, Z2, D1);
        double BX2 = 0.0;
        double BY2 = 0.0;
        double BZ2 = 0.0;
        for (int I = 0; I < 14; I++) {
            BX2 += C1[I] * D1[0][I];
            BY2 += C1[I] * D1[1][I];
            BZ2 += C1[I] * D1[2][I];
        }

        double X_DIFF2 = X2 - X1;
        double Y_DIFF2 = Y2 - Y1;
        double Z_DIFF2 = Z2 - Z1;
        double SS = sqrt(X_DIFF2 * X_DIFF2 + Y_DIFF2 * Y_DIFF2 + Z_DIFF2 * Z_DIFF2);
        double X_DIFF = X - X1;
        double Y_DIFF = Y - Y1;
        double Z_DIFF = Z - Z1;
        double DS = sqrt(X_DIFF * X_DIFF + Y_DIFF * Y_DIFF + Z_DIFF * Z_DIFF);
        double FRAC = DS / SS;
        BXX = BX1 * (1.0 - FRAC) + BX2 * FRAC;
        BYY = BY1 * (1.0 - FRAC) + BY2 * FRAC;
        BZZ = BZ1 * (1.0 - FRAC) + BZ2 * FRAC;
    }

    double BSX, BSY, BSZ;
    birk1shld(X, Y, Z, &BSX, &BSY, &BSZ);
    *BX = BXX + BSX;
    *BY = BYY + BSY;
    *BZ = BZZ + BSZ;
}

void diploop1(double X, double Y, double Z, double D[3][14]) {
    static const double XX[12] = { -11.0, -7.0, -7.0, -3.0, -3.0, 1.0, 1.0, 1.0, 5.0, 5.0, 9.0, 9.0 };
    static const double YY[12] = { 2.0, 0.0, 4.0, 2.0, 6.0, 0.0, 4.0, 8.0, 2.0, 6.0, 0.0, 4.0 };
    const double Z2 = Z * Z;

    for (int I = 0; I < 12; I++) {
        double XX_DIPY = YY[I] * 0.945719;
        double XD = XX[I] * 1.12541;

        const double XXD = X - XD;
        double YXX_DIPY = Y - XX_DIPY;
        double R2 = XXD * XXD + YXX_DIPY * YXX_DIPY + Z2;
        double XMR5 = 30574.0 / (R2 * R2 * sqrt(R2));
        double XMR53 = 3.0 * XMR5;

        const double BX1Z = XMR53 * XXD * Z;
        const double BY1Z = XMR53 * YXX_DIPY * Z;
        const double BZ1Z = XMR5 * (3.0 * Z2 - R2);
        
        double BX2Z = 0.0, BY2Z = 0.0, BZ2Z = 0.0;
        if(fabs(XX_DIPY) > 1.0e-10) {
            YXX_DIPY = Y + XX_DIPY;
            R2 = XXD * XXD + YXX_DIPY * YXX_DIPY + Z2;
            XMR5 = 30574.0 / (R2 * R2 * sqrt(R2));
            XMR53 = 3.0 * XMR5;
            BX2Z = XMR53 * XXD * Z;
            BY2Z = XMR53 * YXX_DIPY * Z;
            BZ2Z = XMR5 * (3.0 * Z2 - R2);
        }

        D[0][I] = BX1Z + BX2Z;
        D[1][I] = BY1Z + BY2Z;
        D[2][I] = BZ1Z + BZ2Z;
    }

    double BXOCT1, BYOCT1, BZOCT1;
    crosslp(X, Y, Z, &BXOCT1, &BYOCT1, &BZOCT1, 2.28397, 1.86106, 1.00891);
    D[0][12] = BXOCT1;
    D[1][12] = BYOCT1;
    D[2][12] = BZOCT1;

    double BX, BY, BZ;
    circle(X + 5.60831, Y, Z, 7.83281, &BX, &BY, &BZ);
    D[0][13] = BX;
    D[1][13] = BY;
    D[2][13] = BZ;
}

void circle(double X, double Y, double Z, double RL, double *BX, double *BY, double *BZ) {
    double RHO2 = X * X + Y * Y;
    double RHO = sqrt(RHO2);
    double RHO_RL = RHO + RL;
    double R22 = Z * Z + RHO_RL * RHO_RL;
    double R2 = sqrt(R22);
    double R12 = R22 - 4.0 * RHO * RL;
    double R32 = 0.5 * (R12 + R22);
    double XK2 = 1.0 - R12 / R22;
    double XK2S = 1.0 - XK2;
    double DL = log(1.0 / XK2S);
    double K = 1.38629436112 + XK2S * (0.09666344259 + XK2S * (0.03590092383 +
         XK2S * (0.03742563713 + XK2S * 0.01451196212))) + DL *
         (0.5 + XK2S * (0.12498593597 + XK2S * (0.06880248576 +
         XK2S * (0.03328355346 + XK2S * 0.00441787012))));
    double E = 1.0 + XK2S * (0.44325141463 + XK2S * (0.0626060122 + XK2S *
         (0.04757383546 + XK2S * 0.01736506451))) + DL *
         XK2S * (0.2499836831 + XK2S * (0.09200180037 + XK2S *
         (0.04069697526 + XK2S * 0.00526449639)));

    double BRHO;
    if (RHO > 1e-6) {
        BRHO = Z / (RHO2 * R2) * (R32 / R12 * E - K);
    } else {
        BRHO = 3.141592654 * RL / R2 * (RL - RHO) / R12 * Z / (R32 - RHO2);
    }

    *BX = BRHO * X;
    *BY = BRHO * Y;
    *BZ = (K - E * (R32 - 2.0 * RL * RL) / R12) / R2;
}

void crosslp(double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double RL, double AL) {
    double CAL = cos(AL);
    double SAL = sin(AL);

    double Y1 = Y * CAL - Z * SAL;
    double Z1 = Y * SAL + Z * CAL;
    double Y2 = Y * CAL + Z * SAL;
    double Z2 = -Y * SAL + Z * CAL;

    double BX1, BY1, BZ1, BX2, BY2, BZ2;
    circle(X - XC, Y1, Z1, RL, &BX1, &BY1, &BZ1);
    circle(X - XC, Y2, Z2, RL, &BX2, &BY2, &BZ2);

    *BX = BX1 + BX2;
    *BY = (BY1 + BY2) * CAL + (BZ1 - BZ2) * SAL;
    *BZ = -(BY1 - BY2) * SAL + (BZ1 + BZ2) * CAL;
}

void dipxyz(double X, double Y, double Z, double *BXX, double *BYX, double *BZX, double *BXY, double *BYY, double *BZY, double *BXZ, double *BYZ, double *BZZ) {
    double X2 = X * X;
    double Y2 = Y * Y;
    double Z2 = Z * Z;
    double R2 = X2 + Y2 + Z2;

    double XMR5 = 30574.0 / (R2 * R2 * sqrt(R2));
    double XMR53 = 3.0 * XMR5;

    *BXX = XMR5 * (3.0 * X2 - R2);
    *BYX = XMR53 * X * Y;
    *BZX = XMR53 * X * Z;

    *BXY = *BYX;
    *BYY = XMR5 * (3.0 * Y2 - R2);
    *BZY = XMR53 * Y * Z;

    *BXZ = *BZX;
    *BYZ = *BZY;
    *BZZ = XMR5 * (3.0 * Z2 - R2);
}

void condip1(double X, double Y, double Z, double D[3][42]) {
    static const double XX[14] = { -10.0, -7.0, -4.0, -4.0, 0.0, 4.0, 4.0, 7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    static const double YY[14] = { 3.0, 6.0, 3.0, 9.0, 6.0, 3.0, 9.0, 6.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    static const double ZZ[14] = { 20.0, 20.0, 4.0, 20.0, 4.0, 4.0, 20.0, 20.0, 20.0, 2.0, 3.0, 4.50, 7.0, 10.0 };

    double XSM = X + 0.16;
    double ZSM = Z;
    double RO2 = XSM * XSM + Y * Y;
    double RO = sqrt(RO2);

    double CF[5];
    double SF[5];

    CF[0] = XSM / RO;
    SF[0] = Y / RO;

    CF[1] = CF[0] * CF[0] - SF[0] * SF[0];
    SF[1] = 2.0 * SF[0] * CF[0];
    CF[2] = CF[1] * CF[0] - SF[1] * SF[0];
    SF[2] = SF[1] * CF[0] + CF[1] * SF[0];
    CF[3] = CF[2] * CF[0] - SF[2] * SF[0];
    SF[3] = SF[2] * CF[0] + CF[2] * SF[0];
    CF[4] = CF[3] * CF[0] - SF[3] * SF[0];
    SF[4] = SF[3] * CF[0] + CF[3] * SF[0];

    double R2 = RO2 + ZSM * ZSM;
    double R = sqrt(R2);
    double C = ZSM / R;
    double S = RO / R;
    double CH = sqrt(0.5 * (1.0 + C));
    double SH = sqrt(0.5 * (1.0 - C));
    double TNH = SH / CH;
    double CNH = 1.0 / TNH;

    for (int M = 0; M < 5; M++) {
        double M_NEXT = M + 1;
        double BT = M_NEXT * CF[M] / (R * S) * (pow(TNH, M_NEXT) + pow(CNH, M_NEXT));
        double BF = -0.5 * M_NEXT * SF[M] / R * (pow(TNH, M) / (CH * CH) - pow(CNH, M) / (SH * SH));
        D[0][M] = BT * C * CF[0] - BF * SF[0];
        D[1][M] = BT * C * SF[0] + BF * CF[0];
        D[2][M] = -BT * S;
    }

    XSM = X;
    ZSM = Z;

    for (int I = 0; I < 9; I++) {
        double XD, YD;

        if (I == 2 || I == 4 || I == 5) {
            XD = XX[I] * 0.08;
            YD = YY[I] * 0.08;
        } else {
            XD = XX[I] * 0.4;
            YD = YY[I] * 0.4;
        }

        double ZD = ZZ[I];

        double BX1X, BY1X, BZ1X, BX1Y, BY1Y, BZ1Y, BX1Z, BY1Z, BZ1Z;
        dipxyz(XSM - XD, Y - YD, ZSM - ZD, &BX1X, &BY1X, &BZ1X, &BX1Y, &BY1Y, &BZ1Y, &BX1Z, &BY1Z, &BZ1Z);
        double BX2X, BY2X, BZ2X, BX2Y, BY2Y, BZ2Y, BX2Z, BY2Z, BZ2Z;
        dipxyz(XSM - XD, Y + YD, ZSM - ZD, &BX2X, &BY2X, &BZ2X, &BX2Y, &BY2Y, &BZ2Y, &BX2Z, &BY2Z, &BZ2Z);
        double BX3X, BY3X, BZ3X, BX3Y, BY3Y, BZ3Y, BX3Z, BY3Z, BZ3Z;
        dipxyz(XSM - XD, Y - YD, ZSM + ZD, &BX3X, &BY3X, &BZ3X, &BX3Y, &BY3Y, &BZ3Y, &BX3Z, &BY3Z, &BZ3Z);
        double BX4X, BY4X, BZ4X, BX4Y, BY4Y, BZ4Y, BX4Z, BY4Z, BZ4Z;
        dipxyz(XSM - XD, Y + YD, ZSM + ZD, &BX4X, &BY4X, &BZ4X, &BX4Y, &BY4Y, &BZ4Y, &BX4Z, &BY4Z, &BZ4Z);

        int IX = (I + 1) * 3 + 2;
        int IY = IX + 1;
        int IZ = IY + 1;

        D[0][IX] = BX1X + BX2X - BX3X - BX4X;
        D[1][IX] = BY1X + BY2X - BY3X - BY4X;
        D[2][IX] = BZ1X + BZ2X - BZ3X - BZ4X;
        D[0][IY] = BX1Y - BX2Y - BX3Y + BX4Y;
        D[1][IY] = BY1Y - BY2Y - BY3Y + BY4Y;
        D[2][IY] = BZ1Y - BZ2Y - BZ3Y + BZ4Y;
        D[0][IZ] = BX1Z + BX2Z + BX3Z + BX4Z;
        D[1][IZ] = BY1Z + BY2Z + BY3Z + BY4Z;
        D[2][IZ] = BZ1Z + BZ2Z + BZ3Z + BZ4Z;
    }

    for (int I = 1; I <= 5; I++) {
        double ZD = ZZ[I + 8];
        double BX1X, BY1X, BZ1X, BX1Y, BY1Y, BZ1Y, BX1Z, BY1Z, BZ1Z;
        dipxyz(XSM, Y, ZSM - ZD, &BX1X, &BY1X, &BZ1X, &BX1Y, &BY1Y, &BZ1Y, &BX1Z, &BY1Z, &BZ1Z);
        double BX2X, BY2X, BZ2X, BX2Y, BY2Y, BZ2Y, BX2Z, BY2Z, BZ2Z;
        dipxyz(XSM, Y, ZSM + ZD, &BX2X, &BY2X, &BZ2X, &BX2Y, &BY2Y, &BZ2Y, &BX2Z, &BY2Z, &BZ2Z);

        int IX = 30 + I * 2;
        int IZ = IX + 1;
        D[0][IX] = BX1X - BX2X;
        D[1][IX] = BY1X - BY2X;
        D[2][IX] = BZ1X - BZ2X;
        D[0][IZ] = BX1Z + BX2Z;
        D[1][IZ] = BY1Z + BY2Z;
        D[2][IZ] = BZ1Z + BZ2Z;
    }
}

void birk1shld(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    double BXX = 0.0;
    double BYY = 0.0;
    double BZZ = 0.0;

    double RPI = 1.0 / 5.303648988;
    double RPI2 = RPI * RPI;
    double RPIY = RPI * Y;
    double SYPI = sin(RPIY);
    double CYPI = cos(RPIY);

    const double RRK_68 = 1.0 / 1.645049286;
    const double ZRRK_68 = Z * RRK_68;
    const double SZRK_68 = sin(ZRRK_68);
    const double CZRK_68 = cos(ZRRK_68);
    const double RRK2_68 = RRK_68 * RRK_68;
    double SQPR = sqrt(RPI2 + RRK2_68);
    double EPR = exp(X * SQPR);
    double AA = 1.174198045 + (double)-1.463820502;
    BXX += (-SQPR * EPR * CYPI * SZRK_68) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_68) * AA;
    BZZ += (-RRK_68 * EPR * CYPI * CZRK_68) * AA;

    const double RRK_69 = 1.0 / 3.825838190;
    const double ZRRK_69 = Z * RRK_69;
    const double SZRK_69 = sin(ZRRK_69);
    const double CZRK_69 = cos(ZRRK_69);
    const double RRK2_69 = RRK_69 * RRK_69;
    SQPR = sqrt(RPI2 + RRK2_69);
    EPR = exp(X * SQPR);
    AA = 4.840161537 + (double)-3.674506864;
    BXX += (-SQPR * EPR * CYPI * SZRK_69) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_69) * AA;
    BZZ += (-RRK_69 * EPR * CYPI * CZRK_69) * AA;

    const double RRK_70 = 1.0 / 11.66675599;
    const double ZRRK_70 = Z * RRK_70;
    const double SZRK_70 = sin(ZRRK_70);
    const double CZRK_70 = cos(ZRRK_70);
    const double RRK2_70 = RRK_70 * RRK_70;
    SQPR = sqrt(RPI2 + RRK2_70);
    EPR = exp(X * SQPR);
    AA = 82.18368896 + (double)-94.94071588;
    BXX += (-SQPR * EPR * CYPI * SZRK_70) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_70) * AA;
    BZZ += (-RRK_70 * EPR * CYPI * CZRK_70) * AA;
    
    const double RRK_71 = 1.0 / 558.9781177;
    const double ZRRK_71 = Z * RRK_71;
    const double SZRK_71 = sin(ZRRK_71);
    const double CZRK_71 = cos(ZRRK_71);
    const double RRK2_71 = RRK_71 * RRK_71;
    SQPR = sqrt(RPI2 + RRK2_71);
    EPR = exp(X * SQPR);
    AA = -4122.331796 + (double)4670.278676;
    BXX += (-SQPR * EPR * CYPI * SZRK_71) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_71) * AA;
    BZZ += (-RRK_71 * EPR * CYPI * CZRK_71) * AA;

    RPI = 1.0 / 10.40368955;
    RPI2 = RPI * RPI;
    RPIY = RPI * Y;
    SYPI = sin(RPIY);
    CYPI = cos(RPIY);

    SQPR = sqrt(RPI2 + RRK2_68);
    EPR = exp(X * SQPR);
    AA = -21.54975037 + (double)26.72661293;
    BXX += (-SQPR * EPR * CYPI * SZRK_68) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_68) * AA;
    BZZ += (-RRK_68 * EPR * CYPI * CZRK_68) * AA;
    
    SQPR = sqrt(RPI2 + RRK2_69);
    EPR = exp(X * SQPR);
    AA = -72.81365728 + (double)44.09887902;
    BXX += (-SQPR * EPR * CYPI * SZRK_69) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_69) * AA;
    BZZ += (-RRK_69 * EPR * CYPI * CZRK_69) * AA;

    SQPR = sqrt(RPI2 + RRK2_70);
    EPR = exp(X * SQPR);
    AA = 40.08073706 + (double)-51.23563510;
    BXX += (-SQPR * EPR * CYPI * SZRK_70) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_70) * AA;
    BZZ += (-RRK_70 * EPR * CYPI * CZRK_70) * AA;

    SQPR = sqrt(RPI2 + RRK2_71);
    EPR = exp(X * SQPR);
    AA = 1955.348537 + (double)-1940.971550;
    BXX += (-SQPR * EPR * CYPI * SZRK_71) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_71) * AA;
    BZZ += (-RRK_71 * EPR * CYPI * CZRK_71) * AA;
    
    RPI = 1.0 / 69.65230348;
    RPI2 = RPI * RPI;
    RPIY = RPI * Y;
    SYPI = sin(RPIY);
    CYPI = cos(RPIY);

    SQPR = sqrt(RPI2 + RRK2_68);
    EPR = exp(X * SQPR);
    AA = 794.0496433 + (double)-982.2441344;
    BXX += (-SQPR * EPR * CYPI * SZRK_68) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_68) * AA;
    BZZ += (-RRK_68 * EPR * CYPI * CZRK_68) * AA;

    SQPR = sqrt(RPI2 + RRK2_69);
    EPR = exp(X * SQPR);
    AA = 1889.837171 + (double)-558.9779727;
    BXX += (-SQPR * EPR * CYPI * SZRK_69) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_69) * AA;
    BZZ += (-RRK_69 * EPR * CYPI * CZRK_69) * AA;

    SQPR = sqrt(RPI2 + RRK2_70);
    EPR = exp(X * SQPR);
    AA = -1260.543238 + (double)1260.063802;
    BXX += (-SQPR * EPR * CYPI * SZRK_70) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_70) * AA;
    BZZ += (-RRK_70 * EPR * CYPI * CZRK_70) * AA;
    
    SQPR = sqrt(RPI2 + RRK2_71);
    EPR = exp(X * SQPR);
    AA = -293.5942373 + (double)344.7250789;
    BXX += (-SQPR * EPR * CYPI * SZRK_71) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_71) * AA;
    BZZ += (-RRK_71 * EPR * CYPI * CZRK_71) * AA;
    
    RPI = 1.0 / 466.5099509;
    RPI2 = RPI * RPI;
    RPIY = RPI * Y;
    SYPI = sin(RPIY);
    CYPI = cos(RPIY);

    SQPR = sqrt(RPI2 + RRK2_68);
    EPR = exp(X * SQPR);
    AA = -773.7002492 + (double)957.0094135;
    BXX += (-SQPR * EPR * CYPI * SZRK_68) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_68) * AA;
    BZZ += (-RRK_68 * EPR * CYPI * CZRK_68) * AA;
    
    SQPR = sqrt(RPI2 + RRK2_69);
    EPR = exp(X * SQPR);
    AA = -1824.143669 + (double)520.7994379;
    BXX += (-SQPR * EPR * CYPI * SZRK_69) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_69) * AA;
    BZZ += (-RRK_69 * EPR * CYPI * CZRK_69) * AA;
    
    SQPR = sqrt(RPI2 + RRK2_70);
    EPR = exp(X * SQPR);
    AA = 1192.484774 + (double)-1192.184565;
    BXX += (-SQPR * EPR * CYPI * SZRK_70) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_70) * AA;
    BZZ += (-RRK_70 * EPR * CYPI * CZRK_70) * AA;
    
    SQPR = sqrt(RPI2 + RRK2_71);
    EPR = exp(X * SQPR);
    AA = 89.15537624 + (double)-98.52042999;
    BXX += (-SQPR * EPR * CYPI * SZRK_71) * AA;
    BYY += (RPI * EPR * SYPI * SZRK_71) * AA;
    BZZ += (-RRK_71 * EPR * CYPI * CZRK_71) * AA;

    *BX = BXX;
    *BY = BYY;
    *BZ = BZZ;
}

void birk2shl(double X, double Y, double Z, double *HX, double *HY, double *HZ) {
    double HXX = 0.0, HYY = 0.0, HZZ = 0.0;
    double AI = 13.85650567;
    double Y_AI = Y / AI;
    double CYPI = cos(Y_AI);

    double AK = 10.21914434;
    double Z_AK = Z / AK;
    double SZRK = sin(Z_AK);
    double CZRK = cos(Z_AK);
    double SYPI = sin(Y_AI);
    double AIAI = 1.0 / (AI * AI);
    double SQPR = sqrt(AIAI + 1.0 / (AK * AK));
    double EPR = exp(X * SQPR);

    double ASUM = -111.6371348 + (double)124.5402702;
    HXX += (-SQPR * EPR * CYPI * SZRK) * ASUM;
    HYY += (EPR / AI * SYPI * SZRK) * ASUM;
    HZZ += (-EPR / AK * CYPI * CZRK) * ASUM;

    AK = 10.09021632;
    Z_AK = Z / AK;
    SZRK = sin(Z_AK);
    CZRK = cos(Z_AK);
    SYPI = sin(Y_AI);
    SQPR = sqrt(AIAI + 1.0 / (AK * AK));
    EPR = exp(X * SQPR);
    
    ASUM = 110.3735178 + (double)-122.0095905;
    HXX += (-SQPR * EPR * CYPI * SZRK) * ASUM;
    HYY += (EPR / AI * SYPI * SZRK) * ASUM;
    HZZ += (-EPR / AK * CYPI * CZRK) * ASUM;

    AI = 14.90554500;
    Y_AI = Y / AI;
    CYPI = cos(Y_AI);

    AK = 10.21914434;
    Z_AK = Z / AK;
    SZRK = sin(Z_AK);
    CZRK = cos(Z_AK);
    SYPI = sin(Y_AI);
    AIAI = 1.0 / (AI * AI);
    SQPR = sqrt(AIAI + 1.0 / (AK * AK));
    EPR = exp(X * SQPR);

    ASUM = 111.9448247 + (double)-129.1957743;
    HXX += (-SQPR * EPR * CYPI * SZRK) * ASUM;
    HYY += (EPR / AI * SYPI * SZRK) * ASUM;
    HZZ += (-EPR / AK * CYPI * CZRK) * ASUM;

    AK = 10.09021632;
    Z_AK = Z / AK;
    SZRK = sin(Z_AK);
    CZRK = cos(Z_AK);
    SYPI = sin(Y_AI);
    SQPR = sqrt(AIAI + 1.0 / (AK * AK));
    EPR = exp(X * SQPR);

    ASUM = -110.7586562 + (double)126.5649012;
    HXX += (-SQPR * EPR * CYPI * SZRK) * ASUM;
    HYY += (EPR / AI * SYPI * SZRK) * ASUM;
    HZZ += (-EPR / AK * CYPI * CZRK) * ASUM;

    *HX = HXX;
    *HY = HYY;
    *HZ = HZZ;
}

void r2_birk(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    double DELARG = 0.03;
    double DELARG1 = 0.015;
    double XKS = xksi(X, Y, Z);
    double BXSM = 0.0, BYY = 0.0, BZSM = 0.0;

    if (XKS < -(DELARG + DELARG1)) {
        r2outer(X, Y, Z, &BXSM, &BYY, &BZSM);
        BXSM *= -0.02;
        BYY *= -0.02;
        BZSM *= -0.02;
    }

    if (XKS >= -(DELARG + DELARG1) && XKS < -DELARG + DELARG1) {
        double BXSM1, BYY1, BZSM1;
        r2outer(X, Y, Z, &BXSM1, &BYY1, &BZSM1);
        double BXSM2, BYY2, BZSM2;
        r2sheet(X, Y, Z, &BXSM2, &BYY2, &BZSM2);

        double F2 = -0.02 * tksi(XKS, -DELARG, DELARG1);
        double F1 = -0.02 - F2;
        BXSM = BXSM1 * F1 + BXSM2 * F2;
        BYY = BYY1 * F1 + BYY2 * F2;
        BZSM = BZSM1 * F1 + BZSM2 * F2;
    }

    if (XKS >= -DELARG + DELARG1 && XKS < DELARG - DELARG1) {
        r2sheet(X, Y, Z, &BXSM, &BYY, &BZSM);
        BXSM *= -0.02;
        BYY *= -0.02;
        BZSM *= - 0.02;
    }

    if (XKS >= DELARG - DELARG1 && XKS < DELARG + DELARG1) {
        double BXSM1, BYY1, BZSM1;
        r2inner(X, Y, Z, &BXSM1, &BYY1, &BZSM1);
        double BXSM2, BYY2, BZSM2;
        r2sheet(X, Y, Z, &BXSM2, &BYY2, &BZSM2);

        double F1 = -0.02 * tksi(XKS, DELARG, DELARG1);
        double F2 = -0.02 - F1;
        BXSM = BXSM1 * F1 + BXSM2 * F2;
        BYY = BYY1 * F1 + BYY2 * F2;
        BZSM = BZSM1 * F1 + BZSM2 * F2;
    }

    if (XKS >= DELARG + DELARG1) {
        r2inner(X, Y, Z, &BXSM, &BYY, &BZSM);
        BXSM = -BXSM * 0.02;
        BYY = -BYY * 0.02;
        BZSM = -BZSM * 0.02;
    }

    *BX = BXSM;
    *BY = BYY;
    *BZ = BZSM;
}

void r2inner(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    double CBX[5];
    double CBY[5];
    double CBZ[5];
    
    bconic(X, Y, Z, CBX, CBY, CBZ);
    double DBX8, DBY8, DBZ8;
    loops4(X, Y, Z, &DBX8, &DBY8, &DBZ8, -8.1902, 6.5239, 5.504, 7.7815, .8573, 3.0986);
    double DBX6, DBY6, DBZ6;
    dipdistr(X - .0774, Y, Z, &DBX6, &DBY6, &DBZ6, 0);
    double DBX7, DBY7, DBZ7;
    dipdistr(X + .038, Y, Z, &DBX7, &DBY7, &DBZ7, 1);

    *BX = 154.185 * CBX[0] + -2.12446 * CBX[1] + 0.601735E-01 * CBX[2] + -.153954E-02 * CBX[3] + .355077E-04 * CBX[4] + 29.9996 * DBX6 + 262.886 * DBX7 + 99.9132 * DBX8;
    *BY = 154.185 * CBY[0] + -2.12446 * CBY[1] + 0.601735E-01 * CBY[2] + -.153954E-02 * CBY[3] + .355077E-04 * CBY[4] + 29.9996 * DBY6 + 262.886 * DBY7 + 99.9132 * DBY8;
    *BZ = 154.185 * CBZ[0] + -2.12446 * CBZ[1] + 0.601735E-01 * CBZ[2] + -.153954E-02 * CBZ[3] + .355077E-04 * CBZ[4] + 29.9996 * DBZ6 + 262.886 * DBZ7 + 99.9132 * DBZ8;
}

void bconic(double X, double Y, double Z, double *CBX, double *CBY, double *CBZ) {
    const double RO2 = X * X + Y * Y;
    const double RO = sqrt(RO2);

    const double CF = X / RO;
    const double SF = Y / RO;

    const double R2 = RO2 + Z * Z;
    const double R = sqrt(R2);
    const double C = Z / R;
    const double S = RO / R;
    const double CH = sqrt(0.5 * (1.0 + C));
    const double SH = sqrt(0.5 * (1.0 - C));
    double TNHM1 = 1.0;
    double CNHM1 = 1.0;
    const double TNH = SH / CH;
    const double CNH = 1.0 / TNH;
    double CFM1 = 1.0;
    double SFM1 = 0.0;

    for (int M = 0; M < 5; M++) {
        const double CFM = CFM1 * CF - SFM1 * SF;
        const double SFM = CFM1 * SF + SFM1 * CF;
        CFM1 = CFM;
        SFM1 = SFM;
        const double TNHM = TNHM1 * TNH;
        const double CNHM = CNHM1 * CNH;
        const double BT = (M + 1) * CFM / (R * S) * (TNHM + CNHM);
        const double BF = -0.5 * (M + 1) * SFM / R * (TNHM1 / (CH * CH) - CNHM1 / (SH * SH));
        TNHM1 = TNHM;
        CNHM1 = CNHM;
        CBX[M] = BT * C * CF - BF * SF;
        CBY[M] = BT * C * SF + BF * CF;
        CBZ[M] = -BT * S;
    }
}

void dipdistr(double X, double Y, double Z, double *BX, double *BY, double *BZ, int MODE) {
    const double X2 = X * X;
    const double RHO2 = X2 + Y * Y;

    if (MODE == 0) {
        const double R2 = RHO2 + Z * Z;
        const double R3 = R2 * sqrt(R2);
        *BX = Z / (RHO2 * RHO2) * (R2 * (Y * Y - X2) - RHO2 * X2) / R3;
        *BY = -X * Y * Z / (RHO2 * RHO2) * (2.0 * R2 + RHO2) / R3;
        *BZ = X / R3;
        return;
    }

    *BX = Z / (RHO2 * RHO2) * (Y * Y - X2);
    *BY = -2.0 * X * Y * Z / (RHO2 * RHO2);
    *BZ = X / RHO2;
}

void r2outer(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    double DBX1, DBY1, DBZ1;
    crosslp(X, Y, Z, &DBX1, &DBY1, &DBZ1, 0.55, 0.694, 0.0031);
    double DBX2, DBY2, DBZ2;
    crosslp(X, Y, Z, &DBX2, &DBY2, &DBZ2, 1.55, 2.8, 0.1375);
    double DBX3, DBY3, DBZ3;
    crosslp(X, Y, Z, &DBX3, &DBY3, &DBZ3, -0.7, 0.2, 0.9625);

    double DBX4, DBY4, DBZ4;
    circle(X + 2.994, Y, Z, 2.925, &DBX4, &DBY4, &DBZ4);

    double DBX5, DBY5, DBZ5;
    loops4(X, Y, Z, &DBX5, &DBY5, &DBZ5, -1.775, 4.3, -0.275, 2.7, 0.4312, 1.55);

    *BX = -34.105 * DBX1 + -2.00019 * DBX2 + 628.639 * DBX3 + 73.4847 * DBX4 + 12.5162 * DBX5;
    *BY = -34.105 * DBY1 + -2.00019 * DBY2 + 628.639 * DBY3 + 73.4847 * DBY4 + 12.5162 * DBY5;
    *BZ = -34.105 * DBZ1 + -2.00019 * DBZ2 + 628.639 * DBZ3 + 73.4847 * DBZ4 + 12.5162 * DBZ5;
}

void loops4(double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double YC, double ZC, double R, double THETA, double PHI) {
    const double ST = sin(THETA);
    const double CT = cos(THETA);
    const double SP = sin(PHI);
    const double CP = cos(PHI);

    double XS = (X - XC) * CP + (Y - YC) * SP;
    double YSS = (Y - YC) * CP - (X - XC) * SP;
    double ZS = Z - ZC;
    double XSS = XS * CT - ZS * ST;
    double ZSS = ZS * CT + XS * ST;

    double BXSS, BYS, BZSS;
    circle(XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS);
    double BXS = BXSS * CT + BZSS * ST;
    const double BZ1 = BZSS * CT - BXSS * ST;
    const double BX1 = BXS * CP - BYS * SP;
    const double BY1 = BXS * SP + BYS * CP;

    XS = (X - XC) * CP - (Y + YC) * SP;
    YSS = (Y + YC) * CP + (X - XC) * SP;
    ZS = Z - ZC;
    XSS = XS * CT - ZS * ST;
    ZSS = ZS * CT + XS * ST;

    circle(XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS);
    BXS = BXSS * CT + BZSS * ST;
    const double BZ2 = BZSS * CT - BXSS * ST;
    const double BX2 = BXS * CP + BYS * SP;
    const double BY2 = -BXS * SP + BYS * CP;
    
    XS = -(X - XC) * CP + (Y + YC) * SP;
    YSS = -(Y + YC) * CP - (X - XC) * SP;
    ZS = Z + ZC;
    XSS = XS * CT - ZS * ST;
    ZSS = ZS * CT + XS * ST;

    circle(XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS);
    BXS=BXSS*CT+BZSS*ST;
    const double BZ3 = BZSS * CT - BXSS * ST;
    const double BX3 = -BXS * CP - BYS * SP;
    const double BY3 = BXS * SP - BYS * CP;

    XS = -(X - XC) * CP - (Y - YC) * SP;
    YSS = -(Y - YC) * CP + (X - XC) * SP;
    ZS = Z + ZC;
    XSS = XS * CT - ZS * ST;
    ZSS = ZS * CT + XS * ST;

    circle(XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS);
    BXS = BXSS * CT + BZSS * ST;

    *BX = BX1 + BX2 + BX3 + (-BXS * CP + BYS * SP);
    *BY = BY1 + BY2 + BY3 + (-BXS * SP - BYS * CP);
    *BZ = BZ1 + BZ2 + BZ3 + (BZSS * CT - BXSS * ST);
}

void r2sheet(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double XKS = xksi(X, Y, Z);
    const double XKS2 = XKS * XKS;
    const double T1X = XKS / sqrt(XKS2 + 0.483750e-01 * 0.483750e-01);
    const double PNONX72 = 0.396953e-01 * 0.396953e-01;
    const double XKS2_PNONX7 = sqrt(XKS2 + PNONX72);
    const double T2X = (PNONX72 * 0.396953e-01) / (XKS2_PNONX7 * XKS2_PNONX7 * XKS2_PNONX7);
    const double T3X = XKS / pow(sqrt(XKS2 + 0.579023e-01 * 0.579023e-01), 5) * 3.493856 * pow(0.579023e-01, 4);
    const double T1Y = XKS / sqrt(XKS2 + 0.478750e-01 * 0.478750e-01);
    const double PNONY72 = 0.363750e-01 * 0.363750e-01;
    const double XKS2_PNONY7 = sqrt(XKS2 + PNONY72);
    const double T2Y = (PNONY72 * 0.363750e-01) / (XKS2_PNONY7 * XKS2_PNONY7 * XKS2_PNONY7);
    const double T3Y = XKS / pow(sqrt(XKS2 + 0.567500e-01 * 0.567500e-01), 5) * 3.493856 * pow(0.567500e-01, 4);
    const double T1Z = XKS / sqrt(XKS2 + 0.355625e-01 * 0.355625e-01);
    const double PNONZ72 = 0.318750e-01 * 0.318750e-01;
    const double XKS2_PNON27 = sqrt(XKS2 + PNONZ72);
    const double T2Z = (PNONZ72 * 0.318750e-01) / (XKS2_PNON27 * XKS2_PNON27 * XKS2_PNON27);
    const double T3Z = XKS / pow(sqrt(XKS2 + 0.538750e-01 * 0.538750e-01), 5) * 3.493856 * pow(0.538750e-01, 4);
    const double RHO2 = X * X + Y * Y;
    const double R = sqrt(RHO2 + Z * Z);
    const double RHO = sqrt(RHO2);
    const double C1P = X / RHO;
    const double S1P = Y / RHO;
    const double S2P = 2.0 * S1P * C1P;
    const double C2P = C1P * C1P - S1P * S1P;
    const double S3P = S2P * C1P + C2P * S1P;
    const double C3P = C2P * C1P - S2P * S1P;
    const double S4P = S3P * C1P + C3P * S1P;
    const double CT = Z / R;
    const double CT2 = CT * CT;
    const double CT2_M1 = CT2 - 1.0;
    const double S_TMP = -2. * 2.718281828459;

    double S1 = sqrt(S_TMP * -19.0969) * CT * exp(-19.0969 * CT2);
    double S2 = sqrt(S_TMP * -9.28828) * CT * exp(-9.28828 * CT2);
    double S3 = sqrt(S_TMP * -0.129687) * CT * exp(-0.129687 * CT2);
    double S4 = CT * exp(5.58594 * CT2_M1);
    double S5 = CT * exp(22.5055 * CT2_M1);

    *BX = S1 * ((8.07190 + -7.39582 * T1X + -7.62341 * T2X + 0.684671 * T3X) 
    + C1P * (-13.5672 + 11.6681 * T1X + 13.1154 * T2X + -0.890217 * T3X) 
    + C2P * (7.78726 + -5.38346 * T1X + -8.08738 * T2X + 0.609385 * T3X) 
    + C3P * (-2.70410 + 3.53741 * T1X + 3.15549 * T2X + -1.11069 * T3X))
    + S2 * ((-8.47555 + 0.278122 * T1X + 2.73514 * T2X + 4.55625 * T3X) 
    + C1P * (13.1134 + 1.15848 * T1X + -3.52648 * T2X + -8.24698 * T3X) 
    + C2P * (-6.85710 + -2.81369 * T1X + 2.03795 * T2X + 4.64383 * T3X) 
    + C3P * (2.49309 + -1.22041 * T1X + -1.67432 * T2X + -0.422526 * T3X))
    + S3 * ((-5.39796 + 7.10326 * T1X + 5.53730 * T2X + -13.1918 * T3X) 
    + C1P * (4.67853 + -7.60329 * T1X + -2.53066 * T2X + 7.76338 * T3X) 
    + C2P * (5.60165 + 5.34816 * T1X + -4.56441 * T2X + 7.05976 * T3X) 
    + C3P * (-2.62723 + -0.529078 * T1X + 1.42019 * T2X + -2.93919 * T3X))
    + S4 * ((55.6338 + -1.55181 * T1X + 39.8311 * T2X + -80.6561 * T3X) 
    + C1P * (-46.9655 + 32.8925 * T1X + -6.32296 * T2X + 19.7841 * T3X) 
    + C2P * (124.731 + 10.4347 * T1X + -30.7581 * T2X + 102.680 * T3X) 
    + C3P * (-47.4037 + -3.31278 * T1X + 9.37141 * T2X + -50.0268 * T3X))
    + S5 * ((-533.319 + 110.426 * T1X + 1000.20 * T2X + -1051.40 * T3X) 
    + C1P * (1619.48 + 589.855 * T1X + -1462.73 * T2X + 1087.10 * T3X) 
    + C2P * (-1994.73 + -1654.12 * T1X + 1263.33 * T2X + -260.210 * T3X) 
    + C3P * (1424.84 + 1255.71 * T1X + -956.733 * T2X + 219.946 * T3X));

    S1 = sqrt(S_TMP * -13.6750) * CT * exp(-13.6750 * CT2);
    S2 = sqrt(S_TMP * -6.70625) * CT * exp(-6.70625 * CT2);
    S3 = CT * exp(2.31875 * CT2_M1);
    S4 = CT * exp(11.4062 * CT2_M1);
    S5 = CT * exp(20.4562 * CT2_M1);

    *BY = S1 * (S1P * (-9.08427 + 10.6777 * T1Y + 10.3288 * T2Y + -0.969987 * T3Y)
    + S2P * (6.45257 + -8.42508 * T1Y + -7.97464 * T2Y + 1.41996 * T3Y)
    + S3P * (-1.92490 + 3.93575 * T1Y + 2.83283 * T2Y + -1.48621 * T3Y)
    + S4P * (0.244033 + -0.757941 * T1Y + -0.386557 * T2Y + 0.344566 * T3Y))
    + S2 * (S1P * (9.56674 + -2.5365 * T1Y + -3.32916 * T2Y + -5.86712 * T3Y)
    + S2P * (-6.19625 + 1.83879 * T1Y + 2.52772 * T2Y + 4.34417 * T3Y)
    + S3P * (1.87268 + -2.13213 * T1Y + -1.69134 * T2Y + -0.176379 * T3Y)
    + S4P * (-0.261359 + 0.566419 * T1Y + 0.3138 * T2Y + -0.134699 * T3Y))
    + S3 * (S1P * (-3.83086 + -8.4154 * T1Y + 4.77005 * T2Y + -9.31479 * T3Y)
    + S2P * (37.5715 + 19.3992 * T1Y + -17.9582 * T2Y + 36.4604 * T3Y)
    + S3P * (-14.9993 + -3.1442 * T1Y + 6.17409 * T2Y + -15.5519 * T3Y)
    + S4P * (2.28621 + -0.00891549 * T1Y + -0.462912 * T2Y + 2.47314 * T3Y))
    + S4 * (S1P * (41.7555 + 208.614 * T1Y + -45.7861 * T2Y + -77.8687 * T3Y)
    + S2P * (239.357 + -67.9226 * T1Y + 66.8743 * T2Y + 238.534 * T3Y)
    + S3P * (-112.136 + 16.2069 * T1Y + -40.4706 * T2Y + -134.328 * T3Y)
    + S4P * (21.56 + -0.201725 * T1Y + 2.21 * T2Y + 32.5855 * T3Y))
    + S5 * (S1P * (-108.217 + -1005.98 * T1Y + 585.753 * T2Y + 323.668 * T3Y)
    + S2P * (-817.056 + 235.750 * T1Y + -560.965 * T2Y + -576.892 * T3Y)
    + S3P * (684.193 + 85.0275 * T1Y + 168.394 * T2Y + 477.776 * T3Y)
    + S4P * (-289.253 + -123.216 * T1Y + 75.6501 * T2Y + -178.605 * T3Y));

    S1 = exp(-16.7125 * CT2);
    S2 = exp(-16.4625 * CT2);
    S3 = exp(-0.1625 * CT2);
    S4 = exp(5.1 * (CT2 - 1.0));
    S5 = exp(23.7125 * (CT2 - 1.0));

    *BZ = S1 * ((1167.61 + -917.782 * T1Z + -1253.2 * T2Z + -274.128 * T3Z)
    + C1P * (-1538.75 + 1257.62 * T1Z + 1745.07 * T2Z + 113.479 * T3Z)
    + C2P * (393.326 + -426.858 * T1Z + -641.1 * T2Z + 190.833 * T3Z)
    + C3P * (-29.9435 + -1.04881 * T1Z + 117.125 * T2Z + -25.7663 * T3Z))
    + S2 * ((-1168.16 + 910.247 * T1Z + 1239.31 * T2Z + 289.515 * T3Z)
    + C1P * (1540.56 + -1248.29 * T1Z + -1727.61 * T2Z + -131.785 * T3Z)
    + C2P * (-394.577 + 426.163 * T1Z + 637.422 * T2Z + -187.965 * T3Z)
    + C3P * (30.0348 + 0.221898 * T1Z + -116.68 * T2Z + 26.0291 * T3Z))
    + S3 * ((12.6804 + 4.84091 * T1Z + 1.18166 * T2Z + -2.75946 * T3Z)
    + C1P * (-17.9822 + -6.80357 * T1Z + -1.47134 * T2Z + 3.02266 * T3Z)
    + C2P * (4.79648 + 0.665255 * T1Z + -0.256229 * T2Z + -0.0857282 * T3Z)
    + C3P * (-0.588997 + 0.0634812 * T1Z + 0.164303 * T2Z + -0.15285 * T3Z))
    + S4 * ((22.2524 + -22.4376 * T1Z + -3.85595 * T2Z + 6.07625 * T3Z)
    + C1P * (-105.959 + -41.6698 * T1Z + 0.378615 * T2Z + 1.55958 * T3Z)
    + C2P * (44.3981 + 18.8521 * T1Z + 3.19466 * T2Z + 5.89142 * T3Z)
    + C3P * (-8.63227 + -2.36418 * T1Z + -1.027 * T2Z + -2.31515 * T3Z))
    + S5 * ((1035.38 + 2040.66 * T1Z + -131.881 * T2Z + -744.533 * T3Z)
    + C1P * (-3274.93 + -4845.61 * T1Z + 482.438 * T2Z + 1567.43 * T3Z)
    + C2P * (1354.02 + 2040.47 * T1Z + -151.653 * T2Z + -845.012 * T3Z)
    + C3P * (-111.723 + -265.343 * T1Z + -26.1171 * T2Z + 216.632 * T3Z));
}

double xksi(double X, double Y, double Z) {
    const double X2 = X * X;
    const double Y2 = Y * Y;
    const double Z2 = Z * Z;
    const double R2 = X2 + Y2 + Z2;
    const double R = sqrt(R2);
    const double XR = X / R;
    const double YR = Y / R;
    const double ZR = Z / R;

    const double DR = 7.50937;
    const double DR2 = DR * DR;
    const double R0 = 1.21563;
    double PR = 0.0;
    if (R >= R0) {
        PR = sqrt((R - R0) * (R - R0) + DR2) - DR;
    }

    const double F = X + PR * (0.305662 + -0.383593 * XR + 0.2677733 * XR * XR + -0.097656 * YR * YR + -0.636034 * ZR * ZR);
    const double G = Y + PR * (-0.359862 * YR + 0.424706 * XR * YR);
    const double H = Z + PR * (-0.126366 * ZR + 0.292578 * XR * ZR);
    const double G2 = G * G;
    
    const double FGH = F * F + G2 + H * H;
    const double SQFGH32 = sqrt(FGH);
    const double FGH32 = SQFGH32 * SQFGH32 * SQFGH32;
    const double FCHSG2 = (F * F) + G2;

    if (FCHSG2 < 1.0e-5) {
        return -1.0;
    }

    const double THETAS = sin(0.3665191 + 0.5 * 0.09599309 * (1.0 - F / sqrt(FCHSG2)));
    return (FCHSG2 / FGH32) - (THETAS * THETAS);
}

double tksi(double XKSI, double XKS0, double DXKSI) {
    double TDZ3 = 2.0 * DXKSI * DXKSI * DXKSI;

    if (XKSI - XKS0 < -DXKSI) {
        return 0.0;
    }
    
    if (XKSI - XKS0 >= DXKSI) {
        return 1.0;
    }
    
    if (XKSI >= XKS0 - DXKSI && XKSI < XKS0) {
        double XKSI_XKS0_DXKSI = XKSI - XKS0 + DXKSI;
        double BR3 = XKSI_XKS0_DXKSI * XKSI_XKS0_DXKSI * XKSI_XKS0_DXKSI;
        return 1.5 * BR3 / (TDZ3 + BR3);
    }
    
    double XKSI_XKS0_DXKSI = XKSI - XKS0 - DXKSI;
    double BR3 = XKSI_XKS0_DXKSI * XKSI_XKS0_DXKSI * XKSI_XKS0_DXKSI;
    return 1.0 + 1.5 * BR3 / (TDZ3 - BR3);
}

void dipole_t96(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double P = X * X;
    const double U = Z * Z;
    const double V = 3.0 * Z * X;
    const double T = Y * Y;
    const double PTU = sqrt(P + T + U);
    const double Q = 30574.0 / (PTU * PTU * PTU * PTU * PTU);

    *BX = Q * ((T + U - 2.0 * P) - V);
    *BY = -3.0 * Y * Q * Z;
    *BZ = Q * (P + T - 2.0 * U);
}
