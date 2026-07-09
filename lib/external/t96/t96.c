#include <math.h>
#include "t96.h"
#include "tail_rc.h"
#include "birk_tot.h"
#include "r2_birk.h"

static void dipole(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
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

void t96_gsm(const ExternalModelContext *ctx, double X_GSM, double Y_GSM, double Z_GSM, double *BX_GSM, double *BY_GSM, double *BZ_GSM) {
    const double PDYN = ctx->pdyn;
    const double DEPR = 0.8 * ctx->dst - 13.0 * sqrt(PDYN);
    const double BYIMF = ctx->by_imf;
    const double BZIMF = ctx->bz_imf;
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
    const double YS = Y_GSM * CT - Z_GSM * ST;
    const double ZS = Z_GSM * CT + Y_GSM * ST;
    const double YS_DELIMFY = YS / 10.0;
    const double FACTIMF = exp(X_GSM / 20.0 - YS_DELIMFY * YS_DELIMFY);
    const double OIMFY = 0.7850 * BYIMF * FACTIMF;
    const double OIMFZ = 0.7850 * BZIMF * FACTIMF;
    const double XAPPA = pow(PDYN / 2.0, 0.14);
    const double XAPPA3 = XAPPA * XAPPA * XAPPA;
    const double AM = 70.0 / XAPPA;
    const double RHO2 = Y_GSM * Y_GSM + Z_GSM * Z_GSM;
    const double ASQ = AM * AM;

    double XMXM = AM + X_GSM - (5.48 / XAPPA);
    if (XMXM < 0.0) {
        XMXM = 0.0;
    }

    const double AXX0 = XMXM * XMXM;
    const double ARO = ASQ + RHO2;
    const double ARO_AXX0 = ARO + AXX0;
    const double SIGMA = sqrt((ARO + AXX0 + sqrt(ARO_AXX0 * ARO_AXX0 - 4.0 * ASQ * AXX0)) / (2.0 * ASQ));

    if (SIGMA >= 1.08 + 0.005) {
        double QX, QY, QZ;
        dipole(X_GSM, Y_GSM, Z_GSM, &QX, &QY, &QZ);

        *BX_GSM = -QX;
        *BY_GSM = OIMFY - QY;
        *BZ_GSM = OIMFZ - QZ;
        return;
    }

    const double FACTEPS = EPS / 3630.7 - 1.0;
    const double FACTPD = sqrt(PDYN / 2.0) - 1.0;
    const double RCAMPL = -1.162 * DEPR;
    const double TAMPL2 = 22.344 + 18.50 * FACTPD + 2.602 * FACTEPS;
    const double TAMPL3 = 6.903 + 5.287 * FACTPD;
    const double B1AMPL = 0.5790 + 0.4462 * FACTEPS;
    const double B2AMPL = 20.0 * B1AMPL;
    const double XX = X_GSM * XAPPA;
    const double YY = Y_GSM * XAPPA;
    const double ZZ = Z_GSM * XAPPA;

    double CFX, CFY, CFZ;
    cylharm(XX, YY, ZZ, &CFX, &CFY, &CFZ);

    double BXRC, BYRC, BZRC, BXT2, BYT2, BZT2, BXT3, BYT3, BZT3;
    tailrc(XX, YY, ZZ, &BXRC, &BYRC, &BZRC, &BXT2, &BYT2, &BZT2, &BXT3, &BYT3, &BZT3);

    double R1X, R1Y, R1Z;
    birk1tot_02(XX, YY, ZZ, &R1X, &R1Y, &R1Z);

    double R2X, R2Y, R2Z;
    r2_birk(XX, YY, ZZ, &R2X, &R2Y, &R2Z);

    double RIMFX, RIMFYS, RIMFZS;
    intercon(XX, YS * XAPPA, ZS * XAPPA, &RIMFX, &RIMFYS, &RIMFZS);

    const double RIMFY = RIMFYS * CT + RIMFZS * ST;
    const double RIMFZ = RIMFZS * CT - RIMFYS * ST;
    const double RIMFAMPL = 0.7850 * Bt;

    const double FX = CFX * XAPPA3 + RCAMPL * BXRC + TAMPL2 * BXT2 + TAMPL3 * BXT3 + B1AMPL * R1X + B2AMPL * R2X + RIMFAMPL * RIMFX;
    const double FY = CFY * XAPPA3 + RCAMPL * BYRC + TAMPL2 * BYT2 + TAMPL3 * BYT3 + B1AMPL * R1Y + B2AMPL * R2Y + RIMFAMPL * RIMFY;
    const double FZ = CFZ * XAPPA3 + RCAMPL * BZRC + TAMPL2 * BZT2 + TAMPL3 * BZT3 + B1AMPL * R1Z + B2AMPL * R2Z + RIMFAMPL * RIMFZ;

    if (SIGMA < 1.08 - 0.005) {
        *BX_GSM = FX;
        *BY_GSM = FY;
        *BZ_GSM = FZ;
        return;
    }

    const double FINT = 0.5 * (1.0 - (SIGMA - 1.08) / 0.005);
    const double FEXT = 0.5 * (1.0 + (SIGMA - 1.08) / 0.005);
    
    double QX, QY, QZ;
    dipole(X_GSM, Y_GSM, Z_GSM, &QX, &QY, &QZ);
    
    *BX_GSM = (FX + QX) * FINT - QX;
    *BY_GSM = (FY + QY) * FINT + OIMFY * FEXT - QY;
    *BZ_GSM = (FZ + QZ) * FINT + OIMFZ * FEXT - QZ;
}
