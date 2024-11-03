#include <math.h>
#include "ts05.h"
#include "deformed.h"
#include "birk_tot.h"
#include "full_rc.h"

void ts05_gsm(const double *PARMOD, double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    static const double A[69] = {1.00000,5.44118,0.891995,9.09684,0.00000,-7.18972,12.2700,
                          -4.89408,0.00000,0.870536,1.36081,0.00000,0.688650,0.602330,
                          0.00000,0.316346,1.22728,-0.363620E-01,-0.405821,0.452536,
                          0.755831,0.215662,0.152759,5.96235,23.2036,11.2994,69.9596,
                          0.989596,-0.132131E-01,0.985681,0.344212E-01,1.02389,0.207867,
                          1.51220,0.682715E-01,1.84714,1.76977,1.37690,0.696350,0.343280,
                          3.28846,111.293,5.82287,4.39664,0.383403,0.648176,0.318752E-01,
                          0.581168,1.15070,0.843004,0.394732,0.846509,0.916555,0.550920,
                          0.180725,0.898772,0.387365,2.26596,1.29123,0.436819,1.28211,
                          1.33199,.405553,1.6229,.699074,1.26131,2.42297,.537116,.619441};

    const double PDYN = PARMOD[0];
    const double DST = PARMOD[1]*0.8 - 13.*sqrt(PDYN);
    const double BYIMF = PARMOD[2];
    const double BZIMF = PARMOD[3];
    const double W1 = PARMOD[4];
    const double W2 = PARMOD[5];
    const double W3 = PARMOD[6];
    const double W4 = PARMOD[7];
    const double W5 = PARMOD[8];
    const double W6 = PARMOD[9];

    const double XAPPA = pow((PDYN / 2.), A[22]);
    const double XAPPA3 = XAPPA*XAPPA*XAPPA;
    const double XX = X * XAPPA;
    const double YY = Y * XAPPA;
    const double ZZ = Z * XAPPA;

    const double AM = 34.586 / XAPPA;
    const double X0 = 3.4397 / XAPPA;

    const double FACTIMF = A[19];
    const double OIMFY = BYIMF * FACTIMF;
    const double OIMFZ = BZIMF * FACTIMF;

    double XMXM = AM + X - X0;
    if (XMXM < 0) {
        XMXM = 0.;
    }

    const double AXX0 = XMXM*XMXM;
    const double ASQ = AM*AM;
    const double RHO2 = Y*Y + Z*Z;
    const double ARO_AXX0 = (ASQ + RHO2) + AXX0;
    const double SIGMA = sqrt((ARO_AXX0 + sqrt(ARO_AXX0*ARO_AXX0 - 4. * ASQ * AXX0)) / (2. * ASQ));

    double QX = 0., QY = 0., QZ = 0.;
    if(SIGMA < 1.1960 + 0.005) {
        double BXCF, BYCF, BZCF;
        double CFX, CFY, CFZ;
        shlcar3x3(XX, YY, ZZ, &CFX, &CFY, &CFZ);
        BXCF = CFX*XAPPA3;
        BYCF = CFY*XAPPA3;
        BZCF = CFZ*XAPPA3;

        double ZNAM, BXT1, BYT1, BZT1, BXT2, BYT2, BZT2;
        double DSTT = -20.;
        if(DST < DSTT) {
            DSTT = DST;
        }
        ZNAM = pow(fabs(DSTT), 0.37);
        const double DXSHIFT1 = A[23] - A[24]/ZNAM;
        const double DXSHIFT2 = A[25] - A[26]/ZNAM;
        const double D = A[35]*exp(-W1/A[36]) + A[68];
        warped(XX, YY, ZZ, DXSHIFT1, DXSHIFT2, D, &BXT1, &BYT1, &BZT1, &BXT2, &BYT2, &BZT2);

        double BXR11, BYR11, BZR11, BXR12, BYR12, BZR12, BXR21, BYR21, BZR21, BXR22, BYR22, BZR22;
        ZNAM = fabs(DST);
        if(DST >= -20.) {
            ZNAM = 20.;
        }
        const double XKAPPA1 = A[31]*pow(ZNAM/20., A[32]);
        const double XKAPPA2 = A[33]*pow(ZNAM/20., A[34]);

        birk_tot(XX, YY, ZZ, XKAPPA1, XKAPPA2,
                    &BXR11, &BYR11, &BZR11, &BXR12, &BYR12, &BZR12, &BXR21, &BYR21, &BZR21, &BXR22, &BYR22, &BZR22);

        double BXSRC, BYSRC, BZSRC, BXPRC, BYPRC, BZPRC;
        const double PHI = A[37];
        ZNAM = fabs(DST);
        if(DST >= -20.) {
            ZNAM = 20.;
        }
        const double SC_SY = A[27]*pow((20./ZNAM), A[28])*XAPPA;
        const double SC_PR = A[29]*pow((20./ZNAM), A[30])*XAPPA;
        full_rc(XX, YY, ZZ, PHI, SC_SY, SC_PR, &BXSRC, &BYSRC, &BZSRC, &BXPRC, &BYPRC, &BZPRC);

        double HXIMF, HYIMF, HZIMF;
        HXIMF = 0.;
        HYIMF = BYIMF;
        HZIMF = BZIMF;

        const double DLP1 = pow((PDYN/2.), A[20]);
        const double DLP2 = pow((PDYN/2.), A[21]);

        const double TAMP1 = A[1] + A[2]*DLP1 + A[3]*(A[38]*W1)/sqrt(W1*W1 + A[38]*A[38]);
        const double TAMP2 = A[5] + A[6]*DLP2 + A[7]*(A[39]*W2)/sqrt(W2*W2 + A[39]*A[39]);

        const double A_SRC = A[9] + A[10]*(A[40]*W3)/sqrt(W3*W3 + A[40]*A[40]) + A[11]*DST;
        const double A_PRC = A[12] + A[13]*(A[41]*W4)/sqrt(W4*W4 + A[41]*A[41]) + A[14]*DST;

        const double A_R11 = A[15] + A[16]*(A[42]*W5)/sqrt(W5*W5 + A[42]*A[42]);
        const double A_R21 = A[17] + A[18]*(A[43]*W6)/sqrt(W6*W6 + A[43]*A[43]);

        const double BBX = A[0]*BXCF + TAMP1*BXT1 + TAMP2*BXT2 + A_SRC*BXSRC + A_PRC*BXPRC + A_R11*BXR11 + A_R21*BXR21 + A[19]*HXIMF;
        const double BBY = A[0]*BYCF + TAMP1*BYT1 + TAMP2*BYT2 + A_SRC*BYSRC + A_PRC*BYPRC + A_R11*BYR11 + A_R21*BYR21 + A[19]*HYIMF;
        const double BBZ = A[0]*BZCF + TAMP1*BZT1 + TAMP2*BZT2 + A_SRC*BZSRC + A_PRC*BZPRC + A_R11*BZR11 + A_R21*BZR21 + A[19]*HZIMF;

        if(SIGMA < 1.1960-0.005) {
            *BX = BBX;
            *BY = BBY;
            *BZ = BBZ;
        } else {
            dipole(X, Y, Z, &QX, &QY, &QZ);
            const double FINT = 0.5*(1. - (SIGMA-1.1960)/0.005);
            const double FEXT = 0.5*(1. + (SIGMA-1.1960)/0.005);
            *BX = (BBX + QX)*FINT - QX;
            *BY = (BBY + QY)*FINT+OIMFY*FEXT - QY;
            *BZ = (BBZ + QZ)*FINT+OIMFZ*FEXT - QZ;
        }
    } else {
        dipole(X, Y, Z, &QX, &QY, &QZ);
        *BX = -QX;
        *BY = OIMFY - QY;
        *BZ = OIMFZ - QZ;
    }
}

void shlcar3x3(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double P1 = 9.620648151;
    const double P12 = 1./(P1*P1);
    const double R1 = 12.44199571;
    const double R12 = 1./(R1*R1);

    const double Z1_R1 = Z/R1;
    const double Z1_R1_CZR = cos(Z1_R1);
    const double Z1_R1_SZR = sin(Z1_R1);
    const double Y_P1 = Y/P1;
    double CYP = cos(Y_P1);
    double SYP = sin(Y_P1);
    double SQPR = sqrt(P12 + R12);
    double EXPR = exp(SQPR*X);
    const double HX1 = -SQPR*EXPR*CYP*Z1_R1_SZR;
    const double HY1 = EXPR/P1*SYP*Z1_R1_SZR;
    const double HZ1 = -EXPR*CYP/R1*Z1_R1_CZR;

    const double R2 = 5.122226936;
    const double R22 = 1./(R2*R2);
    SQPR = sqrt(P12 + R22);
    const double Z1_R2 = Z/R2;
    const double Z1_R2_CZR = cos(Z1_R2);
    const double Z1_R2_SZR = sin(Z1_R2);
    EXPR = exp(SQPR*X);
    const double HX2 = -SQPR*EXPR*CYP*Z1_R2_SZR;
    const double HY2 = EXPR/P1*SYP*Z1_R2_SZR;
    const double HZ2 = -EXPR*CYP/R2*Z1_R2_CZR;

    const double R3 = 6.982039615;
    const double R32 = R3*R3;
    const double R32_1 = 1./R32;
    SQPR = sqrt(P12 + R32_1);
    const double Z1_R3 = Z/R3;
    const double Z1_R3_CZR = cos(Z1_R3);
    const double Z1_R3_SZR = sin(Z1_R3);
    EXPR = exp(SQPR*X);
    const double HX3 = -EXPR*CYP*(SQPR*Z*Z1_R3_CZR + Z1_R3_SZR/R3*(X + 1./SQPR));
    const double HY3 = EXPR/P1*SYP*(Z*Z1_R3_CZR + X/R3*Z1_R3_SZR/SQPR);
    const double HZ3 = -EXPR*CYP*(Z1_R3_CZR*(1. + X/R32/SQPR) - Z1_R3*Z1_R3_SZR);

    const double P2 = 6.082014949;
    const double P22 = 1./(P2*P2);
    SQPR = sqrt(P22 + R12);
    const double Y_P2 = Y/P2;
    CYP = cos(Y_P2);
    SYP = sin(Y_P2);
    EXPR = exp(SQPR*X);
    const double HX4 = -SQPR*EXPR*CYP*Z1_R1_SZR;
    const double HY4 = EXPR/P2*SYP*Z1_R1_SZR;
    const double HZ4 = -EXPR*CYP/R1*Z1_R1_CZR;

    SQPR = sqrt(P22 + R22);
    EXPR = exp(SQPR*X);
    const double HX5 = -SQPR*EXPR*CYP*Z1_R2_SZR;
    const double HY5 = EXPR/P2*SYP*Z1_R2_SZR;
    const double HZ5 = -EXPR*CYP/R2*Z1_R2_CZR;

    SQPR = sqrt(P22 + R32_1);
    EXPR = exp(SQPR*X);
    const double HX6 = -EXPR*CYP*(SQPR*Z*Z1_R3_CZR + Z1_R3_SZR/R3*(X + 1./SQPR));
    const double HY6 = EXPR/P2*SYP*(Z*Z1_R3_CZR + X/R3*Z1_R3_SZR/SQPR);
    const double HZ6 = -EXPR*CYP*(Z1_R3_CZR*(1. + X/R32/SQPR)-Z1_R3*Z1_R3_SZR);

    const double P3 = 27.75216226;
    const double P32 = 1./(P3*P3);
    SQPR = sqrt(P32 + R12);
    const double Y_P3 = Y/P3;
    CYP = cos(Y_P3);
    SYP = sin(Y_P3);
    EXPR = exp(SQPR*X);
    const double HX7 = -SQPR*EXPR*CYP*Z1_R1_SZR;
    const double HY7 = EXPR/P3*SYP*Z1_R1_SZR;
    const double HZ7 = -EXPR*CYP/R1*Z1_R1_CZR;

    SQPR = sqrt(P32 + R22);
    EXPR = exp(SQPR*X);
    const double HX8 = -SQPR*EXPR*CYP*Z1_R2_SZR;
    const double HY8 = EXPR/P3*SYP*Z1_R2_SZR;
    const double HZ8 = -EXPR*CYP/R2*Z1_R2_CZR;

    SQPR = sqrt(P32 + R32_1);
    EXPR = exp(SQPR*X);
    const double HX9 = -EXPR*CYP*(SQPR*Z*Z1_R3_CZR + Z1_R3_SZR/R3*(X + 1./SQPR));
    const double HY9 = EXPR/P3*SYP*(Z*Z1_R3_CZR + X/R3*Z1_R3_SZR/SQPR);
    const double HZ9 = -EXPR*CYP*(Z1_R3_CZR*(1. + X/R32/SQPR) - Z1_R3*Z1_R3_SZR);

    const double A1 = -901.2327248 + (double)895.8011176;
    const double A2 = 817.6208321 - (double)845.5880889;
    const double A3 = -83.73539535 + (double)86.58542841;
    const double A4 = 336.8781402 - (double)329.3619944;
    const double A5 = -311.2947120 + (double)308.6011161;
    const double A6 = 31.94469304 - (double)31.30824526;
    const double A7 = 125.8739681 - (double)372.3384278;
    const double A8 = -235.4720434 + (double)286.7594095;
    const double A9 = 21.86305585 - (double)27.42344605;
    *BX = A1*HX1 + A2*HX2 + A3*HX3 + A4*HX4 + A5*HX5 + A6*HX6 + A7*HX7 + A8*HX8 + A9*HX9;
    *BY = A1*HY1 + A2*HY2 + A3*HY3 + A4*HY4 + A5*HY5 + A6*HY6 + A7*HY7 + A8*HY8 + A9*HY9;
    *BZ = A1*HZ1 + A2*HZ2 + A3*HZ3 + A4*HZ4 + A5*HZ5 + A6*HZ6 + A7*HZ7 + A8*HZ8 + A9*HZ9;
}

void dipole(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double P = X*X;
    const double U = Z*Z;
    const double T = Y*Y;
    const double SQRT_PTU = sqrt(P + T + U);
    const double Q = 30115. / (SQRT_PTU*SQRT_PTU*SQRT_PTU*SQRT_PTU*SQRT_PTU);
    *BX = Q * 3. * Z * X;
    *BY = Q * -3. * Z * Y;
    *BZ = Q * (P + T - 2.*U);
}
