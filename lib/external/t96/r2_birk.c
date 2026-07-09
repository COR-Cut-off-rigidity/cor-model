#include <math.h>
#include "r2_birk.h"
#include "circle.h"

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

    double R2X, R2Y, R2Z;
    birk2shl(X, Y, Z, &R2X, &R2Y, &R2Z);

    *BX = R2X + BXSM;
    *BY = R2Y + BYY;
    *BZ = R2Z + BZSM;
}