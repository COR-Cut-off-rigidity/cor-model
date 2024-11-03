#include <math.h>
#include "deformed.h"

void warped(double X, double Y, double Z,
            double DXSHIFT1, double DXSHIFT2, double D,
            double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2) {

    const double RHO2 = Y*Y + Z*Z;
    const double RHO = sqrt(RHO2);
    double PHI = 0.;
    double CPHI = 1.;
    double SPHI = 0.;

    if(Y != 0. || Z != 0.) {
        PHI = atan2(Z, Y);
        CPHI = Y/RHO;
        SPHI = Z/RHO;
    }

    const double CF = cos(PHI);
    const double SF = sin(PHI);

    double BY_AS1, BZ_AS1, BY_AS2, BZ_AS2;
    unwarped(X, RHO*CF, RHO*SF, DXSHIFT1, DXSHIFT2, D, BX1, &BY_AS1, &BZ_AS1, BX2, &BY_AS2, &BZ_AS2);

    double BRHO_AS = BY_AS1*CF + BZ_AS1*SF;
    double BPHI_AS = -BY_AS1*SF + BZ_AS1*CF;
    *BY1 = BRHO_AS*CPHI - BPHI_AS*SPHI;
    *BZ1 = BRHO_AS*SPHI + BPHI_AS*CPHI;

    BRHO_AS = BY_AS2*CF + BZ_AS2*SF;
    BPHI_AS = -BY_AS2*SF + BZ_AS2*CF;
    *BY2 = BRHO_AS*CPHI - BPHI_AS*SPHI;
    *BZ2 = BRHO_AS*SPHI + BPHI_AS*CPHI;
}

void unwarped(double X, double Y, double Z,
              double DXSHIFT1, double DXSHIFT2, double D0,
              double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2) {

    static const double A1[60] = {
        -25.45869857,  57.35899080, 317.5501869, -2.626756717,
        -93.38053698, -199.6467926, -858.8129729, 34.09192395, 845.4214929,
        -29.07463068, 47.10678547, -128.9797943, -781.7512093, 6.165038619,
        167.8905046, 492.0680410, 1654.724031, -46.77337920, -1635.922669,
        40.86186772, -.1349775602, -.9661991179E-01, -.1662302354,
        .002810467517, .2487355077, .1025565237, -14.41750229, -.8185333989,
        11.07693629, .7569503173, -9.655264745, 112.2446542, 777.5948964,
        -5.745008536, -83.03921993, -490.2278695, -1155.004209, 39.08023320,
        1172.780574, -39.44349797, -14.07211198, -40.41201127, -313.2277343,
        2.203920979, 8.232835341, 197.7065115, 391.2733948, -18.57424451,
        -437.2779053, 23.04976898, 11.75673963, 13.60497313, 4.691927060,
        18.20923547, 27.59044809, 6.677425469, 1.398283308, 2.839005878,
        31.24817706, 24.53577264
    };

    static const double A2[60] = {
        -287187.1962, 4970.499233, 410490.1952, -1347.839052,
        -386370.3240, 3317.983750, -143462.3895, 5706.513767, 171176.2904,
        250.8882750, -506570.8891, 5733.592632, 397975.5842, 9771.762168,
        -941834.2436, 7990.975260, 54313.10318, 447.5388060, 528046.3449,
        12751.04453, -21920.98301, -21.05075617, 31971.07875, 3012.641612,
        -301822.9103, -3601.107387, 1797.577552, -6.315855803, 142578.8406,
        13161.93640, 804184.8410, -14168.99698, -851926.6360, -1890.885671,
        972475.6869, -8571.862853, 26432.49197, -2554.752298, -482308.3431,
        -4391.473324, 105155.9160, -1134.622050, -74353.53091, -5382.670711,
        695055.0788, -916.3365144, -12111.06667, 67.20923358, -367200.9285,
        -21414.14421, 14.75567902, 20.75638190, 59.78601609, 16.86431444,
        32.58482365, 23.69472951, 17.24977936, 13.64902647, 68.40989058,
        11.67828167
    };

    double FX1, FY1, FZ1;
    taildisk(D0 * 1.1, 1., (X - 6. - DXSHIFT1) * 1.1 + 1.2, Y * 1.1, Z * 1.1, &FX1, &FY1, &FZ1);

    double HX1, HY1, HZ1;
    shlcar5x5(A1, X, Y, Z, DXSHIFT1, &HX1, &HY1, &HZ1);

    *BX1 = FX1 + HX1;
    *BY1 = FY1 + HY1;
    *BZ1 = FZ1 + HZ1;

    double FX2, FY2, FZ2;
    taildisk(D0*0.25, 0., (X - 4. - DXSHIFT2)*0.25 - 9., Y*0.25, Z*0.25, &FX2, &FY2, &FZ2);

    double HX2, HY2, HZ2;
    shlcar5x5(A2, X, Y, Z, DXSHIFT2, &HX2, &HY2, &HZ2);

    *BX2 = FX2 + HX2;
    *BY2 = FY2 + HY2;
    *BZ2 = FZ2 + HZ2;
}

void taildisk(double D0, double DELTADX,
              double X, double Y, double Z, double *BX, double *BY, double *BZ) {

    static const double F[5] = { -71.09346626, -1014.308601, -1272.939359, -3224.935936, -44546.86232 };
    static const double B[5] = { 10.90101242, 12.68393898, 13.51791954, 14.86775017, 15.12306404 };
    static const double C[5] = { .7954069972, .6716601849, 1.174866319, 2.565249920, 10.01986790 };

    const double RHO = sqrt(X*X + Y*Y);
    const double DRHODX = X/RHO;
    const double DRHODY = Y/RHO;

    const double DEX = exp(X/7.);
    const double Y_20 = Y / 20.;
    const double D = D0 + 4.7*Y_20*Y_20 + DELTADX*DEX;
    const double DDDY = 4.7*Y*0.005;
    const double DDDX = DELTADX/7.*DEX;

    const double DZETA = sqrt(Z*Z + D*D);

    const double DDZETADX = D*DDDX/DZETA;
    const double DDZETADY = D*DDDY/DZETA;
    const double DDZETADZ = Z/DZETA;

    double DBX = 0.;
    double DBY = 0.;
    double DBZ = 0.;
    double BI, CI, S1, S2, DS1DRHO, DS2DRHO, DS1DDZ, DS2DDZ, DS1DX, DS1DY, DS1DZ, DS2DX, DS2DY, DS2DZ,
            S1TS2, S1PS2, S1PS2SQ, FAC1, AS, DASDS1, DASDS2, DASDX, DASDY, DASDZ;
    
    for(int i = 0; i < 5; i++) {
        BI = B[i];
        CI = C[i];
        const double DZETA_CI = DZETA + CI;
        const double DZETACI2 = DZETA_CI*DZETA_CI;
        const double RHO_BI = RHO + BI;
        const double RHO_DIF_BI = RHO - BI;
        S1 = sqrt(RHO_BI*RHO_BI + DZETACI2);
        S2 = sqrt(RHO_DIF_BI*RHO_DIF_BI + DZETACI2);

        DS1DRHO = (RHO+BI) / S1;
        DS2DRHO = (RHO-BI) / S2;
        DS1DDZ = DZETA_CI / S1;
        DS2DDZ = DZETA_CI / S2;

        DS1DX = DS1DRHO*DRHODX + DS1DDZ*DDZETADX;
        DS1DY = DS1DRHO*DRHODY + DS1DDZ*DDZETADY;
        DS1DZ = DS1DDZ*DDZETADZ;
        DS2DX = DS2DRHO*DRHODX + DS2DDZ*DDZETADX;
        DS2DY = DS2DRHO*DRHODY + DS2DDZ*DDZETADY;
        DS2DZ = DS2DDZ*DDZETADZ;

        S1TS2 = S1*S2;
        S1PS2 = S1 + S2;
        S1PS2SQ = S1PS2*S1PS2;

        const double BI2 = 2. * BI;
        FAC1 = sqrt(S1PS2SQ - BI2*BI2);
        AS = FAC1/(S1TS2*S1PS2SQ);
        DASDS1 = (1./(FAC1*S2) - AS/S1PS2*(S2*S2 + S1*(3.*S1 + 4.*S2))) / (S1*S1PS2);
        DASDS2 = (1./(FAC1*S1) - AS/S1PS2*(S1*S1 + S2*(3.*S2 + 4.*S1))) / (S2*S1PS2);

        DASDX = DASDS1*DS1DX + DASDS2*DS2DX;
        DASDY = DASDS1*DS1DY + DASDS2*DS2DY;
        DASDZ = DASDS1*DS1DZ + DASDS2*DS2DZ;

        DBX = DBX - F[i]*X*DASDZ;
        DBY = DBY - F[i]*Y*DASDZ;
        DBZ = DBZ + F[i]*(2.*AS + X*DASDX + Y*DASDY);
    }

    *BX = DBX;
    *BY = DBY;
    *BZ = DBZ;
}

void shlcar5x5(const double *A, double X, double Y, double Z, double DSHIFT,
               double *HX, double *HY, double *HZ) {

    double RP, CYPI, SYPI, RR, SZRK, CZRK, SQPR, EPR, DBX, DBY, DBZ, COEF;
    int L = 0;
    double DHX = 0.;
    double DHY = 0.;
    double DHZ = 0.;
    for(int i = 1; i <= 5; i++) {
        RP = 1./A[49+i];
        const double YRP = Y*RP;
        double SQPR_PR = RP*RP;
        CYPI = cos(YRP);
        SYPI = sin(YRP);

        for(int k = 1; k <= 5; k++) {
            RR = 1. / A[54+k];
            SQPR = sqrt(SQPR_PR + RR*RR);
            EPR = exp(X * SQPR);
            const double ZRR = Z * RR;
            SZRK = sin(ZRR);
            CZRK = cos(ZRR);

            DBX = -SQPR * EPR * CYPI * SZRK;
            DBY = RP * EPR * SYPI * SZRK;
            DBZ = -RR * EPR * CYPI * CZRK;

            L += 2;
            COEF = A[L-2] + A[L-1] * DSHIFT;

            DHX = DHX + COEF * DBX;
            DHY = DHY + COEF * DBY;
            DHZ = DHZ + COEF * DBZ;
        }
    }

    *HX = DHX;
    *HY = DHY;
    *HZ = DHZ;
}
