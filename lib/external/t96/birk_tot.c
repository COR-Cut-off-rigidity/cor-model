#include <math.h>
#include "birk_tot.h"
#include "circle.h"

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
        if (fabs(XX_DIPY) > 1.0e-10) {
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
