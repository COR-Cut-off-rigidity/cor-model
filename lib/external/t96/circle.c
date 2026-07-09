#include <math.h>
#include "circle.h"

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
