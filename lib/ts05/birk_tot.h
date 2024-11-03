#ifndef BIRK_TOT_H
#define BIRK_TOT_H

void birk_tot(double X, double Y, double Z, double XKAPPA1, double XKAPPA2, double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12, double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22);

void birk_1n2(int NUMB, int MODE, double X, double Y, double Z, double XKAPPA, double *BX, double *BY, double *BZ);

void twocones(const double *A, double X, double Y, double Z, int MODE, double DTHETA, double *BX, double *BY, double *BZ);

void one_cone(const double *A, double X, double Y, double Z, int MODE, double DTHETA, double *BX, double *BY, double *BZ);

double r_s(const double *A, double R, double THETA);

double theta_s(const double *A, double R, double THETA);

void fialcos(double R, double THETA, double PHI, int MODE, double THETA0, double DT, double *BTHETA, double *BPHI);

void birk_shl(const double *A, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ);

#endif // BIRK_TOT_H
