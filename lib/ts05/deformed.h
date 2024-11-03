#ifndef DEFORMED_H
#define DEFORMED_H

void deformed(double X, double Y, double Z, double DXSHIFT1, double DXSHIFT2, double D, double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2);

void warped(double X, double Y, double Z, double DXSHIFT1, double DXSHIFT2, double D, double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2);

void unwarped(double X, double Y, double Z, double DXSHIFT1, double DXSHIFT2, double D0, double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2);

void taildisk(double D0, double DELTADX, double X, double Y, double Z, double *BX, double *BY, double *BZ);

void shlcar5x5(const double *A, double X, double Y, double Z, double DSHIFT, double *HX, double *HY, double *HZ);

#endif // DEFORMED_H
