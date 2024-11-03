#ifndef T05_H
#define T05_H

#include <stdio.h>

void ts05_gsm(const double *PARMOD, double XGSM, double YGSM, double ZGSM, double *BX, double *BY, double *BZ);

void shlcar3x3(double X, double Y, double Z, double *BX, double *BY, double *BZ);

void dipole(double X, double Y, double Z, double *BX, double *BY, double *BZ);

#endif // T05_H
