#ifndef T96_H
#define T96_H

void t96_gsm(const double *PARMOD, double X, double Y, double Z, double *BX, double *BY, double *BZ);
void cylharm(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void cylhar1(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void intercon(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void tailrc_96(double X, double Y, double Z, double *BXRC, double *BYRC, double *BZRC, double *BXT2, double *BYT2, double *BZT2, double *BXT3, double *BYT3, double *BZT3);
void ringcurr_96(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void taildisk_t96(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void tail87(double X, double Z, double *BX, double *BZ);
void shlcar3x3_t96(const double *A, double X, double Y, double Z, double *HX, double *HY, double *HZ);
void birk1tot_02(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void diploop1(double X, double Y, double Z, double D[3][14]);
void circle(double X, double Y, double Z, double RL, double *BX, double *BY, double *BZ);
void crosslp(double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double RL, double AL);
void dipxyz(double X, double Y, double Z, double *BXX, double *BYX, double *BZX, double *BXY, double *BYY, double *BZY, double *BXZ, double *BYZ, double *BZZ);
void condip1(double X, double Y, double Z, double D[3][42]);
void birk1shld(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void birk2shl(double X, double Y, double Z, double *HX, double *HY, double *HZ);
void r2_birk(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void r2inner(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void bconic(double X, double Y, double Z, double *CBX, double *CBY, double *CBZ);
void dipdistr(double X, double Y, double Z, double *BX, double *BY, double *BZ, int MODE);
void r2outer(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void loops4(double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double YC, double ZC, double R, double THETA, double PHI);
void r2sheet(double X, double Y, double Z, double *BX, double *BY, double *BZ);
void dipole_t96(double X, double Y, double Z, double *BX, double *BY, double *BZ);
double tksi(double XKSI, double XKS0, double DXKSI);
double xksi(double X, double Y, double Z);
double bes0(double X);
double bes1(double X);

#endif // T96_H
