#ifndef FULL_RC_H
#define FULL_RC_H

void full_rc(double X, double Y, double Z, double PHI, double SC_SY, double SC_PR, double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC, double *BZPRC);

void src_prc(double SC_SY, double SC_PR, double PHI, double X, double Y, double Z, double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC, double *BZPRC);

void rc_symm(double X, double Y, double Z, double *BX, double *BY, double *BZ);

double ap(double R, double SINT, double COST);

void prc_symm(double X, double Y, double Z, double *BX, double *BY, double *BZ);

double apprc(double R, double SINT, double COST);

void prc_quad(double X, double Y, double Z, double *BX, double *BY, double *BZ);

double br_prc_q(double R, double SINT, double COST);

double bt_prc_q(double R, double SINT, double COST);

void calc_ffs(double A, double A0, double DA, double *F, double *FA, double *FS);

void rc_shield(const double *A, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ);

#endif // FULL_RC_H
