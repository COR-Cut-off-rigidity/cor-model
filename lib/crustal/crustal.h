#ifndef CHAOS_H
#define CHAOS_H

#include <crustal_coefs.h>

#define DEF_RE_KM_D 6371.2
#define DEF_2PI_D 6.283185307179586
#define DEF_PI_D 3.141592653589793

extern const CrustalModelCoefs CHAOS_CRUSTAL_7_16_COEFS;
extern const CrustalModelCoefs CHAOS_CRUSTAL_7_17_COEFS;
extern const CrustalModelCoefs CHAOS_CRUSTAL_7_18_COEFS;
extern const CrustalModelCoefs CHAOS_CRUSTAL_8_1_COEFS;
extern const CrustalModelCoefs CHAOS_CRUSTAL_8_2_COEFS;
extern const CrustalModelCoefs CHAOS_CRUSTAL_8_3_COEFS;

void chaos_geo(double XGEO, double YGEO, double ZGEO, const CrustalModelCoefs *coefs, double *BXGEO, double *BYGEO, double *BZGEO);

#endif // CHAOS_H
