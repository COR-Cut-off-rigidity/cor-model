#ifndef CHAOS_H
#define CHAOS_H

#include "coefs/chaos_coefs.h"

#define DEF_RE_KM_D 6371.2
#define DEF_2PI_D 6.283185307179586
#define DEF_PI_D 3.141592653589793

extern const CHAOSCoefs CHAOS_7_16_COEFS;

void chaos_geo(double XGEO, double YGEO, double ZGEO, const CHAOSCoefs *coefs, double *BXGEO, double *BYGEO, double *BZGEO);

#endif // CHAOS_H
