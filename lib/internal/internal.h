#ifndef INTERNAL_H
#define INTERNAL_H

#include "coefs/coefs.h"

#define DEF_TO_DEG 57.295779513 // (180./PI)
//#define DEF_RE_KM_DF 6371.2
//#define DEF_2PI_DF 6.283185307179586
//#define DEF_PI_DF 3.141592653589793
//#define NM_MAX 14

typedef struct {
    int iy;
    int mes;
    int ide;
    int id;
    int ih;
    int min;
    int is;
} RecalcTime;

extern const InternalModelCoefs IGRF_9_COEFS;
extern const InternalModelCoefs IGRF_10_COEFS;
extern const InternalModelCoefs IGRF_11_COEFS;
extern const InternalModelCoefs IGRF_12_COEFS;
extern const InternalModelCoefs IGRF_13_COEFS;
extern const InternalModelCoefs HIST_COEFS;
extern const InternalModelCoefs CALS10K_2_COEFS;

int internal_recalc(const InternalModelCoefs *coefs, const RecalcTime *time, double *G, double *H, double *REC, double *A1, double *A2, double *A3);
void internal_geo(double XGEO, double YGEO, double ZGEO, const double *G, const double *H, const double *REC, double *BXGEO, double *BYGEO, double *BZGEO);
//void internal_geo2(double XGEO, double YGEO, double ZGEO, const double *G, const double *H, double *BXGEO, double *BYGEO, double *BZGEO);

#endif // INTERNAL_H
