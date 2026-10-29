#ifndef INTERNAL_H
#define INTERNAL_H

#include <internal_coefs.h>
#include <calc_time.h>

typedef struct {
    double *coefs;
    uint32_t nm_max;
} InternalModelContext;

extern const InternalModelCoefs IGRF_9_COEFS;
extern const InternalModelCoefs IGRF_10_COEFS;
extern const InternalModelCoefs IGRF_11_COEFS;
extern const InternalModelCoefs IGRF_12_COEFS;
extern const InternalModelCoefs IGRF_13_COEFS;
extern const InternalModelCoefs IGRF_14_COEFS;
extern const InternalModelCoefs HISTORICAL_COEFS;
extern const InternalModelCoefs CALS10K_2_COEFS;
extern const InternalModelCoefs CHAOS_CORE_7_16_COEFS;
extern const InternalModelCoefs CHAOS_CORE_7_17_COEFS;
extern const InternalModelCoefs CHAOS_CORE_7_18_COEFS;
extern const InternalModelCoefs CHAOS_CORE_8_1_COEFS;
extern const InternalModelCoefs CHAOS_CORE_8_2_COEFS;
extern const InternalModelCoefs CHAOS_CORE_8_3_COEFS;

int internal_recalc(const InternalModelCoefs *model_coefs, const CalcTime *time, InternalModelContext *internal_model_context);
void internal_geo(const InternalModelContext *internal_model_context, double XGEO, double YGEO, double ZGEO, double *BXGEO, double *BYGEO, double *BZGEO);
void internal_clear(InternalModelContext *internal_model_context);

#endif // INTERNAL_H
