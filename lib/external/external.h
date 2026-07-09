#ifndef EXTERNAL_H
#define EXTERNAL_H

#include <calc_time.h>
#include <context.h>

typedef struct {
    double A1[3];
    double A2[3];
    double A3[3];
} ExternalGeoGsmMatrix;

void external_init_transformation_matrix(const CalcTime *calc_time, double G1_0, double G1_1, double H1_1, ExternalGeoGsmMatrix *mtx);
void external_geo_to_gsm(const ExternalGeoGsmMatrix *mtx, double geo_x, double geo_y, double geo_z, double *gsm_x, double *gsm_y, double *gsm_z);
void external_gsm_to_geo(const ExternalGeoGsmMatrix *mtx, double gsm_x, double gsm_y, double gsm_z, double *geo_x, double *geo_y, double *geo_z);

#endif // EXTERNAL_H
