#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

int string_to_uint64(char *str_num, uint64_t *num) {
    char *end_ptr;
    *num = strtoul(str_num, &end_ptr, 10);
    if(end_ptr == str_num || *end_ptr != '\0' || ((*num == LONG_MIN || *num == LONG_MAX) && errno == ERANGE)) {
        return -1;
    }

    return 0;
}

int string_to_double(char *str_num, double *num) {
    char *end_ptr;
    *num = strtod(str_num, &end_ptr);
    if(end_ptr == str_num || *end_ptr != '\0' || ((*num == LONG_MIN || *num == LONG_MAX) && errno == ERANGE)) {
        return -1;
    }

    return 0;
}

void sph_to_car(Vector *sph, Vector *car) {
    const double r = sph->r;
    const double theta = sph->theta;
    const double phi = sph->phi;

    const double SQ = r * sin(theta);
    const double x = SQ * cos(phi);
    const double y = SQ * sin(phi);
    const double z = r * cos(theta);

    car->x = x;
    car->y = y;
    car->z = z;
}

void car_to_sph(Vector *car, Vector *sph) {
    const double x = car->x;
    const double y = car->y;
    const double z = car->z;

    double SQ = x*x + y*y;
    double r = sqrt(SQ + z*z);
    
    if(SQ != 0.) {
        double phi = atan2(y, x);
        double theta = atan2(sqrt(SQ), z);
        
        if(phi < 0.) {
            phi += DEF_2PI;
        }

        sph->r = r;
        sph->theta = theta;
        sph->phi = phi;
        return;
    }
    
    double theta = 0.;
    if(z < 0.) {
        theta = DEF_PI;
    }

    sph->r = r;
    sph->theta = theta;
    sph->phi = 0.;
}

void geo_to_gsm(Vector *geo, Vector *gsm, const double *A1, const double *A2, const double *A3) {
    const double x = geo->x;
    const double y = geo->y;
    const double z = geo->z;

    const double nx = A1[0]*x + A1[1]*y + A1[2]*z;
    const double ny = A2[0]*x + A2[1]*y + A2[2]*z;
    const double nz = A3[0]*x + A3[1]*y + A3[2]*z;

    gsm->x = nx;
    gsm->y = ny;
    gsm->z = nz;
}

void gsm_to_geo(Vector *gsm, Vector *geo, const double *A1, const double *A2, const double *A3) {
    const double x = gsm->x;
    const double y = gsm->y;
    const double z = gsm->z;
    
    const double nx = A1[0]*x + A2[0]*y + A3[0]*z;
    const double ny = A1[1]*x + A2[1]*y + A3[1]*z;
    const double nz = A1[2]*x + A2[2]*y + A3[2]*z;

    geo->x = nx;
    geo->y = ny;
    geo->z = nz;
}
