#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

int string_to_uint64(char *str_num, uint64_t *num) {
    char *end_ptr;
    *num = strtoul(str_num, &end_ptr, 10);
    if (end_ptr == str_num || *end_ptr != '\0' || ((*num == ULONG_MAX ) && errno == ERANGE)) {
        return -1;
    }

    return 0;
}

int string_to_double(char *str_num, double *num) {
    char *end_ptr;
    *num = strtod(str_num, &end_ptr);
    if (end_ptr == str_num || *end_ptr != '\0' || ((*num == HUGE_VAL || *num == -HUGE_VAL) && errno == ERANGE)) {
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
    
    if (SQ != 0.) {
        double phi = atan2(y, x);
        double theta = atan2(sqrt(SQ), z);
        
        if (phi < 0.) {
            phi += DEF_2PI;
        }

        sph->r = r;
        sph->theta = theta;
        sph->phi = phi;
        return;
    }
    
    double theta = 0.;
    if (z < 0.) {
        theta = DEF_PI;
    }

    sph->r = r;
    sph->theta = theta;
    sph->phi = 0.;
}
