#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>
#include "common.h"

int string_to_uint64(char *str_num, uint64_t *num);
int string_to_double(char *str_num, double *num);
void sph_to_car(Vector *sph, Vector *car);
void car_to_sph(Vector *car, Vector *sph);

#endif // UTILS_H
