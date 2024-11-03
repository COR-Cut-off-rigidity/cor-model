#ifndef COEFS_H
#define COEFS_H

#include <stdint.h>
#include <stdbool.h>

#define INTERNAL_COEFS_SIZE 105

typedef struct {
    const double **coefs;
    uint32_t coefs_size;
    uint32_t years_step;
    uint32_t min_year;
    uint32_t max_year;
    bool extrapolate;
} InternalModelCoefs;

#endif // COEFS_H