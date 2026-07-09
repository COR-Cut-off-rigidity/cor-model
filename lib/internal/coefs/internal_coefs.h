#ifndef COEFS_H
#define COEFS_H

#include <stdint.h>
#include <stdbool.h>

typedef struct {
    const double *G;
    const double *H;
    double year;
} GHTime;

typedef struct {
    const GHTime **coefs;
    int coefs_size;
    int nm_max;
    bool extrapolate;
} InternalModelCoefs;

#endif // COEFS_H