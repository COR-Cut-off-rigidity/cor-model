#ifndef CHAOS_COEFS
#define CHAOS_COEFS

#include <stdint.h>

typedef struct {
    const void *G;
    const void *H;
    uint32_t nm_max;
} CHAOSCoefs;

#endif // CHAOS_COEFS
