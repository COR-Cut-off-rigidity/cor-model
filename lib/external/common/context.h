#ifndef CONTEXT_H
#define CONTEXT_H

typedef struct {
    double pdyn;
    double dst;
    double by_imf;
    double bz_imf;
    double W[6]; // for TS05
} ExternalModelContext;

#endif // CONTEXT_H
