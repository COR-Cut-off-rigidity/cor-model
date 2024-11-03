#ifndef MODELS_H
#define MODELS_H

#include <stdio.h>
#include <stdint.h>
#include <internal/internal.h>
#include "common.h"

typedef struct {
    double PARMOD[10];
    double G[INTERNAL_COEFS_SIZE];
    double H[INTERNAL_COEFS_SIZE];
    double REC[INTERNAL_COEFS_SIZE];
    double A1[3];
    double A2[3];
    double A3[3];
} ModelsParams;

typedef struct {
    const char *name;
    const char *version;
    const void *setup_data;
    int (*setup_fun) (FILE* input_file, const InfileHeader *infile_header, ModelsParams *model_params, const void *setup_data);
    void (*calc_fun) (Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params);
} Model;

extern const Model *external_models[];
extern const Model *internal_models[];
extern const Model *crustal_models[];

const Model* get_model(const Model **models, const char *version);
void calculate_field(Vector gsm_pos, Vector *total_geo_field, const Model **models, ModelsParams *params);

#endif // MODELS_H
