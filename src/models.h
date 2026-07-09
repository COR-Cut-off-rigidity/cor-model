#ifndef MODELS_H
#define MODELS_H

#include <stdio.h>
#include <stdint.h>
#include <internal.h>
#include <external.h>

#include "common.h"

typedef struct {
    InternalModelContext internal_model_context;
    ExternalModelContext external_model_context;
    ExternalGeoGsmMatrix external_geogsm_matrix;
} ModelsParams;

typedef struct {
    const char *name;
    const char *version;
    const void *setup_data;
    int (*setup_fun) (FILE* input_file, const InfileHeader *infile_header, ModelsParams *model_params, const void *setup_data);
    void (*calc_fun) (Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params);
    void (*destroy_fun) (ModelsParams *model_params);
} Model;

extern const Model *external_models[];
extern const Model *internal_models[];
extern const Model *crustal_models[];

const Model* get_model(const Model **models, const char *version);
void calculate_field(Vector gsm_pos, Vector *total_geo_field, const Model **models, const ModelsParams *params);

#endif // MODELS_H
