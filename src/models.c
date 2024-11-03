#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "common.h"
#include <ts05/ts05.h>
#include <t96/t96.h>
#include <chaos/chaos.h>
#include "models.h"
#include "io.h"
#include "utils.h"

int internal_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data);
int ts05_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data);
int t96_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data);
int chaos_7_16_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data);

void ts05_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void t96_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void internal_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void chaos_7_16_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);

static const Model igrf_9 = {
    .name = "IGRF-9",
    .version = "igrf9",
    .setup_data = &IGRF_9_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc
};

static const Model igrf_10 = {
    .name = "IGRF-10",
    .version = "igrf10",
    .setup_data = &IGRF_10_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc
};

static const Model igrf_11 = {
    .name = "IGRF-11",
    .version = "igrf11",
    .setup_data = &IGRF_11_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc
};

static const Model igrf_12 = {
    .name = "IGRF-12",
    .version = "igrf12",
    .setup_data = &IGRF_12_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc
};

static const Model igrf_13 = {
    .name = "IGRF-13",
    .version = "igrf13",
    .setup_data = &IGRF_13_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc
};

static const Model hist = {
    .name = "Historical (0-1968 CE)",
    .version = "hist",
    .setup_data = &HIST_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc
};

static const Model cals10k_2 = {
    .name = "Historical (CALS10k.2)",
    .version = "cals10k.2",
    .setup_data = &CALS10K_2_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc
};

static const Model ts_05 = {
    .name = "Tsyganenko-Sitnov 05",
    .version = "ts05",
    .setup_data = NULL,
    .setup_fun = ts05_setup,
    .calc_fun = ts05_calc
};

static const Model t_96 = {
    .name = "Tsyganenko 96",
    .version = "t96",
    .setup_data = NULL,
    .setup_fun = t96_setup,
    .calc_fun = t96_calc
};

static const Model chaos_7_16 = {
    .name = "CHAOS-7.16",
    .version = "chaos7.16",
    .setup_data = NULL,
    .setup_fun = NULL,
    .calc_fun = chaos_7_16_calc
};

const Model *internal_models[] = { &igrf_13, &igrf_12, &igrf_11, &igrf_10, &igrf_9, &hist, &cals10k_2, NULL };
const Model *external_models[] = { &ts_05, &t_96, NULL };
const Model *crustal_models[] = { &chaos_7_16, NULL };

const Model* get_model(const Model **models, const char *version) {
    for(int i = 0; models[i] != NULL; i++) {
        if(strcmp(models[i]->version, version) == 0) {
            return models[i];
        }
    }

    return NULL;
}

int internal_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data) {
    const RecalcTime time = {
        .iy = infile_header->iy,
        .mes = infile_header->mes,
        .ide = infile_header->ide,
        .id = infile_header->id,
        .ih = infile_header->ih,
        .min = infile_header->min,
        .is = infile_header->is
    };

    return internal_recalc((InternalModelCoefs*) setup_data, &time, params->G, params->H, params->REC, params->A1, params->A2, params->A3);
}

int ts05_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data) {
    return read_ts05_params(input_file, params);
}

int t96_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data) {
    return read_t96_params(input_file, params);
}

void internal_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    gsm_to_geo(&gsm_pos, &geo_pos, params->A1, params->A2, params->A3);
    
    Vector geo_field;
    internal_geo(geo_pos.x, geo_pos.y, geo_pos.z, params->G, params->H, params->REC, &geo_field.x, &geo_field.y, &geo_field.z);
    //internal_geo2(pos.x, pos.y, pos.z, params->G, params->H, &HX, &HY, &HZ);

    Vector gsm_field;
    geo_to_gsm(&geo_field, &gsm_field, params->A1, params->A2, params->A3);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void ts05_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    double BBX, BBY, BBZ;
    
    ts05_gsm(params->PARMOD, gsm_pos.x, gsm_pos.y, gsm_pos.z, &BBX, &BBY, &BBZ);
    total_gsm_field->x += BBX;
    total_gsm_field->y += BBY;
    total_gsm_field->z += BBZ;
}

void t96_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    double BBX, BBY, BBZ;
    t96_gsm(params->PARMOD, gsm_pos.x, gsm_pos.y, gsm_pos.z, &BBX, &BBY, &BBZ);
    total_gsm_field->x += BBX;
    total_gsm_field->y += BBY;
    total_gsm_field->z += BBZ;
}

void chaos_7_16_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    gsm_to_geo(&gsm_pos, &geo_pos, params->A1, params->A2, params->A3);
    
    Vector geo_field;
    chaos_geo(geo_pos.x, geo_pos.y, geo_pos.z, &CHAOS_7_16_COEFS, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    geo_to_gsm(&geo_field, &gsm_field, params->A1, params->A2, params->A3);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void calculate_field(Vector gsm_pos, Vector *total_geo_field, const Model **models, ModelsParams *params) {
    Vector gsm_pos_re = {
        .x = gsm_pos.x / DEF_RE,
        .y = gsm_pos.y / DEF_RE,
        .z = gsm_pos.z / DEF_RE
    };

    Vector total_gsm_field = {0};
    for(int i = 0; models[i]; i++) {
        models[i]->calc_fun(gsm_pos_re, &total_gsm_field, params);
    }

    gsm_to_geo(&total_gsm_field, total_geo_field, params->A1, params->A2, params->A3);
}
