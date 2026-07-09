#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <crustal.h>
#include <ts05/ts05.h>
#include <t96/t96.h>

#include "models.h"
#include "io.h"
#include "utils.h"

int internal_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data);
int ts05_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data);
int t96_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data);

void ts05_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void t96_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void internal_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void chaos_crustal_7_16_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void chaos_crustal_7_17_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void chaos_crustal_7_18_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void chaos_crustal_8_1_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void chaos_crustal_8_2_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);
void chaos_crustal_8_3_calc(Vector gsm_pos, Vector *field, const ModelsParams *params);

void internal_destroy(ModelsParams *params);

static const Model igrf_9 = {
    .name = "IGRF-9",
    .version = "igrf9",
    .setup_data = &IGRF_9_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model igrf_10 = {
    .name = "IGRF-10",
    .version = "igrf10",
    .setup_data = &IGRF_10_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model igrf_11 = {
    .name = "IGRF-11",
    .version = "igrf11",
    .setup_data = &IGRF_11_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model igrf_12 = {
    .name = "IGRF-12",
    .version = "igrf12",
    .setup_data = &IGRF_12_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model igrf_13 = {
    .name = "IGRF-13",
    .version = "igrf13",
    .setup_data = &IGRF_13_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model igrf_14 = {
    .name = "IGRF-14",
    .version = "igrf14",
    .setup_data = &IGRF_14_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model chaos_core_7_16 = {
    .name = "CHAOS-7.16 (core field)",
    .version = "chaos7.16",
    .setup_data = &CHAOS_CORE_7_16_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model chaos_core_7_17 = {
    .name = "CHAOS-7.17 (core field)",
    .version = "chaos7.17",
    .setup_data = &CHAOS_CORE_7_17_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model chaos_core_7_18 = {
    .name = "CHAOS-7.18 (core field)",
    .version = "chaos7.18",
    .setup_data = &CHAOS_CORE_7_18_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model chaos_core_8_1 = {
    .name = "CHAOS-8.1 (core field)",
    .version = "chaos8.1",
    .setup_data = &CHAOS_CORE_8_1_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model chaos_core_8_2 = {
    .name = "CHAOS-8.2 (core field)",
    .version = "chaos8.2",
    .setup_data = &CHAOS_CORE_8_2_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model chaos_core_8_3 = {
    .name = "CHAOS-8.3 (core field)",
    .version = "chaos8.3",
    .setup_data = &CHAOS_CORE_8_3_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model historical = {
    .name = "Historical (0-1968 CE)",
    .version = "hist",
    .setup_data = &HISTORICAL_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model cals10k_2 = {
    .name = "CALS10k.2",
    .version = "cals10k.2",
    .setup_data = &CALS10K_2_COEFS,
    .setup_fun = internal_setup,
    .calc_fun = internal_calc,
    .destroy_fun = internal_destroy
};

static const Model ts_05 = {
    .name = "Tsyganenko-Sitnov 05",
    .version = "ts05",
    .setup_data = NULL,
    .setup_fun = ts05_setup,
    .calc_fun = ts05_calc,
    .destroy_fun = NULL
};

static const Model t_96 = {
    .name = "Tsyganenko 96",
    .version = "t96",
    .setup_data = NULL,
    .setup_fun = t96_setup,
    .calc_fun = t96_calc,
    .destroy_fun = NULL
};

static const Model chaos_crustal_7_16 = {
    .name = "CHAOS-7.16 (crustal field)",
    .version = "chaos7.16",
    .setup_data = NULL,
    .setup_fun = NULL,
    .calc_fun = chaos_crustal_7_16_calc,
    .destroy_fun = NULL
};

static const Model chaos_crustal_7_17 = {
    .name = "CHAOS-7.17 (crustal field)",
    .version = "chaos7.17",
    .setup_data = NULL,
    .setup_fun = NULL,
    .calc_fun = chaos_crustal_7_17_calc,
    .destroy_fun = NULL
};

static const Model chaos_crustal_7_18 = {
    .name = "CHAOS-7.18 (crustal field)",
    .version = "chaos7.18",
    .setup_data = NULL,
    .setup_fun = NULL,
    .calc_fun = chaos_crustal_7_18_calc,
    .destroy_fun = NULL
};

static const Model chaos_crustal_8_1 = {
    .name = "CHAOS-8.1 (crustal field)",
    .version = "chaos8.1",
    .setup_data = NULL,
    .setup_fun = NULL,
    .calc_fun = chaos_crustal_8_1_calc,
    .destroy_fun = NULL
};

static const Model chaos_crustal_8_2 = {
    .name = "CHAOS-8.2 (crustal field)",
    .version = "chaos8.2",
    .setup_data = NULL,
    .setup_fun = NULL,
    .calc_fun = chaos_crustal_8_2_calc,
    .destroy_fun = NULL
};

static const Model chaos_crustal_8_3 = {
    .name = "CHAOS-8.3 (crustal field)",
    .version = "chaos8.3",
    .setup_data = NULL,
    .setup_fun = NULL,
    .calc_fun = chaos_crustal_8_3_calc,
    .destroy_fun = NULL
};

const Model *internal_models[] = { &igrf_14, &igrf_13, &igrf_12, &igrf_11, &igrf_10, &igrf_9, &chaos_core_8_3, &chaos_core_8_2, &chaos_core_8_1, &chaos_core_7_18, &chaos_core_7_17, &chaos_core_7_16, &historical, &cals10k_2, NULL };
const Model *external_models[] = { &ts_05, &t_96, NULL };
const Model *crustal_models[] = { &chaos_crustal_8_3, &chaos_crustal_8_2, &chaos_crustal_8_1, &chaos_crustal_7_18, &chaos_crustal_7_17, &chaos_crustal_7_16, NULL };

const Model* get_model(const Model **models, const char *version) {
    for (int i = 0; models[i] != NULL; i++) {
        if (strcmp(models[i]->version, version) == 0) {
            return models[i];
        }
    }

    return NULL;
}

int internal_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data) {
    int result = internal_recalc((InternalModelCoefs*) setup_data, &infile_header->time, &params->internal_model_context);
    if (result == -1) {
        return result;
    }

    const double G1_0 = params->internal_model_context.coefs[3];
    const double G1_1 = params->internal_model_context.coefs[6];
    const double H1_1 = params->internal_model_context.coefs[7];
    external_init_transformation_matrix(&infile_header->time, G1_0, G1_1, H1_1, &params->external_geogsm_matrix);

    return result;
}

int ts05_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data) {
    int result = read_ts05_params(input_file, &params->external_model_context);
    if (result != 0) {
        return result;
    }

    return 0;
}

int t96_setup(FILE* input_file, const InfileHeader *infile_header, ModelsParams *params, const void *setup_data) {
    return read_t96_params(input_file, &params->external_model_context);
}

void internal_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    external_gsm_to_geo(&params->external_geogsm_matrix, gsm_pos.x, gsm_pos.y, gsm_pos.z, &geo_pos.x, &geo_pos.y, &geo_pos.z);
    
    Vector geo_field;
    internal_geo(&params->internal_model_context, geo_pos.x, geo_pos.y, geo_pos.z, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    external_geo_to_gsm(&params->external_geogsm_matrix, geo_field.x, geo_field.y, geo_field.z, &gsm_field.x, &gsm_field.y, &gsm_field.z);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void ts05_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    double BBX, BBY, BBZ;
    
    ts05_gsm(&params->external_model_context, gsm_pos.x, gsm_pos.y, gsm_pos.z, &BBX, &BBY, &BBZ);
    total_gsm_field->x += BBX;
    total_gsm_field->y += BBY;
    total_gsm_field->z += BBZ;
}

void t96_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    double BBX, BBY, BBZ;

    t96_gsm(&params->external_model_context, gsm_pos.x, gsm_pos.y, gsm_pos.z, &BBX, &BBY, &BBZ);
    total_gsm_field->x += BBX;
    total_gsm_field->y += BBY;
    total_gsm_field->z += BBZ;
}

void chaos_crustal_7_16_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    external_gsm_to_geo(&params->external_geogsm_matrix, gsm_pos.x, gsm_pos.y, gsm_pos.z, &geo_pos.x, &geo_pos.y, &geo_pos.z);
    
    Vector geo_field;
    chaos_geo(geo_pos.x, geo_pos.y, geo_pos.z, &CHAOS_CRUSTAL_7_16_COEFS, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    external_geo_to_gsm(&params->external_geogsm_matrix, geo_field.x, geo_field.y, geo_field.z, &gsm_field.x, &gsm_field.y, &gsm_field.z);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void chaos_crustal_7_17_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    external_gsm_to_geo(&params->external_geogsm_matrix, gsm_pos.x, gsm_pos.y, gsm_pos.z, &geo_pos.x, &geo_pos.y, &geo_pos.z);
    
    Vector geo_field;
    chaos_geo(geo_pos.x, geo_pos.y, geo_pos.z, &CHAOS_CRUSTAL_7_17_COEFS, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    external_geo_to_gsm(&params->external_geogsm_matrix, geo_field.x, geo_field.y, geo_field.z, &gsm_field.x, &gsm_field.y, &gsm_field.z);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void chaos_crustal_7_18_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    external_gsm_to_geo(&params->external_geogsm_matrix, gsm_pos.x, gsm_pos.y, gsm_pos.z, &geo_pos.x, &geo_pos.y, &geo_pos.z);
    
    Vector geo_field;
    chaos_geo(geo_pos.x, geo_pos.y, geo_pos.z, &CHAOS_CRUSTAL_7_18_COEFS, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    external_geo_to_gsm(&params->external_geogsm_matrix, geo_field.x, geo_field.y, geo_field.z, &gsm_field.x, &gsm_field.y, &gsm_field.z);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void chaos_crustal_8_1_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    external_gsm_to_geo(&params->external_geogsm_matrix, gsm_pos.x, gsm_pos.y, gsm_pos.z, &geo_pos.x, &geo_pos.y, &geo_pos.z);
    
    Vector geo_field;
    chaos_geo(geo_pos.x, geo_pos.y, geo_pos.z, &CHAOS_CRUSTAL_8_1_COEFS, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    external_geo_to_gsm(&params->external_geogsm_matrix, geo_field.x, geo_field.y, geo_field.z, &gsm_field.x, &gsm_field.y, &gsm_field.z);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void chaos_crustal_8_2_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    external_gsm_to_geo(&params->external_geogsm_matrix, gsm_pos.x, gsm_pos.y, gsm_pos.z, &geo_pos.x, &geo_pos.y, &geo_pos.z);
    
    Vector geo_field;
    chaos_geo(geo_pos.x, geo_pos.y, geo_pos.z, &CHAOS_CRUSTAL_8_2_COEFS, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    external_geo_to_gsm(&params->external_geogsm_matrix, geo_field.x, geo_field.y, geo_field.z, &gsm_field.x, &gsm_field.y, &gsm_field.z);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void chaos_crustal_8_3_calc(Vector gsm_pos, Vector *total_gsm_field, const ModelsParams *params) {
    Vector geo_pos;
    external_gsm_to_geo(&params->external_geogsm_matrix, gsm_pos.x, gsm_pos.y, gsm_pos.z, &geo_pos.x, &geo_pos.y, &geo_pos.z);
    
    Vector geo_field;
    chaos_geo(geo_pos.x, geo_pos.y, geo_pos.z, &CHAOS_CRUSTAL_8_3_COEFS, &geo_field.x, &geo_field.y, &geo_field.z);

    Vector gsm_field;
    external_geo_to_gsm(&params->external_geogsm_matrix, geo_field.x, geo_field.y, geo_field.z, &gsm_field.x, &gsm_field.y, &gsm_field.z);

    total_gsm_field->x += gsm_field.x;
    total_gsm_field->y += gsm_field.y;
    total_gsm_field->z += gsm_field.z;
}

void internal_destroy(ModelsParams *model_params) {
    internal_clear(&model_params->internal_model_context);
}

void calculate_field(Vector gsm_pos, Vector *total_geo_field, const Model **models, const ModelsParams *params) {
    Vector gsm_pos_re = {
        .x = gsm_pos.x / DEF_RE,
        .y = gsm_pos.y / DEF_RE,
        .z = gsm_pos.z / DEF_RE
    };

    Vector total_gsm_field = {0};
    for (int i = 0; models[i]; i++) {
        models[i]->calc_fun(gsm_pos_re, &total_gsm_field, params);
    }

    external_gsm_to_geo(&params->external_geogsm_matrix, total_gsm_field.x, total_gsm_field.y, total_gsm_field.z, &total_geo_field->x, &total_geo_field->y, &total_geo_field->z);
}
