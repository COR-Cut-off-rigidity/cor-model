#ifndef ARGS_H
#define ARGS_H

#include <stdint.h>

#define OPT_PARALLEL 1

typedef struct {
    char *input_file_name;
    char *output_file_name;
    char *external_model_name;
    char *internal_model_name;
    char *crustal_model_name;
    char *steps_count;
    char *tu_angle;
    char *mode;
    char *particle_type;
    uint32_t options;
} CLIArguments;

void args_parse(int argc, char **argv, CLIArguments *arguments);

#endif // ARGS_H
