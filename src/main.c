#include <stdio.h>
#include <stdint.h>
#include "common.h"
#include <string.h>
#include "args.h"
#include "models.h"
#include "utils.h"
#include "io.h"
#include "trajectory.h"
#include "field.h"

int cutoff(const CLIArguments *args, const Model **models) {
    uint64_t steps_count = 25000;
    if(args->steps_count && string_to_uint64(args->steps_count, &steps_count) == -1) {
        fprintf(stderr, "Error: Invalid steps count '%s'\n", args->steps_count);
        return 1;
    }

    double tu_angle = .01;
    if(args->tu_angle && string_to_double(args->tu_angle, &tu_angle) == -1) {
        fprintf(stderr, "Error: Invalid TU angle '%s'\n", args->tu_angle);
        return 1;
    }

    double particle_mass = DEF_HM0_P;
    if(strcmp(args->particle_type, "proton") == 0) {
        particle_mass = DEF_HM0_P;
    } else if(strcmp(args->particle_type, "electron") == 0) {
        particle_mass = DEF_HM0_E;
    }

    FILE *input_file = fopen(args->input_file_name, "r");
    if(!input_file) {
        perror("Error: Cannot open input file for reading");
        return 1;
    }

    InfileHeader infile_header;
    if(read_cutoff_infile_header(input_file, &infile_header) != 0) {
        fprintf(stderr, "Error: Invalid input file\n");
        fclose(input_file);
        return 1;
    }

    if(infile_header.del <= 0.) {
        fprintf(stderr, "Error: Rigidity step size must be bigger than 0\n");
        fclose(input_file);
        return 1;
    }

    ModelsParams models_params;
    for(int i = 0; models[i]; i++) {
        if(!models[i]->setup_fun) {
            continue;
        }
        const void* setup_data = models[i]->setup_data;
        if(models[i]->setup_fun(input_file, &infile_header, &models_params, setup_data) != 0) {
            fprintf(stderr, "Error: Invalid input file\n");
            fclose(input_file);
            return 1;
        }
    }

    if(fclose(input_file) != 0) {
        perror("Error: Cannot close input file");
        return 1;
    }

    FILE *output_file = stdout;
    if(args->output_file_name) {
        output_file = fopen(args->output_file_name, "w");
        if(!output_file) {
            perror("Error: Cannot open output file for writing");
            return 1;
        }
    }

    if(write_cutoff_outfile_header(output_file, &infile_header, models, steps_count) != 0) {
        fprintf(stderr, "Error: Cannot write output file header\n");
        if(output_file != stdout) {
            fclose(output_file);
        }
        return 1;
    }

    trajectory_simulation_cutoff(&infile_header, models, &models_params, output_file, steps_count, tu_angle, particle_mass, args->options & OPT_PARALLEL);

    if(output_file != stdout && fclose(output_file) != 0) {
        perror("Error: Cannot close output file");
        return 1;
    }

    return 0;
}

int field(const CLIArguments *args, const Model **models) {
    FILE *input_file = fopen(args->input_file_name, "r");
    if(!input_file) {
        perror("Error: Cannot open input file for reading");
        return 1;
    }

    InfileHeader infile_header;
    if(read_field_infile_header(input_file, &infile_header) != 0) {
        fprintf(stderr, "Error: Invalid input file\n");
        if(input_file != stdin) {
            fclose(input_file);
        }
        return 1;
    }

    ModelsParams models_params;
    for(int i = 0; models[i]; i++) {
        if(!models[i]->setup_fun) {
            continue;
        }
        const void* setup_data = models[i]->setup_data;
        if(models[i]->setup_fun(input_file, &infile_header, &models_params, setup_data) != 0) {
            fprintf(stderr, "Error: Invalid input file\n");
            fclose(input_file);
            return 1;
        }
    }

    if(fclose(input_file) != 0) {
        perror("Error: Cannot close input file");
        return 1;
    }

    FILE *output_file = stdout;
    if(args->output_file_name) {
        output_file = fopen(args->output_file_name, "w");
        if(!output_file) {
            perror("Error: Cannot open output file for writing");
            return 1;
        }
    }

    if(write_field_outfile_header(output_file, &infile_header, models) != 0) {
        fprintf(stderr, "Error: Cannot write output file header\n");
        if(output_file != stdout) {
            fclose(output_file);
        }
        return 1;
    }

    field_simulation(&infile_header, &models_params, models, output_file, args->options & OPT_PARALLEL);

    if(output_file != stdout && fclose(output_file) != 0) {
        perror("Error: Cannot close output file");
        return 1;
    }

    return 0;
}

int trajectory(const CLIArguments *args, const Model **models) {
    if(!args->output_file_name) {
        fprintf(stderr, "Error: Output file name is required for trajectory calculation mode\n");
        return 1;
    }

    uint64_t steps_count = 25000;
    if(args->steps_count && string_to_uint64(args->steps_count, &steps_count) == -1) {
        fprintf(stderr, "Error: Invalid steps count '%s'\n", args->steps_count);
        return 1;
    }

    double tu_angle = .01;
    if(args->tu_angle && string_to_double(args->tu_angle, &tu_angle) == -1) {
        fprintf(stderr, "Error: Invalid TU angle '%s'\n", args->tu_angle);
        return 1;
    }

    double particle_mass = DEF_HM0_P;
    if(strcmp(args->particle_type, "proton") == 0) {
        particle_mass = DEF_HM0_P;
    } else if(strcmp(args->particle_type, "electron") == 0) {
        particle_mass = DEF_HM0_E;
    }

    FILE *input_file = fopen(args->input_file_name, "r");
    if(!input_file) {
        perror("Error: Cannot open input file for reading");
        return 1;
    }

    InfileHeader infile_header;
    if(read_cutoff_infile_header(input_file, &infile_header) != 0) {
        fprintf(stderr, "Error: Invalid input file\n");
        fclose(input_file);
        return 1;
    }

    if(infile_header.rig != infile_header.rk) {
        fprintf(stderr, "Error: Invalid input file - starting and ending rigidity must be the same for trajectory calculation mode\n");
        fclose(input_file);
        return 1;
    }

    ModelsParams models_params;
    for(int i = 0; models[i]; i++) {
        if(!models[i]->setup_fun) {
            continue;
        }
        const void* setup_data = models[i]->setup_data;
        if(models[i]->setup_fun(input_file, &infile_header, &models_params, setup_data) != 0) {
            fprintf(stderr, "Error: Invalid input file\n");
            fclose(input_file);
            return 1;
        }
    }

    if(fclose(input_file) != 0) {
        perror("Error: Cannot close input file");
        return 1;
    }

    size_t name_len = strlen(args->output_file_name);
    char tragsm_file_name[name_len + 8];
    strcpy(tragsm_file_name, "tragsm_");
    strcat(tragsm_file_name, args->output_file_name);

    char trasph_file_name[name_len + 8];
    strcpy(trasph_file_name, "trasph_");
    strcat(trasph_file_name, args->output_file_name);

    FILE *output_file_tragsm = fopen(tragsm_file_name, "w");
    if(!output_file_tragsm) {
        perror("Error: Cannot open tragsm output file for writing");
        return 1;
    }
    
    FILE *output_file_trasph = fopen(trasph_file_name, "w");
    if(!output_file_trasph) {
        fclose(output_file_tragsm);
        perror("Error: Cannot open trasfer output file for writing");
        return 1;
    }

    trajectory_simulation(&infile_header, models, &models_params, output_file_tragsm, output_file_trasph, steps_count, tu_angle, particle_mass);

    if(fclose(output_file_tragsm) != 0) {
        perror("Error: Cannot close tragsm output file");
        fclose(output_file_trasph);
        return 1;
    }

    if(fclose(output_file_trasph) != 0) {
        perror("Error: Cannot close tragsm output file");
        return 1;
    }

    return 0;
}

int main(int argc, char **argv) {
    CLIArguments args = { 0 };
    args.mode = "cutoff";
    args.particle_type = "proton";

    args_parse(argc, argv, &args);

    if(!args.internal_model_name && !args.external_model_name && !args.crustal_model_name) {
        args.internal_model_name = "igrf13";
    }

    int index = 0;
    const Model *models[4] = { NULL };
    
    if(args.internal_model_name) {
        const Model *model = get_model(internal_models, args.internal_model_name);
        if(!model) {
            fprintf(stderr, "Error: Invalid internal geomagnetic field model: '%s'\n", args.internal_model_name);
            return 1;
        }
        
        models[index++] = model;
    }

    if(args.external_model_name) {
        const Model *model = get_model(external_models, args.external_model_name);
        if(!model) {
            fprintf(stderr, "Error: Invalid external geomagnetic field model: '%s'\n", args.external_model_name);
            return 1;
        }

        models[index++] = model;
    }

    if(args.crustal_model_name) {
        const Model *model = get_model(crustal_models, args.crustal_model_name);
        if(!model) {
            fprintf(stderr, "Error: Invalid crustal geomagnetic field model: '%s'\n", args.crustal_model_name);
            return 1;
        }

        models[index++] = model;
    }

    if(strcmp(args.mode, "field") == 0) {
        return field(&args, models);
    }

    if(strcmp(args.mode, "cutoff") == 0) {
        return cutoff(&args, models);
    }

    if(strcmp(args.mode, "trajectory") == 0) {
        return trajectory(&args, models);
    }
    
    fprintf(stderr, "Error: Invalid calculation mode\n");
    return 1;
}
