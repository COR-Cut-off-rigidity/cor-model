#include <argp.h>
#include <string.h>
#include "args.h"

#ifndef DEF_COR_ENGINE_AUTHOR_EMAIL
#define DEF_COR_ENGINE_AUTHOR_EMAIL NULL
#endif

#define ARGP_KEY_EPSILON_ANGLE 500
#define ARGP_KEY_STEPS_LIMIT 501
#define ARGP_KEY_NK_LIMIT 502

const char *argp_program_version = "cor-engine " DEF_COR_ENGINE_VERSION;
const char *argp_program_bug_address = DEF_COR_ENGINE_AUTHOR_EMAIL;

static const char doc[] = "Trajectory simulation of cosmic ray particle in Earth's magnetosphere";
static const char args_doc[] = "INFILE";

static struct argp_option options[] = {
    { "par", 'p', 0, 0, "Parallel computation" },
    { "output", 'o', "FILE", 0, "Output to FILE instead of standard output. This option is required for trajectory calculation mode." },
    { "steps-limit", ARGP_KEY_STEPS_LIMIT, "LIMIT", 0, "Single trajectory steps limit (default: 25000)" },
    { "epsilon", ARGP_KEY_EPSILON_ANGLE, "RAD", 0, "Epsilon angle, smaller values yield higher trajectory precision (default: 0.01 rad)" },
    { "nk-limit", ARGP_KEY_NK_LIMIT, "LIMIT", 0, "Circle division limit; larger values yield higher trajectory precision (default: 50000)" },
    { "mode", 'm', "MODE", 0, "Calculation modes: cutoff, field, trajectory (default: cutoff)"},
    { "particle-type", 't', "TYPE", 0, "Particle type: electron or proton (default: proton)" },
    { 0, 0, 0, 0, "Geomagnetic field models:" },
    { "internal-field", 'i', "MODEL", 0, "Use internal field model. Supported internal field models: igrf9-14, chaos7.16-7.18, chaos8.1-8.3, hist, cals10k.2 (default: igrf14)" },
    { "external-field", 'e', "MODEL", 0, "Use external field model. Supported external field models: ts05, t96 (default: none)" },
    { "crustal-field", 'c', "MODEL", 0, "Use crustal field model. Supported crustal field models: chaos7.16-7.18, chaos8.1-8.3 (default: none)"},
    { 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    CLIArguments* arguments = state->input;

    switch (key) {
        case 'p':
            arguments->options |= OPT_PARALLEL;
            break;
        case 'o':
            arguments->output_file_name = arg;
            break;
        case 'm':
            arguments->mode = arg;
            break;
        case 't':
            arguments->particle_type = arg;
            break;
        case 'i':
            arguments->internal_model_name = arg;
            break;
        case 'e':
            arguments->external_model_name = arg;
            break;
        case 'c':
            arguments->crustal_model_name = arg;
            break;
        case ARGP_KEY_STEPS_LIMIT:
            arguments->steps_limit = arg;
            break;
        case ARGP_KEY_EPSILON_ANGLE:
            arguments->epsilon_angle = arg;
            break;
        case ARGP_KEY_NK_LIMIT:
            arguments->nk_limit = arg;
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= 1) {
                argp_usage(state);
            }

            arguments->input_file_name = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 1) {
                argp_usage(state);
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

void args_parse(int argc, char **argv, CLIArguments *arguments) {
    struct argp argp = { options, parse_opt, args_doc, doc };
    argp_parse(&argp, argc, argv, 0, NULL, arguments);
}
