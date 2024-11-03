#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdint.h>
#include "common.h"
#include "models.h"

int read_cutoff_infile_header(FILE *input_file, InfileHeader *infile_header);
int write_cutoff_outfile_header(FILE *output_file, InfileHeader *infile_header, const Model **models, uint64_t steps_count);
int read_field_infile_header(FILE *input_file, InfileHeader *infile_header);
int write_field_outfile_header(FILE *output_file, InfileHeader *infile_header, const Model **models);
int read_ts05_params(FILE *input_file, ModelsParams *params);
int read_t96_params(FILE *input_file, ModelsParams *params);

#endif // IO_H
