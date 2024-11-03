#include <string.h>
#include "io.h"

int read_cutoff_infile_header(FILE *input_file, InfileHeader *infile_header) {
    int items = fscanf(
        input_file,
        "%lf %lf %lf\n%lf %lf %lf\n%lf %lf\n%d %d %d %d %d %d %d\n%d %d %d %lf",
        &infile_header->rig, &infile_header->zn, &infile_header->rk, &infile_header->r0,
        &infile_header->the0, &infile_header->fi0, &infile_header->the1, &infile_header->fi1,
        &infile_header->iy, &infile_header->mes, &infile_header->ide, &infile_header->id,
        &infile_header->ih, &infile_header->min, &infile_header->is, &infile_header->nk1,
        &infile_header->iopt, &infile_header->ist, &infile_header->del
    );

    if(items < 0) {
        return -1;
    }

    return 0;
}

int write_cutoff_outfile_header(FILE *output_file, InfileHeader *infile_header, const Model **models, uint64_t steps_count) {
    if(fprintf(output_file, "\n\n\n                ASYMPTOTIC COORDINATES\n   calculated by model(s): ") < 0) {
        return -1;
    }

    const char *pad = "";
    for(int i = 0; models[i]; i++) {
        if(fprintf(output_file, "%s%s", pad, models[i]->name) < 0) {
            return -1;
        }
        pad = ", ";
    }

    int res = fprintf(output_file, "\n Station with geo.latitude:%9.3f  & longitude:%9.3f & radius : %8.5f\n"
        " Direction of trajectory with latitude: %9.3f & longitude: %9.3f\n"
        " Date: %4d %2d %2d  time: %2d hod %2d min %2d sec\n"
        " Starting rigidity : %8.4f GV Epsilon=%7.4f\n Limit of total number of steps : %lu \n\n"
        " rig : v : rad : eth : efi : ath : afi : time : length\n",
        infile_header->the0, infile_header->fi0, infile_header->r0, infile_header->the1,
        infile_header->fi1, infile_header->iy, infile_header->mes, infile_header->ide,
        infile_header->ih, infile_header->min, infile_header->is, infile_header->rig,
        infile_header->del, steps_count);

    if(res < 0) {
        return -1;
    }

    return 0;
}

int read_field_infile_header(FILE *input_file, InfileHeader *infile_header) {
    int items = fscanf(
        input_file,
        "%lf %lf %lf\n%lf %lf %lf %lf\n%d %d %d %d %d %d %d\n",
        &infile_header->r0, &infile_header->the0, &infile_header->fi0, &infile_header->lat_step,
        &infile_header->lon_step, &infile_header->the1, &infile_header->fi1, &infile_header->iy,
        &infile_header->mes, &infile_header->ide, &infile_header->id, &infile_header->ih,
        &infile_header->min, &infile_header->is
    );

    if(items < 0) {
        return -1;
    }

    return 0;
}

int write_field_outfile_header(FILE *output_file, InfileHeader *infile_header, const Model **models) {
    if(fprintf(output_file, "\n\n\n                GEOMAGNETIC FIELD\n Geomagnetic Model(s): ") < 0) {
        return -1;
    }

    const char *pad = "";
    for(int i = 0; models[i]; i++) {
        if(fprintf(output_file, "%s%s", pad, models[i]->name) < 0) {
            return -1;
        }
        pad = ", ";
    }

    int items = fprintf(output_file, 
        "\n lat step: %.3f, lon step: %.3f\n"
        " Date: %4d %2d %2d  time: %2d h %2d min %2d sec\n\n"
        " lat : lon : B [nT] : BXGEO [nT] : BYGEO [nT] : BZGEO [nT]\n",
        infile_header->lat_step, infile_header->lon_step,
        infile_header->iy, infile_header->mes, infile_header->ide,
        infile_header->ih, infile_header->min, infile_header->is
    );
    
    if(items < 0) {
        return -1;
    }

    return 0;
}

int read_ts05_params(FILE *input_file, ModelsParams *params) {
    int items = fscanf(
        input_file,
        "%lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n",
        &params->PARMOD[1], &params->PARMOD[0], &params->PARMOD[2],
        &params->PARMOD[3], &params->PARMOD[4], &params->PARMOD[5],
        &params->PARMOD[6], &params->PARMOD[7], &params->PARMOD[8],
        &params->PARMOD[9]
    );

    if(items != 10) {
        return -1;
    }

    return 0;
}

int read_t96_params(FILE *input_file, ModelsParams *params) {
    int items = fscanf(
        input_file,
        "%lf %lf %lf %lf\n",
        &params->PARMOD[1], &params->PARMOD[0], &params->PARMOD[2], &params->PARMOD[3]
    );

    if(items != 4) {
        return -1;
    }

    return 0;
}
