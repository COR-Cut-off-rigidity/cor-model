#include "field.h"
#include "utils.h"
#include <math.h>
#include <omp.h>

void field_simulation(InfileHeader *header, ModelsParams *params, const Model **models, FILE *outfile, uint32_t parallel) {
    double del_the = header->the1 - header->the0;
    double del_phi = header->fi1 - header->fi0;

    int n_the = 0, n_phi = 0;
    
    if(header->lat_step > 0.) {
        n_the = (int) (del_the / header->lat_step);
    }

    if(header->lon_step > 0.) {
        n_phi = (int) (del_phi / header->lon_step);
    }

    Vector sph_pos = {
        .r = header->r0 * DEF_RE,
        .theta = .0, 
        .phi = .0
    };

    for(int j = 0; j <= n_the; j++) {
        double the0a = header->the0 + (double)(j*header->lat_step);

        #pragma omp parallel for ordered schedule(static, 1) if (parallel)
        for(int i = 0; i < n_phi; i++) {
            double phi0a = header->fi0 + (double)(i*header->lon_step);
            
            Vector new_sph_pos = {
                .r = sph_pos.r,
                .theta = (90. - the0a) * DEF_TO_RAD,
                .phi = phi0a * DEF_TO_RAD
            };

            Vector pos;
            sph_to_car(&new_sph_pos, &pos);

            Vector gsm_pos;
            geo_to_gsm(&pos, &gsm_pos, params->A1, params->A2, params->A3);

            Vector geo_field;
            calculate_field(gsm_pos, &geo_field, models, params);

            double b = sqrt(geo_field.x*geo_field.x + geo_field.y*geo_field.y + geo_field.z*geo_field.z);

            #pragma omp ordered
            fprintf(outfile, "%9.2f%9.2f%11.2f%11.2f%11.2f%11.2f\n", the0a, phi0a, b, geo_field.x, geo_field.y, geo_field.z);
        }
    }
}
