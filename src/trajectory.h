#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <stdio.h>
#include "common.h"
#include <stdint.h>
#include "models.h"

void calculate_position_and_speed(double h, double d, Vector geo_field, Vector *car_vel, Vector *geo_pos);
void trajectory_simulation_cutoff(InfileHeader *header, const Model **models, ModelsParams *models_params, FILE *outfile, uint64_t step_limit, double tu_angle, double particle_mass, uint32_t parallel);
void trajectory_simulation(InfileHeader *header, const Model **models, ModelsParams *models_params, FILE *tragsm, FILE *trasph, uint64_t step_limit, double tu_angle, double particle_mass);

#endif // TRAJECTORY_H