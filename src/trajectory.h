#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <stdio.h>
#include "common.h"
#include <stdint.h>
#include "models.h"

typedef struct {
    uint64_t steps_limit;
    uint64_t nk_limit;
    double epsilon_angle;
    double particle_mass;
    uint32_t parallel;
} TrajectorySimulationParams;

void calculate_position_and_speed(double h, double d, Vector geo_field, Vector *car_vel, Vector *geo_pos);
void trajectory_simulation_cutoff(const TrajectorySimulationParams *params, const InfileHeader *header, const ModelsParams *model_params, const Model **models, FILE *outfile);
void trajectory_simulation(const TrajectorySimulationParams *params, const InfileHeader *header, const ModelsParams *model_params, const Model **models, FILE *tragsm, FILE *trasph);

#endif // TRAJECTORY_H