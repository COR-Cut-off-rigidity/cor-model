#ifndef FIELD_H
#define FIELD_H

#include "common.h"
#include "models.h"

void field_simulation(InfileHeader *header, ModelsParams *models_params, const Model **models, FILE *outfile, uint32_t parallel);

#endif // FIELD_H
