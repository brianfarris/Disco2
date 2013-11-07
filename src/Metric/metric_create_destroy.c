#define METRIC_PRIVATE_DEFS
#include <stdlib.h>
#include "../Headers/Sim.h"
#include "../Headers/Metric.h"

struct Metric* metric_create(struct Sim *theSim, double t, double r, double p, double z)
{
    struct Metric *g = malloc(sizeof(struct Metric));
    return g;
}

void metric_destroy(struct Metric *g)
{
    free(g);
}
