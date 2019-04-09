#include <unistd.h>
#include "compute.h"

#define REF_TIME 10

void compute_edge_A() {
    usleep(REF_TIME * REF_TIME);
}

void compute_edge_B() {
    usleep(REF_TIME * REF_TIME);
}

void compute_middle_A() {
    usleep(REF_TIME * REF_TIME * REF_TIME);
}

void compute_middle_B() {
    usleep(REF_TIME * REF_TIME * REF_TIME);
}

