#include <unistd.h>
#include "compute.h"

#define REF_TIME 10

void compute_edge_A(int tea) {
    usleep(tea);
}

void compute_edge_B(int teb) {
    usleep(teb);
}

void compute_middle_A(int tma) {
  usleep(tma);
}

void compute_middle_B(int tmb) {
  usleep(tmb);
}

