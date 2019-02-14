#include <mpi.h>

#define MAX_REQUESTS 27

typedef struct comms_status_t {
    //pattern_t pattern;
    MPI_Request requests_odd_sends[MAX_REQUESTS];
    MPI_Request requests_odd_recvs[MAX_REQUESTS];
    MPI_Request requests_even_sends[MAX_REQUESTS];
    MPI_Request requests_even_recvs[MAX_REQUESTS];
} comms_status_t;

void comms_preloop(comms_status_t* status);
void comms_postloop(comms_status_t* status);

void send_edge_A_start(comms_status_t* status);
void recv_halo_A_start(comms_status_t* status);

void send_edge_A_wait(comms_status_t* status);
void recv_halo_A_wait(comms_status_t* status);

void send_edge_B_start(comms_status_t* status);
void recv_halo_B_start(comms_status_t* status);

void send_edge_B_wait(comms_status_t* status);
void recv_halo_B_wait(comms_status_t* status);

