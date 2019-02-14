#include <mpi.h>
#include "comms_p2p.h"

void comms_preloop(comms_status_t* status) {}

void send_edge_A_start(comms_status_t* status) {
    MPI_Startall(MAX_REQUESTS, status->requests_odd_sends);
}

void recv_halo_A_start(comms_status_t* status) {
    MPI_Startall(MAX_REQUESTS, status->requests_odd_recvs);
}

void send_edge_A_wait(comms_status_t* status) {
    MPI_Waitall(MAX_REQUESTS, status->requests_odd_sends, MPI_STATUSES_IGNORE);
}

void recv_halo_A_wait(comms_status_t* status) {
    MPI_Waitall(MAX_REQUESTS, status->requests_odd_recvs, MPI_STATUSES_IGNORE);
}

void send_edge_B_start(comms_status_t* status) {
    MPI_Startall(MAX_REQUESTS, status->requests_even_sends);
}

void recv_halo_B_start(comms_status_t* status) {
    MPI_Startall(MAX_REQUESTS, status->requests_even_recvs);
}

void send_edge_B_wait(comms_status_t* status) {
    MPI_Waitall(MAX_REQUESTS, status->requests_even_sends, MPI_STATUSES_IGNORE);
}

void recv_halo_B_wait(comms_status_t* status) {
    MPI_Waitall(MAX_REQUESTS, status->requests_even_recvs, MPI_STATUSES_IGNORE);
}

void comms_postloop(comms_status_t* status) {
    for (int r=0;r<MAX_REQUESTS;++r) {
        MPI_Request_free(&status->requests_odd_sends[r]);
        MPI_Request_free(&status->requests_odd_recvs[r]);
        MPI_Request_free(&status->requests_even_sends[r]);
        MPI_Request_free(&status->requests_even_recvs[r]);
    }
}
