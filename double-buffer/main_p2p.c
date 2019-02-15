#include <mpi.h>

#define MAX_REQUESTS 27

void compute_edge_A() {};
void compute_edge_B() {};
void compute_middle_A() {};
void compute_middle_B() {};

int main() {

    /* read data from storage into A */
    double *bufA;
    double *bufB;


    MPI_Datatype dtA_edges[MAX_REQUESTS];
    MPI_Datatype dtA_halos[MAX_REQUESTS];
    MPI_Datatype dtB_edges[MAX_REQUESTS];
    MPI_Datatype dtB_halos[MAX_REQUESTS];

    MPI_Request  requests[4 * MAX_REQUESTS];

    MPI_Request *requests_even_recvs = &requests[0 * MAX_REQUESTS];
    MPI_Request *requests_odd_sends  = &requests[1 * MAX_REQUESTS];
    MPI_Request *requests_odd_recvs  = &requests[2 * MAX_REQUESTS];
    MPI_Request *requests_even_sends = &requests[3 * MAX_REQUESTS];

    MPI_Request *requests_for_edge_A = &requests[0 * MAX_REQUESTS];
    MPI_Request *requests_pre_loop   = &requests[1 * MAX_REQUESTS];
    MPI_Request *requests_for_edge_B = &requests[2 * MAX_REQUESTS];

    int tagA = 1, tagB = 2;
    MPI_Comm comm = MPI_COMM_WORLD;
    int neighbours[MAX_REQUESTS];

    //comms_preloop(&status);
    for (int neighbour=0;neighbour<MAX_REQUESTS;++neighbour) {
        MPI_Send_init(bufA, 1, dtA_edges[neighbour], neighbours[neighbour], tagA, comm, &requests_odd_sends[neighbour]);
        MPI_Recv_init(bufA, 1, dtA_halos[neighbour], neighbours[neighbour], tagA, comm, &requests_odd_recvs[neighbour]);
        MPI_Send_init(bufB, 1, dtB_edges[neighbour], neighbours[neighbour], tagB, comm, &requests_even_sends[neighbour]);
        MPI_Recv_init(bufB, 1, dtB_halos[neighbour], neighbours[neighbour], tagB, comm, &requests_even_recvs[neighbour]);
    }

#if INPUT_EXCLUDES_HALO
        //send_edge_A_start(&status);
        //recv_halo_A_start(&status);
        MPI_Startall(MAX_REQUESTS * 2, requests_pre_loop);
#endif

    for (int i=0;i<1000;++i) {
        //send_edge_B_wait(&status);
        //recv_halo_A_wait(&status);
        MPI_Waitall(MAX_REQUESTS * 2, requests_for_edge_B, MPI_STATUSES_IGNORE);
        compute_edge_B();
        MPI_Startall(MAX_REQUESTS * 2, requests_for_edge_B);
        //send_edge_B_start(&status);
        //recv_halo_A_start(&status);

        compute_middle_B();

        //send_edge_A_wait(&status);
        //recv_halo_B_wait(&status);
        MPI_Waitall(MAX_REQUESTS * 2, requests_for_edge_A, MPI_STATUSES_IGNORE);
        compute_edge_A();
        //send_edge_A_start(&status);
        //recv_halo_B_start(&status);
        MPI_Startall(MAX_REQUESTS * 2, requests_for_edge_A);

        compute_middle_A();
    }

    //send_edge_B_wait(&status);
    //recv_halo_A_wait(&status);
    //send_edge_A_wait(&status);
    //recv_halo_B_wait(&status);
    MPI_Waitall(MAX_REQUESTS * 4, requests, MPI_STATUSES_IGNORE);

    //comms_postloop(&status);
    for (int r=0;r<MAX_REQUESTS*4;++r) {
        MPI_Request_free(&requests[r]);
    }

    /* write data to storage from A */

    return 0;
}
