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

    MPI_Request requests_odd_sends[MAX_REQUESTS];
    MPI_Request requests_odd_recvs[MAX_REQUESTS];
    MPI_Request requests_even_sends[MAX_REQUESTS];
    MPI_Request requests_even_recvs[MAX_REQUESTS];

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
        MPI_Startall(MAX_REQUESTS, requests_odd_sends);
        //recv_halo_A_start(&status);
        MPI_Startall(MAX_REQUESTS, requests_odd_recvs);
#endif

    for (int i=0;i<1000;++i) {
        //send_edge_B_wait(&status);
        MPI_Waitall(MAX_REQUESTS, requests_even_sends, MPI_STATUSES_IGNORE);
        //recv_halo_A_wait(&status);
        MPI_Waitall(MAX_REQUESTS, requests_odd_recvs, MPI_STATUSES_IGNORE);
        compute_edge_B();
        //send_edge_B_start(&status);
        MPI_Startall(MAX_REQUESTS, requests_even_sends);
        //recv_halo_A_start(&status);
        MPI_Startall(MAX_REQUESTS, requests_odd_recvs);

        compute_middle_B();

        //send_edge_A_wait(&status);
        MPI_Waitall(MAX_REQUESTS, requests_odd_sends, MPI_STATUSES_IGNORE);
        //recv_halo_B_wait(&status);
        MPI_Waitall(MAX_REQUESTS, requests_even_recvs, MPI_STATUSES_IGNORE);
        compute_edge_A();
        //send_edge_A_start(&status);
        MPI_Startall(MAX_REQUESTS, requests_odd_sends);
        //recv_halo_B_start(&status);
        MPI_Startall(MAX_REQUESTS, requests_even_recvs);

        compute_middle_A();
    }

    //send_edge_B_wait(&status);
    MPI_Waitall(MAX_REQUESTS, requests_even_sends, MPI_STATUSES_IGNORE);
    //recv_halo_A_wait(&status);
    MPI_Waitall(MAX_REQUESTS, requests_odd_recvs, MPI_STATUSES_IGNORE);
    //send_edge_A_wait(&status);
    MPI_Waitall(MAX_REQUESTS, requests_odd_sends, MPI_STATUSES_IGNORE);
    //recv_halo_B_wait(&status);
    MPI_Waitall(MAX_REQUESTS, requests_even_recvs, MPI_STATUSES_IGNORE);

    //comms_postloop(&status);
    for (int r=0;r<MAX_REQUESTS;++r) {
        MPI_Request_free(&requests_odd_sends[r]);
        MPI_Request_free(&requests_odd_recvs[r]);
        MPI_Request_free(&requests_even_sends[r]);
        MPI_Request_free(&requests_even_recvs[r]);
    }

    /* write data to storage from A */

    return 0;
}
