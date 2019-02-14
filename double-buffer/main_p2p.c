#include "comms_p2p.h"

int main() {

    comms_status_t status;
    comms_preloop(&status);

    /* read data from storage into A */

#if INPUT_EXCLUDES_HALO
        send_edge_A_start(&status);
        recv_halo_A_start(&status);
#endif

    for (int i=0;i<1000;++i) {
        send_edge_B_wait(&status);
        recv_halo_A_wait(&status);
        compute_edge_B();
        send_edge_B_start(&status);
        recv_halo_A_start(&status);

        compute_middle_B();

        send_edge_A_wait(&status);
        recv_halo_B_wait(&status);
        compute_edge_A();
        send_edge_A_start(&status);
        recv_halo_B_start(&status);

        compute_middle_A();
    }

    send_edge_B_wait(&status);
    recv_halo_A_wait(&status);
    send_edge_A_wait(&status);
    recv_halo_B_wait(&status);

    comms_postloop(&status);

    /* write data to storage from A */

    return 0;
}
