#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#include "compute.h"

#define indx(i,j) (j)*(bx+2)+(i)

//#define MAX_REQUESTS 27 //26?

// we want a request per neighbour, in a 2D grid we've got 4
#define MAX_REQUESTS 4

int main(int argc, char** argv) {
  
  int i,j;
  int N; 
  int rank,size;
  MPI_Init(&argc, &argv);
  

  if(argc < 2) N = 120;
  else N = atoi(argv[1]);
  
  
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  int dims[2] = {0,0}, periods[2] = {1,1}, coords[2];
  MPI_Comm cart;
  int north,south,east,west;
 
  
  // who are my neighbours?
  
  /* Create cartesian communicator */
  MPI_Dims_create (size, 2, dims);
  
  //printf("Procs in X %d and Y %d\n",dims[0], dims[1]);
  MPI_Cart_create (MPI_COMM_WORLD, 2, dims, periods, 0, &cart);
  
  /* Store neighbors in the grid */
  MPI_Cart_shift (cart, 0, 1, &west,  &east);
  MPI_Cart_shift (cart, 1, 1, &north, &south);
  
  /* Create partitioning of overall grid on processes */
  MPI_Cart_coords (cart, rank, 2, coords); /*My coordinate*/
  
#ifdef DEBUG
  printf("rank: %d-> N=%d; S=%d; E=%d; W=%d\n",north,south,east,west);
#endif

  // put neighbours in the same buffer
  int neighbours[MAX_REQUESTS];


  /* decompose the domain */
  
  int bx = N/dims[0]; // block size in x
  int by = N/dims[1]; // block size in y
  
#ifdef DEBUG
  printf("rank: %d-> N = %d , Xlocal = %d, Ylocal= %d and processes= %d\n", 
	 rank, N, bx, by, size);
#endif

  /* allocate space for buffers */
  double *bufA = (double*) calloc(1,(bx+2)*(by+2)*sizeof(double));
  double *bufB = (double*) calloc(1,(bx+2)*(by+2)*sizeof(double));
  
  
  MPI_Datatype dtA_edges[MAX_REQUESTS];
  MPI_Datatype dtA_halos[MAX_REQUESTS];
  MPI_Datatype dtB_edges[MAX_REQUESTS];
  MPI_Datatype dtB_halos[MAX_REQUESTS];

  /* create north-south datatype */
  MPI_Datatype north_south_dt;
  MPI_Type_contiguous(bx, MPI_DOUBLE, &north_south_dt);
  MPI_Type_commit(&north_south_dt);
  /* create east-west datatype */
  MPI_Datatype east_west_dt;
  MPI_Type_vector(by , 1, bx+2, MPI_DOUBLE, &east_west_dt);
  MPI_Type_commit(&east_west_dt);

  /* (0)north -> (1)south -> (2)east -> (3)west */
  neighbours[0] = north; neighbours[1] = south;
  neighbours[2] = east;  neighbours[3] = west; 


  dtA_halos[0] = dtA_halos[1] = dtA_edges[0] =  dtA_edges[1] = north_south_dt;
  dtA_halos[2] = dtA_halos[3] = dtA_edges[2] =  dtA_edges[3] = east_west_dt;
  dtB_halos[0] = dtB_halos[1] = dtB_edges[0] =  dtB_edges[1] = north_south_dt;
  dtB_halos[2] = dtB_halos[3] = dtB_edges[2] =  dtB_edges[3] = east_west_dt;

  MPI_Request  requests[4 * MAX_REQUESTS];

  MPI_Request *requests_even_recvs = &requests[0 * MAX_REQUESTS];
  MPI_Request *requests_odd_sends  = &requests[1 * MAX_REQUESTS];
  MPI_Request *requests_odd_recvs  = &requests[2 * MAX_REQUESTS];
  MPI_Request *requests_even_sends = &requests[3 * MAX_REQUESTS];

  MPI_Request *requests_for_edge_A = &requests[0 * MAX_REQUESTS];
  MPI_Request *requests_pre_loop   = &requests[1 * MAX_REQUESTS];
  MPI_Request *requests_for_edge_B = &requests[2 * MAX_REQUESTS];

  int tagA = 1, tagB = 2;
  
  
  /* initial conditions */
  srand ( time ( NULL));
  for (i = 0 ;i < N; i++)
    for (j = 0 ;j < N; j++)
      bufA[indx(i,j)] =  (double) rand();
  

  int sIndx[MAX_REQUESTS];
  int rIndx[MAX_REQUESTS];
  sIndx[0] = indx(1,1); sIndx[1] = indx(1,by); sIndx[2] = indx(bx,1);sIndx[3]= indx(1,1);
  rIndx[0] = indx(1,0); rIndx[1] = indx(1,by+1); rIndx[2] = indx(bx+1,1);rIndx[3]= indx(0,1);

  double s_preloop = MPI_Wtime();
  //comms_preloop(&status);
  for (int neighbour=0 ;neighbour < MAX_REQUESTS; ++neighbour) {
    
    /*MPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest,
      int tag, MPI_Comm comm, MPI_Request *request) */
    MPI_Send_init(&bufA[sIndx[neighbour]], 1, dtA_edges[neighbour], neighbours[neighbour], tagA,
                  comm, &requests_odd_sends[neighbour]);
    /* MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, 
                 int tag, MPI_Comm comm, MPI_Request *request) */
    MPI_Recv_init(&bufA[rIndx[neighbour]], 1, dtA_halos[neighbour], neighbours[neighbour], tagA,
                  comm, &requests_odd_recvs[neighbour]);
 
    MPI_Send_init(&bufB[sIndx[neighbour]], 1, dtB_edges[neighbour], neighbours[neighbour], tagB,
                  comm, &requests_even_sends[neighbour]);
    MPI_Recv_init(&bufB[rIndx[neighbour]], 1, dtB_halos[neighbour], neighbours[neighbour], tagB,
                  comm, &requests_even_recvs[neighbour]);
  }
  double e_preloop = MPI_Wtime();

  double s_comp_mid_A=0.0;
  double s_comp_mid_B=0.0;
  
  double s_comp_edge_A=0.0;
  double s_comp_edge_B=0.0;
  double tmp=0.0;
  
  double s_mainloop = MPI_Wtime();

#if INPUT_EXCLUDES_HALO
  //send_edge_A_start(&status);
  //recv_halo_A_start(&status);
  MPI_Startall(MAX_REQUESTS * 2, requests_pre_loop);
#endif
  
  for (int i=0;i < 1000; ++i) {

    //send_edge_B_wait(&status);
    //recv_halo_A_wait(&status);
    MPI_Waitall(MAX_REQUESTS * 2, requests_for_edge_B, MPI_STATUSES_IGNORE);
    
    tmp =  MPI_Wtime();
    compute_edge_B();
    s_comp_edge_B +=  MPI_Wtime() - tmp;
    
    MPI_Startall(MAX_REQUESTS * 2, requests_for_edge_B);
    //send_edge_B_start(&status);
    //recv_halo_A_start(&status);
    tmp =  MPI_Wtime();
    compute_middle_B();
    s_comp_mid_B +=  MPI_Wtime() - tmp;
    //send_edge_A_wait(&status);
    //recv_halo_B_wait(&status);
    MPI_Waitall(MAX_REQUESTS * 2, requests_for_edge_A, MPI_STATUSES_IGNORE);
    
    tmp =  MPI_Wtime();
    compute_edge_A();
    s_comp_edge_A +=  MPI_Wtime() - tmp;
    
    //send_edge_A_start(&status);
    //recv_halo_B_start(&status);
    MPI_Startall(MAX_REQUESTS * 2, requests_for_edge_A);
    
    tmp =  MPI_Wtime();
    compute_middle_A();
    s_comp_mid_A +=  MPI_Wtime() - tmp;
  }

  //send_edge_B_wait(&status);
  //recv_halo_A_wait(&status);
  //send_edge_A_wait(&status);
  //recv_halo_B_wait(&status);
  MPI_Waitall(MAX_REQUESTS * 4, requests, MPI_STATUSES_IGNORE);
  
  double e_mainloop = MPI_Wtime() - s_mainloop;
  

  if(!rank)
    printf("Ranks: %d; main_loop= %f; TotComp= %f; comp_edge_B= %f; comp_edge_A= %f; comp_mid_B= %f;"
	   "comp_mid_A= %f\n",
	   size,e_mainloop,(s_comp_edge_B+s_comp_edge_A+s_comp_mid_B+s_comp_mid_A),
	   s_comp_edge_B, s_comp_edge_A, s_comp_mid_B, s_comp_mid_A);
 

  //comms_postloop(&status);
  for (int r=0;r<MAX_REQUESTS*4;++r) {
    MPI_Request_free(&requests[r]);
  }

  /* write data to storage from A */
  if(!rank)
    printf("End of Benchmark\n");
  /* finalize MPI */
  MPI_Finalize();
  
  return 0;
}
