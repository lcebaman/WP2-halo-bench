#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

#include "compute.h"



#define KB 1024
#define MB 1024 * KB
#define GB 1024 * MB

//#define MAX_REQUESTS 27 //26?

// we want a request per neighbour, in a 2D grid we've got 4
#define MAX_REQUESTS 4



void print_usage() {
  printf("Usage:  [ap] -tedgeA num -tedgeB num -tmiddleA num"
         " -tmiddleB num -iter num -N num \n");
}

void print_parameters(int a, int b, int c) {
  printf("Iterations: %d\n",a);
  printf("Message size: %d bytes\n",b*sizeof(double));
  printf("Average samples: %d\n",c);
}

int main(int argc, char** argv) {

  int i,j;
  int N;
  int rank,size;
  int opt;
  int loopIter,avg;
  int tedgeA,tedgeB,tmiddleA,tmiddleB;

  MPI_Init(&argc, &argv);

  /* benchmark default values */
  loopIter = 1000;
  tedgeA = 100;
  tedgeB = 100;
  tmiddleA = 1000;
  tmiddleB = 1000;
  N = 120;
  avg = 10;

  opterr = 0;


  /*****************************************************/
  /*     Command line parameters                       */
  /*                                                   */
  /*****************************************************/
  /*

    Input Parameters:

    .  tedgeA   -  Time (useconds) for computing the edges of array A

    .  tedgeB   -  Time (useconds) for computing the edges of array B

    .  tmiddleA -  Time (useconds) for computing the middle of array A

    .  tmiddleB -  Time (useconds) for computing the middle of array B

    .  iter     -  Number of loop iterations

    .  N        -  Message size in halo (comm)

    .  avg      -  Number of iterations of to get average values

  */

  /* Specifying the expected options */
  static struct option long_options[] = {
    {"tedgeA",    required_argument, 0,  'a' },
    {"tedgeB",    required_argument, 0,  'b' },
    {"tmiddleA",  required_argument, 0,  'c' },
    {"tmiddleB",  required_argument, 0,  'd' },
    {"iter",      required_argument, 0,  't' },
    {"N",         required_argument, 0,  'N' },
    {"avg",       required_argument, 0,  'v' }
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
                                 long_options, &long_index )) != -1) {
    switch (opt) {
    case 'a' : tedgeA = atoi(optarg);
      break;
    case 'b' : tedgeB = atoi(optarg);
      break;
    case 'c' : tmiddleA = atoi(optarg);
      break;
    case 'd' : tmiddleB = atoi(optarg);
      break;
    case 't' : loopIter = atoi(optarg);
      break;
    case 'N' : N = atoi(optarg);
      break;
    case 'v' : avg = atoi(optarg);
      break;
    default:
      if(!rank)
        print_usage();
      exit(EXIT_FAILURE);
    }
  }

  /* if(rank ==0){ */
  /*   printf("Benchmark will run:\n"); */
  /*   printf("Iterations: %d\n",loopIter); */
  /*   printf("Message size: %d bytes: %d\n",N*sizeof(double)); */
  /*   printf("Average samples: %d\n",avg); */
  /*   fflush(stdout); */
  /* } */


  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  int dims[2] = {0,0}, periods[2] = {1,1}, coords[2];
  MPI_Comm cart;
  int up,down,left,right;

  /* Create cartesian communicator */
  MPI_Dims_create (size, 2, dims);

  MPI_Cart_create (MPI_COMM_WORLD, 2, dims, periods, 1, &cart);

#ifdef DEBUG
  if(!rank)
    printf("dims[0] = %d; dims[1] = %d\n",dims[0], dims[1]);
#endif

  /* Create partitioning of overall grid on processes */
  MPI_Cart_coords (cart, rank, 2, coords); /*My coordinate*/
  int my_cart_rank;
  /* use my coords to find my rank in cartesian group*/
  MPI_Cart_rank(cart, coords, &my_cart_rank);

  /* Store neighbors in the grid */
  MPI_Cart_shift (cart, 0, 1, &left, &right);
  MPI_Cart_shift (cart, 1, 1, &up, &down);


#ifdef DEBUG
  printf("rank: %d-> U=%d; D=%d; L=%d; R=%d\n",rank, up, down, left, right);
  fflush(stdout);
#endif

  // put neighbours in the same buffer
  int neighbours[MAX_REQUESTS];


  /* decompose the domain */

  int bx = N; // block size in x
  int by = N; // block size in y


  /* allocate space for buffers */
  double *bufA = (double*) calloc(1,(bx+2)*(by+2)*sizeof(double));
  double *bufB = (double*) calloc(1,(bx+2)*(by+2)*sizeof(double));


  MPI_Datatype dtA_edges[MAX_REQUESTS];
  MPI_Datatype dtA_halos[MAX_REQUESTS];
  MPI_Datatype dtB_edges[MAX_REQUESTS];
  MPI_Datatype dtB_halos[MAX_REQUESTS];

  /* create up-down datatype */
  MPI_Datatype up_down_dt;
  MPI_Type_contiguous(bx, MPI_DOUBLE, &up_down_dt);
  MPI_Type_commit(&up_down_dt);
  /* create right-left datatype */
  MPI_Datatype right_left_dt;
  MPI_Type_vector(by , 1, bx+2, MPI_DOUBLE, &right_left_dt);
  MPI_Type_commit(&right_left_dt);

  /* (0)up -> (1)down -> (2)right -> (3)left */
  neighbours[0] = up; neighbours[1] = down;
  neighbours[2] = right;  neighbours[3] = left;


  dtA_halos[0] = dtA_halos[1] = dtA_edges[0] =  dtA_edges[1] = up_down_dt;
  dtA_halos[2] = dtA_halos[3] = dtA_edges[2] =  dtA_edges[3] = right_left_dt;
  dtB_halos[0] = dtB_halos[1] = dtB_edges[0] =  dtB_edges[1] = up_down_dt;
  dtB_halos[2] = dtB_halos[3] = dtB_edges[2] =  dtB_edges[3] = right_left_dt;

  MPI_Request  requests[4 * MAX_REQUESTS];

#ifdef OVERLAP
  MPI_Request *requests_odd_recvs  = &requests[0 * MAX_REQUESTS];
  MPI_Request *requests_odd_sends  = &requests[1 * MAX_REQUESTS];
  MPI_Request *requests_even_recvs = &requests[2 * MAX_REQUESTS];
  MPI_Request *requests_even_sends = &requests[3 * MAX_REQUESTS];
#else
  MPI_Request *requests_even_recvs = &requests[0 * MAX_REQUESTS];
  MPI_Request *requests_odd_sends  = &requests[1 * MAX_REQUESTS];
  MPI_Request *requests_odd_recvs  = &requests[2 * MAX_REQUESTS];
  MPI_Request *requests_even_sends = &requests[3 * MAX_REQUESTS];
#endif

  MPI_Request *requests_for_edge_A = &requests[0 * MAX_REQUESTS];

#ifdef OVERLAP
  MPI_Request *requests_pre_loop   = &requests[0 * MAX_REQUESTS];
#else
  MPI_Request *requests_pre_loop   = &requests[1 * MAX_REQUESTS];
#endif

  MPI_Request *requests_for_edge_B = &requests[2 * MAX_REQUESTS];

  int tagA = 1, tagB = 2;

  if(!rank){
    printf("Benchmark run on %d MPI processes\n\n",size);
    print_parameters(loopIter, N, avg);
  }


  /* initial conditions */
  srand ( time ( NULL));
  for (i = 0; i < bx+2; i++)
    for (j = 0; j < by+2; j++)
      bufA[indx(i,j)] =  (double) rand();


  int sIndx[MAX_REQUESTS];
  int rIndx[MAX_REQUESTS];
  sIndx[0] = indx(1,1); sIndx[1] = indx(1,by); sIndx[2] = indx(bx,1);sIndx[3]= indx(1,1);
  rIndx[0] = indx(1,0); rIndx[1] = indx(1,by+1); rIndx[2] = indx(bx+1,1);rIndx[3]= indx(0,1);

  double s_preloop = MPI_Wtime();

  //comms_preloop(&status);
  for (int neighbour=0; neighbour < MAX_REQUESTS; ++neighbour) {

    /*MPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest,
      int tag, MPI_Comm comm, MPI_Request *request) */
    MPI_Send_init(&bufA[sIndx[neighbour]], 1, dtA_edges[neighbour], neighbours[neighbour],
                  tagA, comm, &requests_odd_sends[neighbour]);
    /* MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source,
       int tag, MPI_Comm comm, MPI_Request *request) */
    MPI_Recv_init(&bufA[rIndx[neighbour]], 1, dtA_halos[neighbour], neighbours[neighbour],
                  tagA, comm, &requests_odd_recvs[neighbour]);

    MPI_Send_init(&bufB[sIndx[neighbour]], 1, dtB_edges[neighbour], neighbours[neighbour],
                  tagB, comm, &requests_even_sends[neighbour]);
    MPI_Recv_init(&bufB[rIndx[neighbour]], 1, dtB_halos[neighbour], neighbours[neighbour],
                  tagB, comm, &requests_even_recvs[neighbour]);
  }


  double e_preloop = MPI_Wtime();

  double s_comp_mid_A=0.0;
  double s_comp_mid_B=0.0;

  double s_comp_edge_A=0.0;
  double s_comp_edge_B=0.0;
  double tmp=0.0;

  double s_mainloop, mainloop;
  double Tloop=0.0,Tcomp=0.0,Tcomm=0.0;
  double iTcomp, iTcomm;

  /* timer warm up */
  for(int i=0; i < 1000; i++)
    MPI_Wtime();



  /* Run the test several times to get better average */
  for(int a = 0; a < avg; a++){

    s_mainloop = MPI_Wtime();

#if INPUT_EXCLUDES_HALO
    //send_edge_A_start(&status);
    //recv_halo_A_start(&status);
    MPI_Startall(MAX_REQUESTS * 2, requests_pre_loop);
#endif


    for (int i=0; i < loopIter; ++i) {

      //send_edge_B_wait(&status);
      //recv_halo_A_wait(&status);
      MPI_Waitall(MAX_REQUESTS * 2, requests_for_edge_B, MPI_STATUSES_IGNORE);

#ifndef NOCOMPUTE
      tmp =  MPI_Wtime();
      //compute_edge_B(tedgeB);
      compute_edgeB(bufB, bx, by);
      s_comp_edge_B +=  MPI_Wtime() - tmp;
#endif

      MPI_Startall(MAX_REQUESTS * 2, requests_for_edge_B);
      //send_edge_B_start(&status);
      //recv_halo_A_start(&status);

#ifndef NOCOMPUTE
      tmp =  MPI_Wtime();
      //compute_middle_B(tmiddleB);
      compute_middleB(bufB, bx, by);
      s_comp_mid_B +=  MPI_Wtime() - tmp;
#endif

      //send_edge_A_wait(&status);
      //recv_halo_B_wait(&status);
      MPI_Waitall(MAX_REQUESTS * 2, requests_for_edge_A, MPI_STATUSES_IGNORE);

#ifndef NOCOMPUTE
      tmp =  MPI_Wtime();
      //compute_edge_A(tedgeA);
      compute_edgeA(bufA, bx, by);
      s_comp_edge_A +=  MPI_Wtime() - tmp;
#endif

      //send_edge_A_start(&status);
      //recv_halo_B_start(&status);
      MPI_Startall(MAX_REQUESTS * 2, requests_for_edge_A);

#ifndef NOCOMPUTE
      tmp =  MPI_Wtime();
      //compute_middle_A(tmiddleA);
      compute_middleA(bufA, bx, by);
      s_comp_mid_A +=  MPI_Wtime() - tmp;
#endif

    }

    //send_edge_B_wait(&status);
    //recv_halo_A_wait(&status);
    //send_edge_A_wait(&status);
    //recv_halo_B_wait(&status);
    MPI_Waitall(MAX_REQUESTS * 4, requests, MPI_STATUSES_IGNORE);

    /* Tloop = MPI_Wtime() - s_mainloop; */
    /* Tcomp = (s_comp_edge_B+s_comp_edge_A+s_comp_mid_B+s_comp_mid_A); */
    /* Tcomm = Tloop - Tcomp; */

    /* iteration timers */
    mainloop = MPI_Wtime() - s_mainloop;

#ifndef NOCOMPUTE
    iTcomp = (s_comp_edge_B+s_comp_edge_A+s_comp_mid_B+s_comp_mid_A);
#else
    iTcomp = 0.0;
#endif
    iTcomm = mainloop - iTcomp;

    /* Global timers */
    Tloop = Tloop + mainloop;
    Tcomp = Tcomp + iTcomp;
    Tcomm = Tcomm + iTcomm;
    /* reset timers */
    s_comp_edge_B = s_comp_edge_A = s_comp_mid_B = s_comp_mid_A = 0.0;

  }// end of avg loop

  /* Get average values */
  Tloop = Tloop/(double)avg;
  Tcomp = Tcomp/(double)avg;;
  Tcomm = Tcomm/(double)avg;

  if(!rank){
    printf("\n\n", Tloop);
    printf("Total Loop : %f\n", Tloop);
    printf("Total Computation : %f\n", Tcomp);
    printf("Total Communication : %f\n\n", Tcomm);
  }


  //comms_postloop(&status);
  for (int r=0;r<MAX_REQUESTS*4;++r) {
    MPI_Request_free(&requests[r]);
  }

  /* finalize MPI */
  MPI_Finalize();

  return 0;
}
