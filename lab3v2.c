#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <math.h>

#define VECSIZE 65536
#define ITERATIONS 1

struct vecData {
  double val;
  int   rank;
}; 

double When() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

//pass these in with &
void reduceMaxArray(struct vecData *arr1, struct vecData *arr2) {
//printf("in rma\n");
  int i;
  for(i = 0; i < VECSIZE; i++) {
//    printf("val rank1: %f %d\n",arr1[i].val,arr1[i].rank);
//    printf("val rank2: %f %d\n",arr2[i].val,arr2[i].rank);

    if(arr1[i].val < arr2[i].val) {
      arr1[i] = arr2[i];
//      printf("swapped so value is %f %f\n",arr1[i].val,arr2[i].val);
    }
  } 
  return;
}

void bCast(int numDimensions, int rank, struct vecData *daVector) {
 
  int size = VECSIZE;
  int notParticipating = pow(2,(numDimensions-1))-1;
  int bitmask = pow(2,numDimensions-1);
  int curDimension;

  for(curDimension = 0; curDimension < numDimensions; curDimension++) {

    if( (rank & notParticipating) == 0) {

      if( (rank & bitmask) == 0) {
        int msg_dest = rank ^ bitmask;
        //printf("%db about to send to %d\n",rank,msg_dest);
        MPI_Send(daVector,VECSIZE,MPI_DOUBLE_INT,msg_dest,1,MPI_COMM_WORLD);
       // printf("%db sent so moving on\n",rank);
      } else {
        int msg_src = rank ^ bitmask;
        //printf("%db about to recv from %d\n",rank,msg_src);
        MPI_Recv(daVector,VECSIZE,MPI_DOUBLE_INT,msg_src,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //printf("%db recved so moving on\n",rank);
      }

    }
    notParticipating >>= 1;
    bitmask >>= 1;
  }
}

void reduce(int numDimensions,int rank,  struct vecData *daVector) {

  MPI_Status *status;

  int notParticipating = 0;
  int bitmask = 1;
  float newValue;
  struct vecData newVector[VECSIZE];

  int curDimension;
  for(curDimension = 0; curDimension < numDimensions; curDimension++) {
    if( (rank & notParticipating) == 0) {
      if( (rank & bitmask) != 0) {
/*
printf("%dr data to send: ",rank);
int i;
for(i=0;i<VECSIZE;i++){
printf("%f ",daVector[i].val);
}
printf("\n");
*/
        int msg_dest = rank ^ bitmask;
        //printf("%dr about to send to %d\n",rank,msg_dest);
        MPI_Send(daVector,VECSIZE,MPI_DOUBLE_INT,msg_dest,1,MPI_COMM_WORLD);
        //printf("%dr sent so moving on\n",rank);

      } else {
       int msg_src = rank ^ bitmask;
        //printf("%dr about to recv from %d\n",rank,msg_src);
        MPI_Recv(newVector,VECSIZE,MPI_DOUBLE_INT,msg_src,1,MPI_COMM_WORLD,status);
        reduceMaxArray(daVector,newVector);
        //printf("%dr recved so moving on\n",rank); 
/*
printf("%dr data recved: ",rank);
int i;
for(i=0;i<VECSIZE;i++){
printf("%f ",newVector[i].val);
}
printf("\n");
*/

      }
    }
    notParticipating = notParticipating ^ bitmask;
    bitmask <<= 1;
  }

}

main(int argc, char *argv[]) {
  int iproc, nproc,i, iter;
  char host[255];
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  gethostname(host,253);
  printf("I am proc %d of %d running on %s\n", iproc, nproc,host);

  // each process has an array of VECSIZE double: ain[VECSIZE]
  double ain[VECSIZE];
  int ind[VECSIZE];
  struct vecData in[VECSIZE];
  int myrank, root = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // Start time here
  srand(myrank+5);
  double start = When();

  int numDim = log10(nproc)/log10(2);

//do the broadcast and reduce ITERATIONS times
  for(iter = 0; iter < ITERATIONS; iter++) {

    //randomly initialize vectors
//    for(i = 0; i < VECSIZE; i++) {
//      ain[i] = rand();
//    printf("init proc %d [%d]=%f\n",myrank,i,ain[i]);
//    }
    for (i=0; i<VECSIZE; ++i) {
      in[i].val = rand() % 20;
      in[i].rank = myrank;
      //printf("%f ",in[i].val);
    }
    //printf("\n");

    reduce(numDim,myrank,in);

    // At this point, the answer resides on process root
//    if (myrank == root) {
      /* read ranks out*/
//      for (i=0; i<VECSIZE; ++i) {
//        printf("root out[%d] = %f from %d\n",i,in[i].val,in[i].rank);
//      }
//    }

    //now broadcast
    bCast(numDim,myrank,in);
//    for(i = 0; i < VECSIZE; i++) {
//      printf("final proc %d [%d]=%f from %d\n",myrank,i,in[i].val,in[i].rank);
//    }

  }
  MPI_Finalize();
  double end = When();
  if(myrank == root) {
    printf("Time %f\n",end-start);
  }
}

