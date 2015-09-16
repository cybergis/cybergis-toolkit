#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char **argv )
 {
  //argv[1] is the size of the message to send
  //argv[2] is the number of messages to send
  if (argc != 3) {
    //printf("Usage: fence_pp message_size message_number\n");
    printf("Usage: pscw_pp message_size message_number\n");
    return -1;
  }
  int numprocs, myid, j, i;

  int msize = atoi(argv[1]);
  int mnum  = atoi(argv[2]);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Win win;

  if (myid == 0) printf("One-sided ping-pong "
    //"using fence for syncs\n");
    "using post/start-complete/wait for syncs\n");

  //We'll set up arrays of doubles for messaging
  double *send;
  double *recv;
  send = (double *) calloc(msize, sizeof (double));
  recv = (double *) calloc(msize, sizeof (double));

  //Initialize the send vector to some random numbers
  srand(time(NULL));
  for (j=0; j < msize; j++) {
    send[j] = myid + (double)rand()/(double)(RAND_MAX-1);
    recv[j] = 0;
  }

  //Identify the buffers to be used for getting -
  //Remember, both ranks are doing all of this...
  int soi = sizeof(int);
  MPI_Win_create(send, soi*msize, soi,
      MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  //Sneaky way to set the target rank
  int target_rank = 1 - myid;

  //Set up the group
  MPI_Group comm_group, group;
  MPI_Comm_group(MPI_COMM_WORLD, &comm_group);
  //For fence: set up a group with 0 and 1 in it
  //int ranks[2];
  //ranks[0] = 0;
  //ranks[1] = 1;
  //MPI_Group_incl(comm_group, 2, ranks, &group);
  //For PSCW: group includes ONLY the other process
  MPI_Group_incl(comm_group, 1, &target_rank, &group);

  double start, stop, diff;
  start = MPI_Wtime();
  printf("Rank %d RECV should start out as zero: %8.8f\n",
         myid, recv[0]);
  for (i=0; i < mnum; i++) {
    //MPI_Win_fence(0, win);
    //Start exposure epoch
    MPI_Win_post(group, 0, win);
    //Start access epoch
    MPI_Win_start(group, 0, win);
    MPI_Get(recv, msize, MPI_DOUBLE, target_rank,
               0, msize, MPI_DOUBLE, win);
    //Sync to make sure the get is complete
    MPI_Win_complete(win);
    MPI_Win_wait(win);
    //MPI_Win_fence(0, win);
  }
  stop = MPI_Wtime();
  diff = stop - start;
  printf("Rank %d RECV should no longer be zero: %8.8f\n",
         myid, recv[0]);
  printf("Rank %d, %d %d-double transactions took %8.8fs\n",
         myid, mnum, msize, diff);
  MPI_Group_free(&group);
  MPI_Win_free(&win);
  MPI_Finalize();
  return 0;
 }
