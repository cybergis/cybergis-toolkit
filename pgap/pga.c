/******************************************************************************
 * pga.h: PGA                                                                 *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef PGA_C
#define PGA_C

#ifdef PGA_THREAD_MODE
#include <pthread.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "pga.h"
#include "ga.h"

int migrate_freq=50; // export interval as number of iter b/w two exports
int immigrate_freq=25; // import interval as number of iter b/w two imports

// message format: chrom_type, origin, send_seq, solution array
int mig_msglen = 0; // a msg's lenth as number of integers, set in psearch() 

// emigration buffer. size: emi_size * n
int emi_buffer_size=0; // in bytes
int emi_size=2; // number of solutions to export. command line param
// snd_parallelism: we count the sending of one emigrated data to all neighbors
// as one group send. snd_parallelism defines the max number of group sends
// to be performed before sync (MPI_Waitall()). This param determines how much
// overlapping of migration and GA computation is. Min value is 1. Bigger value
// requires more memory, but also increase the overlapping of comm/comp. Max
// value must be set to not cause MPI outgoing message buffer overflow. 
int snd_parallelism = 5; // default value. command line parameter
// emi_queue: collection of emi buffers for non-blocking sending. Function
// emigrate() fills one emi buffer; this buffer is handed to Ibsend; before 
// MPI_Waitall() is called, nobody should write this buffer again, based on MPI
// standard 3.0. emi_queue maintains a cyclic queue of emi buffers so that no
// Ibsend and emigrate() will access the same buffer. Since in one emi, sends
// to all neighbors share the same emi buffer, the size of emi_queue is then
// (snd_parallelism * emi_buffer_size)
int **emi_queue = NULL;
int emi_queue_index = 0;
int emi_queue_size = 0; // in bytes, must be multiple of emi_buffer_size
// remember all Ibsend requests
MPI_Request *sndreq_list = NULL;
// just for MPI_Waitall()
MPI_Status *sndreq_status_list = NULL;
int sndreq_size = 0; // should be (neighbor_count * snd_parallelism)
int sndreq_index = 0; // current avail monitoring holder space
#ifdef PGA_THREAD_MODE
int emi_ready = 0; // control b/w GA and comm thread. GA sets, comm reads 
pthread_mutex_t emi_mutex; // ensure m-e on shared emi data
#endif
// immigrant pool
// the capacity of pool: multiplier of neighbor_count
#define IMI_BUFFER_CAP 8
int imi_size=10; // the size of immigrant pool
int **imi_buffer=NULL; // immigrant pool: a circular queue
int *imi_chrom_index = NULL; 
int *imi_temp;
int imi_buffer_head = 0; // queue head

#ifdef PGA_THREAD_MODE
int imi_ready = 0; // control b/w GA and comm thread. comm reads, GA sets
pthread_mutex_t imi_mutex; // ensure m-e on shared imi data
#endif

// network info
int net_size = 0;
int island_topo=0; // topology, e.g., grid
int conn_degree=4; // connectivity: number of neighbors
int * neighbor=NULL; // neighbor id storage
int neighbor_count = 4;

// communication
MPI_Comm topoComm; // the grid topology MPI communicator
int myrank, mpi_rcode;
// Capacity of user-provided buffer for MPI non-blocking buffered sending:
// Because we know the msg size (from problem size), we can guide MPI to use 
// sufficient amount of buffer.
// Note: this is not the send buf parameter in Ibsend() function interface.
// This buffer is supplied for MPI to use.
// Guideline: 
// - Set to be the multiplier of (chrom_size * neighbor_count * emi_size)
// - set to be sufficiently large to facilitate hetero comm scenarios
// - set to be small to test Ibsend() and non-blocking receive
//int mpi_user_buffer_cap = 5; // how many emi data to hold. will be adjusted

// communication stat
int send_round = 0;
int recv_round = 0;
int send_seq = 0; // sequence number to record number of emi 
// not used: unless we solve the dealock problem b/w communicators
MPI_Comm migComm; // MPI communicator that includes itself and neighbors

int pdebug = 0;

// communication config
int comm_sync = 0;
int comm_blocking = 0;

// performance study: track immigrants' history
long long imi_hist_chrom[IMI_HIST_BUFFER_SIZE]; // chrom buffer
int imi_hist_origin[IMI_HIST_BUFFER_SIZE]; // chrom origin
int imi_hist_index = 0; // a loop-around index

// check to see an export is possible. if yes, return emi_buffer
// it's not if send queue is full and not complete.
// TODO: a receiver might not get turn to receive if multi processes run on
// a single node. In this case, we skip. But might be better if we can use
// MPI_Testsome() to allow delayed receive from some procs. But that will
// complicate the sending queue policy.
int * send_req() {
	if (sndreq_index < sndreq_size) return emi_queue[emi_queue_index]; // queue has space
	int flag;
	// we shall NEVER use MPI_Waitall(). In isend, it can cause deadlock.
	mpi_rcode = MPI_Testall(sndreq_size, sndreq_list, &flag, sndreq_status_list);
	if (flag) {
		sndreq_index = 0;
		emi_queue_index = 0;
		return emi_queue[emi_queue_index];
	}
	return NULL;
}

// send emigrants. must be paired with send_check()
int send_emi()
{
	int i;
	// try to fix MPI_ERR_BUFFER. 
	//int j;
	//for (j=1; j<neighbor_count; j++) {
	//	memcpy(emi_queue[emi_queue_index + j], emi_buffer, emi_buffer_size);
	//}
	// end try
	// send solutions to neighbors
	for (i=1; i<=neighbor_count; i++) {
#ifdef DEBUG_COMM
		fprintf(myout, "ss%d->%d: %d %d \n", myrank, neighbor[i], send_seq, sndreq_index);
		fflush(myout);
#endif 
#ifdef PGA_USE_IBSEND
		mpi_rcode = MPI_Ibsend(emi_queue[emi_queue_index], emi_size*mig_msglen, MPI_INT, neighbor[i], 999, topoComm, &(sndreq_list[sndreq_index]));
#else
		mpi_rcode = MPI_Isend(emi_queue[emi_queue_index], emi_size*mig_msglen, MPI_INT, neighbor[i], 999, topoComm, &(sndreq_list[sndreq_index]));
#endif
		//mpi_rcode = MPI_Ibsend(emi_queue[emi_queue_index + i - 1], emi_size*mig_msglen, MPI_INT, neighbor[i], 999, topoComm, &(sndreq_list[sndreq_index]));
		if (mpi_rcode != MPI_SUCCESS) {
			fprintf(stderr, ">>>%d: wrong MPI_Ibsend to %d, return=%d\n", myrank, neighbor[i], mpi_rcode);
			fflush(stderr);
		}
		sndreq_index = sndreq_index + 1;
	}
	emi_queue_index = emi_queue_index + 1;
	send_round ++;
/* using MPI_Waitall() is  deprecated 
	// throttle control: sync prev Ibsends if buffer use is full
	if (sndreq_index >= sndreq_size) {
		// wait for all to finish
		double t1, t2;
		t1 = get_ga_time();
		mpi_rcode = MPI_Waitall(sndreq_size, sndreq_list, sndreq_status_list);
		t2 = get_ga_time();
#ifdef DEBUG_COMM
		fprintf(myout, "w: %d %d %.6lf - %.6lf %d/%d\n", myrank, send_seq, t1, t2, sndreq_index, sndreq_size);
		for (i=0; i<sndreq_size; i++) {
			fprintf(myout, "%d ", sndreq_status_list[i].MPI_ERROR);
		}
		fprintf(myout, "\n");
		fflush(myout);
#endif 
		if (mpi_rcode != MPI_SUCCESS) {
			fprintf(stderr, ">>>%d: wrong MPI_Waitall, return=%d\n", myrank, mpi_rcode);
			fflush(stderr);
		}
		if (t2-t1 > 0.2) {
			fprintf(stdout, ">>>%d: check all sends used %.3lf seconds\n", myrank, t2-t1);
			fflush(stdout);
		}
		sndreq_index = 0;
		emi_queue_index = 0;
	} else {
		// next emi uses the next emi_buffer
		emi_queue_index ++;
		//emi_queue_index += neighbor_count ;
	}
	emi_buffer = emi_queue[emi_queue_index];
*/
#ifdef DEBUG_COMM
	fprintf(myout, "s: %d %d %d\n", myrank, send_seq, sndreq_index);
	fflush(myout);
#endif 
	return 1;
}

// receive immigrants
int recv_imi()
{
	int received = 0, countb, counts;
	//MPI_Comm_rank(topoComm, &myrank);
	MPI_Status status;
	int flag;
	mpi_rcode = MPI_Iprobe(MPI_ANY_SOURCE, 999, topoComm, &flag, &status);  
	if (mpi_rcode != MPI_SUCCESS) {
		fprintf(stderr, ">>>%d: wrong MPI_Iprobe, return=%d\n", myrank, mpi_rcode);
		fflush(stderr);
	}
	MPI_Get_count(&status, MPI_INT, &countb);
	counts = countb / (emi_size * mig_msglen);
#ifdef DEBUG_COMM
	if (counts * emi_size * mig_msglen != countb) {
		fprintf(myout, "rc[%lld]: %d: %d int in, %d solns and ...", iterations, myrank, countb, counts);
	}
	fflush(myout);
#endif
	while (flag) {
		// recv
		// imi_temp is a global variable for message sharing b/w recv_imi() and fill_imi()
		mpi_rcode = MPI_Recv(imi_temp, emi_size*mig_msglen, MPI_INT, MPI_ANY_SOURCE, 999, topoComm, &status );
		if (mpi_rcode != MPI_SUCCESS) {
			fprintf(stderr, ">>>%d: wrong MPI_Recv, return=%d\n", myrank, mpi_rcode);
			fflush(stderr);
		}
#ifdef DEBUG_COMM
		// r: myrank sender sequence
		fprintf(myout, "rr: %d %d %d\n", myrank, *(imi_temp+1), *(imi_temp+2));
		fflush(myout);
#endif
		// fill into immigrant pool
		if (!fill_imi()) {
			fprintf(myout, ">>>%d: Err type=%d from=%d seq=%d\n", myrank, *imi_temp,  *(imi_temp+1), *(imi_temp+2));
			fflush(myout);
		}

		received ++;
		// check remaining messages
		MPI_Iprobe(MPI_ANY_SOURCE, 999, topoComm, &flag, &status);
		MPI_Get_count(&status, MPI_INT, &countb);
		counts = countb / (emi_size * mig_msglen);
#ifdef DEBUG_COMM
		if (counts * emi_size * mig_msglen != countb) {
			fprintf(myout, "rc[%lld]: %d: %d int in, %d solns and ...", iterations, myrank, countb, counts);
		}
		fflush(myout);
#endif
	}
	recv_round ++;
#ifdef DEBUG_COMM
	fprintf(myout, "r: %d round %d got %d solns\n", myrank, recv_round, received);
	fflush(myout);
#endif
	return received;
}
// copy received immigrant to immigrant pool
int fill_imi()
{
	if (*imi_temp != MIG_ELITE && *imi_temp != MIG_RANDOM) {
		return 0;
	}
	// it's going to override in a forced way
	if (imi_chrom_index[imi_buffer_head] >= 0)
		imi_buffer_head = (imi_buffer_head + 1) % imi_size;
	// copy to buffer
	memcpy(imi_buffer[imi_buffer_head], imi_temp, sizeof(int) * emi_size * mig_msglen);
	imi_chrom_index[imi_buffer_head] = 0;	
	// debug
	//if (debug) printq('+', *imi_temp);
	return 1;
}

// fetch immigrant to be processed by immigrate() 
// returns - the type of first chrom
//		 - copy immigrant to buffer, buffer size must be n
// input buffer is usually Chrom->solution
//TODO: in multithread version, this method is mutual-exclusive w/ fill_imi()
int fetch_imi(int * buffer, int *origin)
{
	int v, newhead;
	if (imi_chrom_index[imi_buffer_head] >= 0) {
		// copy to buffer
		int i;
		int *start = imi_buffer[imi_buffer_head]+imi_chrom_index[imi_buffer_head] * mig_msglen;
		v = start[0]; // v is chrom type, non-zero
		if (v==0) {
			fprintf(stderr, ">>>%d: fetch_imi() got wrong type chrom data\n", myrank);
			fflush(stderr);
		}
		*origin = start[1];
		// copy chrom data, skip header 
		for (i=3; i<mig_msglen; i++)
			buffer[i-3] = start[i];
		imi_chrom_index[imi_buffer_head] ++;
		if (imi_chrom_index[imi_buffer_head] >= emi_size) {
			imi_chrom_index[imi_buffer_head] = -1; // reset
			newhead = (imi_buffer_head==0)?imi_size-1:imi_buffer_head-1;
			if (imi_chrom_index[newhead] >= 0)
				imi_buffer_head = newhead;
		}
	} else {
		v = 0; *origin = -1;
	}
	return v;
}

// parallel search
void psearch()
{
	int i, j;

	// init: mpi. note: MPI_Init and MPI_Finalize are in main()
	int globalRank;
	MPI_Comm_rank( MPI_COMM_WORLD, &globalRank );
	MPI_Comm_size( MPI_COMM_WORLD, &net_size);
	// set msglen for migration communications: chrom_type, origin, sequence number, chrom data
	mig_msglen = 3 + n;
	// create grid topology
	int dims[2], periods[2];
	// find dim sizes closest to square root of np
	j = get_best_dim(net_size, &dims[0], &dims[1]); //dim0:Y; dim1:X
	if (globalRank == 0) {
		fprintf(stdout, "np: %d; topo: %d x %d. cost: %d loops\n", net_size, dims[0], dims[1], j);
		fflush(stdout);
	}
	// make sure left,right,up,below are diff. procs 
	periods[0] = (dims[0]>=3)?1:0;
	periods[1] = (dims[1]>=3)?1:0;
	//periods[0] = periods[1] = 1; // loop around
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &topoComm);
	MPI_Comm_rank(topoComm, &myrank);
	int degree0 = 4; // adjacent neighbor count
	int neighbor_probe[degree0];
	// grid topology neighborhood
	MPI_Cart_shift(topoComm, 0, 1, &neighbor_probe[0], &neighbor_probe[1]); //up&down 'cause direction 0 means Y
	MPI_Cart_shift(topoComm, 1, 1, &neighbor_probe[2], &neighbor_probe[3]); //left&right 'cause direction 1 means X
	// myself + adjacent neighbors + random remote neighbors
	int rn = 0; //TODO: set random remote neighbor info. 
	neighbor = (int *)malloc(sizeof(int) * (degree0+rn+1));
	for (i=0; i<degree0+rn+1; i++)
		neighbor[i] = -1;
	neighbor[0] = myrank; // myself
	neighbor_count = 0;
	for (i=0; i<degree0; i++) {
		if (neighbor_probe[i] >= 0 && neighbor_probe[i] != globalRank ) {
			if ( (i % 2 == 0) || ((i % 2 == 1) && (neighbor_probe[i] != neighbor_probe[i-1]))) {
				neighbor[neighbor_count + 1] = neighbor_probe[i];
				neighbor_count ++;
			}
		}
	}
	//TODO: add remote neighbor to the list 

/*
	// create migration communication world
	MPI_Group topoGroup, migGroup;
	MPI_Comm_group(topoComm, &topoGroup);
	MPI_Group_incl(topoGroup, neighbor_count+1, neighbor, &migGroup);
	MPI_Comm_create(topoComm, migGroup, &migComm); 
*/
	// set global id
	globalId = (IDTYPE) myrank;

	// init: export-related data structure
	// number of solutions to send in one export
	emi_size = 2;
	// snd_parallelism is set at command line or using default value
	emi_buffer_size = emi_size * mig_msglen * sizeof(int);
	emi_queue_size = snd_parallelism;
	//emi_queue_size = snd_parallelism * neighbor_count;
	emi_queue = (int **)malloc(emi_queue_size * sizeof(int *));
	emi_queue_index = 0;
	for (i=0; i<emi_queue_size; i++) {
		emi_queue[i] = (int *)malloc(emi_buffer_size);
		if (emi_queue[i] == NULL) {
			fprintf(myout, "ERROR: could not malloc %d bytes of memory for emi_queue[%d]\n", emi_buffer_size, i);
			fflush(myout);
			exit(1);
		}
		memset(emi_queue[i], 0, emi_buffer_size);
	}
	sndreq_size = neighbor_count * snd_parallelism;
	sndreq_list = (MPI_Request *)malloc(sizeof(MPI_Request) * sndreq_size);
	sndreq_status_list = (MPI_Status *)malloc(sizeof(MPI_Status) * sndreq_size);
	sndreq_index = 0;
#ifdef DEBUG_COMM
	if (myrank == 0) {
		fprintf(myout, "emi_queue size: %d buffers, set to hold %d exports in a raw, each buffer is %dx%dx%dx%lu bytes.\n", emi_queue_size, snd_parallelism, neighbor_count, emi_size, mig_msglen, sizeof(int));
		fflush(myout);
	}
#endif

	// init: import-related data structure
	imi_size = neighbor_count * IMI_BUFFER_CAP;
	imi_buffer = (int **)malloc(sizeof(int *) * imi_size);
	imi_buffer_head = 0;
	for (i=0; i<imi_size; i++) {
		imi_buffer[i] = (int *)malloc(sizeof(int) * emi_size * mig_msglen);
		memset(imi_buffer[i], 0, sizeof(int) * emi_size * mig_msglen);
	}
	imi_chrom_index = (int *)malloc(sizeof(int) * imi_size);
	for (i=0; i<imi_size; i++)
		imi_chrom_index[i] = -1; // -1 means not used
	imi_temp = (int *)malloc(sizeof(int) * emi_size * mig_msglen);
	memset(imi_temp, 0, sizeof(int) * emi_size * mig_msglen);
#ifdef PGA_USE_IBSEND
	// Buffered non-blocking send:
	// provide our own outgoing message buffer for MPI to use. This is
	// because we know problem size and the buffer size needed better.
	// NOTE: ibsend's robustness differs across diff MPI. NOT RECOMMENDED
	int mpi_buffer_size;
	MPI_Pack_size (snd_parallelism*neighbor_count*emi_size*mig_msglen, MPI_INT, topoComm, &mpi_buffer_size);
	mpi_buffer_size += (snd_parallelism*neighbor_count*MPI_BSEND_OVERHEAD);
// constant as the multiplier of basic outgoing message buffer size.
// useful when there are more MPI processes than number of cores on a node.
// don't know the right value yet. But it seems must be at least 2
#define MY_MPI_SNDBUF_FACTOR 100
	mpi_buffer_size *= MY_MPI_SNDBUF_FACTOR;
	int *mpi_buffer = (int *)malloc(mpi_buffer_size);
	if (mpi_buffer == NULL) {
		fprintf(myout, "ERROR: could not malloc %d bytes of memory\n", mpi_buffer_size);
		fflush(myout);
		exit(1);
	}
	MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);
#endif


	// performance study: init history buffers
	memset(imi_hist_chrom, 0, sizeof(long long) * IMI_HIST_BUFFER_SIZE);
	memset(imi_hist_origin, 0, sizeof(int) * IMI_HIST_BUFFER_SIZE);

	// print out pga configuration
	fprintf(myout, ">>>[pga_config]: topo=");
	for (i=0; i<neighbor_count+1; i++)
		fprintf(myout, "%d ", neighbor[i]);
	fprintf(myout, "\n");
	if (myrank == 0) {
		fprintf(myout, ">>>[pga_config]: neighbor_count=%d\n", neighbor_count);
		fprintf(myout, ">>>[pga_config]: emi_size=%d solutions\n", emi_size);
		fprintf(myout, ">>>[pga_config]: mig_msglen=%d solution variables\n", mig_msglen);
		fprintf(myout, ">>>[pga_config]: emi_buffer_size=%d bytes\n", emi_buffer_size);
		fprintf(myout, ">>>[pga_config]: snd_parallelism=%d export buffers\n", snd_parallelism);
		fprintf(myout, ">>>[pga_config]: emi_queue_size=%d bytes\n", emi_queue_size);
		fprintf(myout, ">>>[pga_config]: sndreq_size=%d Ibsends/Isends\n", sndreq_size);
#ifdef PGA_USE_IBSEND
		fprintf(myout, ">>>[pga_config]: MPI outgoing msg buf size=%d bytes (%dx%dx(%dx%dx%d+%d)x%d)\n", mpi_buffer_size, neighbor_count, snd_parallelism, emi_size, mig_msglen, (int)sizeof(int), MPI_BSEND_OVERHEAD, MY_MPI_SNDBUF_FACTOR);
		fprintf(myout, ">>>[pga_config]: MPI outgoing msg buf size multiplier=%d \n", MY_MPI_SNDBUF_FACTOR);
#endif
		fprintf(myout, ">>>[pga_config]: IMI_BUFFER_CAP=%d imports\n", IMI_BUFFER_CAP);
		fprintf(myout, ">>>[pga_config]: imi_size=%d single-neighbor imports\n", imi_size);
		fprintf(myout, ">>>[pga_config]: IMI_HIST_BUFFER_SIZE=%d solutions\n", IMI_HIST_BUFFER_SIZE);
	}
#ifdef HETERO
	fprintf(myout, ">>>[pga_config]: HETERO %d %d\n", myrank, HETERO);
#else
	fprintf(myout, ">>>[pga_config]: HETERO %d %d\n", myrank, 1);
#endif
	fflush(myout);

	// start of PGA
	MPI_Barrier(MPI_COMM_WORLD);

#ifdef PGA_THREAD_MODE
	// init: pthread mutex
	pthread_mutex_init(&emi_mutex, NULL);
	pthread_mutex_init(&imi_mutex, NULL);

	// start GA thread
	pthread_t ga_thread;
	pthread_create(&ga_thread, NULL, search, NULL);

	// go into comm loop
	while (1) {
		if (ga_done) break;
		// check if emi_buffer is ready. non-blocking
		if (pthread_mutex_trylock(&emi_mutex) == 0) {
			if (emi_ready) {
				for (i=0; i<conn_degree; i++) {
					//MPI_Ibsend(emi_buffer, ...);
				}
				emi_ready = 0;
			}
			pthread_mutex_unlock(&emi_mutex);
		}
		if (pthread_mutex_trylock(&imi_mutex) == 0) {
			if (!imi_ready) {
				for (i=0; i<conn_degree; i++) {
					//MPI_Irecv(imi_buffer[i], ...);
				}
				imi_ready = 1;
			}
			pthread_mutex_unlock(&imi_mutex);
		}
		//sched_yield(); // yield on single cpu arch
	}
	// destroy mutex
	pthread_mutex_destroy(&emi_mutex);
	pthread_mutex_destroy(&imi_mutex);
#else
	search(NULL);
#endif 

	// free resources
	free(neighbor);
	//free(emi_buffer); 
	for (i=0; i<emi_queue_size; i++) {
		free(emi_queue[i]);
	}
	free(emi_queue); 
	for (i=0; i<imi_size; i++) {
		free(imi_buffer[i]);
	}
	free(imi_buffer);
	free(imi_chrom_index);
#ifdef PGA_USE_IBSEND
	int actual_freed_size;
	MPI_Buffer_detach(&mpi_buffer, &actual_freed_size);
	free(mpi_buffer);
#endif
}
// find the rectangle size closest to the square root of np
int get_best_dim(int np, int *rowSize, int *colSize)
{
	int i = 0; // loop counter
	int p = (int) sqrt(np);
	int q = np / p;
	while (p > 1 &&(p * q != np || p > q)) {
		p --;
		q = np / p;
		i++;
	}
	*rowSize = p;
	*colSize = q;
	return i;
}
// search in immigrant history for a chrom's origin
int imi_find_origin(long long fitv)
{
	int i;
	for (i=0; i<IMI_HIST_BUFFER_SIZE; i++) {
		if (imi_hist_chrom[i]) {
			if (imi_hist_chrom[i] == fitv)
				return imi_hist_origin[i];
		} else break;
	}
	return -1;
}
#endif
