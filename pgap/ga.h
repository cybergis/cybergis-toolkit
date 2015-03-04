/******************************************************************************
 * ga.h: header file for the main program                                     *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef GA_H
#define GA_H

#include <stdio.h>
#ifdef PGAMODE
#include "mpi.h"
#include "pga.h"
#ifdef PGA_THREAD_MODE
#include <pthread.h>
#endif
#endif
#include "data.h"
#include "log.h"
#include "addr.h"

#include "myrng.h"

#define EV_COUNT  3
#define EV_OBJV   0
#define EV_FITV   1
#define EV_UFITV  2 
typedef struct {
	int * solution; // a vector of size of items
	// measure set: obj, fitness, and unfitness
	long long ev[EV_COUNT]; // obj func value
} Chrom;

// population statistics
typedef struct {
	int total_improve;
#ifdef PGAMODE
	int mig_improve;
#endif
	double startT;
	double bestT;
	double endT;
	long long min_fitv;
	long long max_fitv;
	long long avg_fitv;
	double commT;
	double comm_sendT;
	double comm_recvT;
#ifdef T_PROFILING
	double ioT;
	double immigrateT;
	double selectionT;
	double injectT;
	double crossmutT; // crossmut or direct copy mig
	double feasibT;
	double improveT;
	double evalT;
	double replaceT;
#endif
} STATTYPE;

// GA routines 
// general routines 
void gen_init_pop(int seeding);
void gen_init_chrom_constraint(Chrom * chrom);
void improve_feasibility_chu(Chrom * chrom);
void improve_quality_chu(Chrom * chrom);
void get_cap(Chrom * chrom, long long *cap);
int is_duplicate(int size, Chrom chrom, int inclusive, int checkAll);
void eval_chrom(Chrom * chrom);
void solution_verify(Chrom child);
int eval_pop(void);
void eval_pop_update(Chrom newc, Chrom oldc);
void rank_pop(int);
void quicksort (int lo, int hi, int which);
void quicksort_verify(void);
void swap(int i, int j, int which);
void insertsort (int which);
int get_wheel_selection(void);
void selection(int *p1, int *p2);
void rank_tournament(int *p1, int *p2, int which);
void binary_tournament(int *p1, int *p2);
void crossover(Chrom parent1, Chrom parent2, Chrom *child);
void mutate(Chrom * child);
int replacewith(Chrom child);
void roulette_init(void);
void chrom_copy(Chrom c1, Chrom *c2);
// utils
double get_ga_time(void);
void print_pop(int howmany);
void print_chrom(Chrom chrom, char *title);
void print_stat(FILE * myout);
void print_rank(int ev_i);
void print_config(void);
void output_result(int * solution, long long objv, double exec_time);
void * search(void *args);// thread-type def
#ifdef PGAMODE
int emigrate(int * emi_buffer, int newimprove, int noimprove);
int immigrate(Chrom * imi, int *origin);
//void print_imi(void);
#endif
// read config
void config(int argc, char **argv);
void gastat_init();

// data.c
extern VTYPE *v;
extern VTYPE *vt;
extern WTYPE *w;
extern WTYPE *wt;
extern WTYPE *b;
extern int n;
extern int m;

// log.c
extern IDTYPE glogId, globalId;
extern char glogDir[];
extern char glogPrefix[];
extern int debug;

#ifdef PGAMODE
#ifdef PGA_THREAD_MODE
extern pthread_mutex_t emi_mutex, imi_mutex;
extern int emi_ready, imi_ready;
#endif
extern int emi_size, imi_size, **imi_buffer, migrate_freq, immigrate_freq;
extern int myrank, net_size;
extern MPI_Comm topoComm, migComm;
extern int send_seq; 
extern int snd_parallelism; 
// parallel version debug
extern int pdebug;
// performance study
extern long long imi_hist_chrom[IMI_HIST_BUFFER_SIZE]; // chrom buffer
extern int imi_hist_origin[IMI_HIST_BUFFER_SIZE]; // origin
extern int imi_hist_index; // a loop-around index
#endif
#endif
