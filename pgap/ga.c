/******************************************************************************
 * ga.c: main program                                                         *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef GA_ASYNC_C
#define GA_ASYNC_C

//#define PGAMODE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include "ga.h"


// GA data structure
Chrom * population;
int * rank[EV_COUNT]; // population rank memory to sort chromosomes
char *ev_text[EV_COUNT]={"OBJV", "FITV", "UFITV"};
int pop_size = 1000; // population size 

// statistics
STATTYPE gastat;

// control parameters
// selection strategy
char strategy_selection = 'b'; // b-binary; r-rank
// replacement strategy
char strategy_replacement = 'u'; // u - unfittest, w - worst
// init pop improvement strategy
char strategy_init_pop = 'n'; // n-none; c-constraint; f-feasibility improve
// stopping rules 
char strategy_stop = 'i'; // fixed number of iterations; i-w/o improve
long long iterations; // loop control: number of iterations 
int max_iterations = 1000000; // maximum iterations 
int improve_threshold = 500; // upper-limit of non-improved iterations 
// walltime is not accurate: each process may reach time T at different time.
// this option is used to capture output before job is killed by cluster sched
int walltime = 3600; // upper-limit of exec time in seconds. 
// get other heuristic's result	
int seeding = 0; // seeding or not
// always keep the best solution
int elitism = 0;
float Pcross = 0.8; // probability of crossover, not used 
float Pmut = 0.2; // probability of mutation, not used
// rank-based selection wheel	
double *wheel; 

// stopping rule: solution quality
long long stopping_quality = 0;
// dataset 
char dataset[256];
// output file descriptor
FILE * myout;
char myoutputdir[256];
// random number generator
long global_seed = -1; // default
#ifdef GSL_SPRNG
gsl_rng * gsl_sprng_r = NULL;
#endif

#ifdef PGAMODE
// indicate main thread search() is done
int ga_done = 0;
#endif

// GA routines 
// general routines 
void gen_init_pop(int seeding)
{
	int i, j;
	FILE *f;
	// good-quality init chromosome from some approximate algorithm
	char * file1 = "./seed.solution";
	
	i = 0; 
	while (i<pop_size) {
		if (strategy_init_pop == 'n') {
			for (j=0; j<n; j++) {
				population[i].solution[j] = MYRANDI(m);
			}
		} else if (strategy_init_pop == 'c') {
			gen_init_chrom_constraint(population + i);
		}
		improve_feasibility_chu(population + i);
		eval_chrom(population + i);
		// by considering uniqueness, the process could go very long
		if (!is_duplicate(i, population[i], -1, 0))
			i++;
	}
	// if use other heuristic's result, put it as 1st individual
	if (seeding == 1) { 
		if ( (f=fopen(file1, "r")) == NULL ) {
			fprintf(stderr, "Error on reading ETC matrix\n");
			exit(1);
		}
		for (j=0; j<n; j++) {
			fscanf(f, "%d ", &(population[0].solution[j]));
		}
		fclose(f);
	}
}
// Feltl's constraint heuristic
void gen_init_chrom_constraint(Chrom * chrom)
{
	// Note: chrom must point to valid memory
	int i, j, bi, bin, posj, posi;
	int cand_bins[m]; long long cap[m];
	memset(cap, 0, sizeof(long long) * m);
	posj = MYRANDI(n);
	for (j=0; j<n; j++) {
		bi = 0;
		// find bins that can hold item j
		posi = MYRANDI(m);
		for (i=0; i<m; i++) {
			if (M(w, posi, posj) + cap[posi] <= b[posi]) {
				cand_bins[bi] = posi;
				bi ++;
			}
			posi = (posi + 1) % m;
		}
		// randomly select one that fits j
		if (bi == 0) {
			bin = MYRANDI(m);
		} else {
			bin = cand_bins[MYRANDI(bi)];
		}
		// assign
		(chrom->solution)[posj] = bin;
		// update weigh capacity
		cap[bin] += M(w, bin, posj);
		posj = (posj + 1) % n;
	}
}
// A modified version of Chu's feasibility improvement method
// Original method finds the 1st replacement agent when an agent's
// capacity is exceeded. But we get replacement agent list first,
// then randomly assign an agent. Purpose: increase randomness.
void improve_feasibility_chu(Chrom * chrom)
{
	long long cap[m];
	int items[n], item, ii, i, j, k, kk, loop;
	char mark[n];
	get_cap(chrom, cap);
	// start from a random bin
	i = MYRANDI(m);
	for (loop=0; loop<m; loop++) {
		// get the list of items in the bin
		ii = 0;
		for (j=0; j<n; j++) {
			if ((chrom->solution)[j] == i) {
				items[ii] = j; ii ++;
			}
		}
		for (k=0; k<ii; k++)
			mark[k] = ' '; // not marked
		k = 0;
		while (cap[i] > b[i]) {
			// if cap is exceeded, find a random unmarked item
			k = MYRANDI(ii); j=0;
			while (mark[k] == 'x' && j<ii) {
				k = (k+1) % ii;
				j++;
			}
			if (j==ii) break; // failed to improve feasibility
			item = items[k];
			mark[k] = 'x';
			// try to assign it to another bin
			// TODO: If moving the item does not work, how about
			//	   exchanging it with some item in another bin?
			//kk = MYRANDI(m);
			for (k=0; k<m-1; k++) {
				kk = (k+i+1) % m;
				if (cap[kk] + M(w, kk, item) <= b[kk]) {
					(chrom->solution)[item] = kk;
					cap[kk] += M(w, kk, item);
					cap[i] -= M(w, i, item);
					break;
				}
			}
		}
		i ++; if (i >= m) i = 0;
	}
	if (debug) eval_chrom(chrom);
}
void improve_quality_chu(Chrom * chrom)
{
	//TODO: avoid repeated get_cap() calls in search()
	long long cap[m];
	int i, j, min, minv, bin, loop;
	get_cap(chrom, cap);
	//increase randomness by starting at a random j
	j = MYRANDI(n);
	for (loop=0; loop<n; loop++) {
		// foreach item, find a reassignment
		// that has better fitv and still satisfied
		// TODO: How about exchange, instead of reassignment?
		bin = (chrom->solution[j]);
		minv = -1;
		for (i=0; i<m; i++) { // use trasnposed M to speed up
			if (i!=bin && (minv==-1 || minv>MT(vt, j, i)) && \
				   MT(vt, j, i) < MT(vt, j, bin) && \
				   MT(wt, j, i) + cap[i] <= b[i] ) { 
				min = i;
				minv = MT(vt, j, i);
			}
		}
		if (minv >= 0) {
			// reassign item to the better bin
			(chrom->solution)[j] = min;
			cap[bin] -= M(w, bin, j);
			cap[min] += M(w, min, j);
		}
		j ++; if (j >= n) j = 0;
	}
	if (debug) eval_chrom(chrom);
}
// util to get each bin's capacity in a chrom
void get_cap(Chrom * chrom, long long *cap)
{
	int j;
	memset(cap, 0, sizeof(long long) * m);
	for (j=0; j<n; j++) {
		cap[chrom->solution[j]] += M(w, (chrom->solution)[j], j);
	}
}
// check if a chrom is unique in history. this is important
// because duplicate chroms lead to premature solutions
// TODO: need a memory-efficient implementation to truly decide
// if a chrom is unique among all temporally generated chroms
// not only current chroms
int is_duplicate(int size, Chrom chrom, int inclusive, int checkAll)
// size may not be pop_size
// requirement: pop must be evaluated already.
// return value: 0-unique, 1-duplicate, >1-converged(all the same)
{
	int i, j, duplicate=0;
	if (size < 0) return 0;
	for (i=0; i<size; i++) {
		if (population[i].ev[EV_FITV] == chrom.ev[EV_FITV] && (inclusive != i)) {
			for (j=0; j<n; j++) {
				if (population[i].solution[j] != chrom.solution[j])
					break;
			}
			if (j==n) {
				if (checkAll)
					duplicate++;
				else
					return 1;
			}
		}
	}
	return duplicate;
}

// evaluation function
void eval_chrom(Chrom *chrom)
{
	long long ov = 0, u, capacity[m];
	int j, i;
	// calc obj func value and capacity value 
	memset(capacity, 0, m * sizeof(long long));
	for (j=0; j<n; j++) {
		ov += M(v, (chrom->solution)[j], j);
		capacity[(chrom->solution)[j]] += M(w, (chrom->solution)[j], j);
	}
	// obj func value
	chrom->ev[EV_OBJV] = ov;
	// fitness value
	chrom->ev[EV_FITV] = ov;
	// unfitness value
	u = 0;
	for (i=0; i<m; i++) {
		// feasible solution always fits
		u += ((capacity[i] - b[i]<=0)?0:capacity[i]-b[i]);
	}
	chrom->ev[EV_UFITV] = u;
}
void solution_verify(Chrom chrom)
{
	int i, j, faulty = 0;
	long long ov;
	long long cap[m];
	get_cap(&chrom, cap);
	ov = 0;
	for (j=0; j<n; j++) {
		ov +=  M(v, (chrom.solution)[j], j);
	}
	if (ov != chrom.ev[EV_OBJV]) {
		fprintf(stderr, "wrong: objv(%lld != %lld) ", ov,  chrom.ev[EV_OBJV]);
		faulty = 1;
	}

	for (i=0; i<m; i++) {
		if (cap[i] > b[i]) {
			fprintf(stderr, "wrong: cap violation[%d]: %lld > %d", i, cap[i], b[i]);
			faulty = 1;
		}
	}
	if (faulty) fprintf(stderr, "\n"); 
}
// evaluate the whole population
int eval_pop()
{
	int i, elite;
	long long min_fitv, max_fitv1, sum_fitv1;
	min_fitv=0;
	elite = 0;
	max_fitv1 = 0;
	sum_fitv1 = 0;
	
	// calc fitness value for each chromosome, keep elite 
	for (i=0; i<pop_size; i++) {
		//eval_chrom(population+i); // assume all individuals are evaluated
		if (max_fitv1 < population[i].ev[EV_FITV]) max_fitv1 = population[i].ev[EV_FITV];
		sum_fitv1 += population[i].ev[EV_FITV];
		if (min_fitv == 0 || population[i].ev[EV_FITV] < min_fitv) {
			min_fitv = population[i].ev[EV_FITV];
			elite = i;
		}
	}
	gastat.min_fitv = min_fitv;	
	gastat.max_fitv = max_fitv1;
	gastat.avg_fitv = sum_fitv1 / pop_size;

	return elite;
}
// update a pop's stat
void eval_pop_update(Chrom newc, Chrom oldc)
{
	if (newc.ev[EV_FITV] < gastat.min_fitv)
		gastat.min_fitv = newc.ev[EV_FITV];
	if (newc.ev[EV_FITV] > gastat.max_fitv)
		gastat.max_fitv = newc.ev[EV_FITV];
	gastat.avg_fitv = gastat.avg_fitv + (newc.ev[EV_FITV] - oldc.ev[EV_FITV])/pop_size;
}

// rand whole population
void rank_pop(int which) 
{
	int j;
	for (j=0;j<pop_size;j++)
		rank[which][j] = j;
	quicksort(0, pop_size - 1, which);
}
// sort chromosome, i.e., rearrange index in rank[]
// pivotIndex IS NOT USED, always keep it -1
void quicksort (int lo, int hi, int which)
{
	int mid, i;
	long long pivot;
	if (lo >= hi)
		return;
	pivot = population[rank[which][lo]].ev[which];
	mid = lo;
	for (i = lo + 1; i <= hi; i++)
		if (population[rank[which][i]].ev[which] < pivot) {
			mid++;
			swap (mid, i, which);
		}
	swap (lo, mid, which);
	quicksort (lo, mid - 1, which);
	quicksort (mid + 1, hi, which);
}
void quicksort_verify()
{
	int i, j;
	long long v1 , v2;
	for (i=0; i<EV_COUNT; i++) {
		v1 = 0; // assumer all values are >=0
		for (j=0; j<pop_size; j++) {
			v2 = population[rank[i][j]].ev[i];
			if (v1 > v2) {
				fprintf(stderr, "wrong [%s,%d,%d]: %lld > %lld\n", ev_text[i], j-1, j, v1, v2);
				return;
			}
			v1 = v2;
		}
	}
}

void swap(int i, int j, int which)
{
	int tmp;
	if (i==j) return;
	tmp = rank[which][i]; 
	rank[which][i] = rank[which][j]; 
	rank[which][j] = tmp;
}
// sort when the Chrom in rank[randIndex] is changed
void insertsort(int which)
{
	int tmp;
	int pivot; // wtb: where_to_be
	// find the pivot (the one that screwed order)
	pivot = pop_size - 1; // 'cause usually we replace the last one for one of the measures
	while (pivot > 0 && population[rank[which][pivot]].ev[which] >= population[rank[which][pivot - 1]].ev[which])
		pivot --;
	if (pivot == 0) return;
	// now Xi-1>Xi, find out which one is the real pivot
	swap(pivot, pivot-1, which);
	if (pivot < pop_size - 1 && population[rank[which][pivot]].ev[which] > population[rank[which][pivot+1]].ev[which]) {
		// Xi-1 (now Xi) is the pivot
		tmp = rank[which][pivot];
		while (pivot < pop_size - 1 && population[tmp].ev[which] > population[rank[which][pivot+1]].ev[which]) {
			rank[which][pivot] = rank[which][pivot + 1];
			pivot ++;
		}
		rank[which][pivot] = tmp;
	}else {
		// Xi (now Xi-1) is the pivot
		pivot = pivot - 1;
		tmp = rank[which][pivot];
		while (pivot > 0 && population[tmp].ev[which] < population[rank[which][pivot-1]].ev[which]) {
			rank[which][pivot] = rank[which][pivot - 1];
			pivot --;
		}
		rank[which][pivot] = tmp;
	}

/* Following impl doesn't apply for all measures
	// move to end
	while (wtb < pop_size-1 && population[rank[which][rankIndex]].ev[which] > population[rank[which][wtb+1]].ev[which]) 
		wtb ++;
	if (wtb > rankIndex) {
		// one big swap
		tmp = rank[which][rankIndex];
		while (rankIndex < wtb) {
			rank[which][rankIndex] = rank[which][rankIndex + 1];
			rankIndex ++;
		}
		rank[which][wtb] = tmp;
	} else { // move to front
		while (wtb > 0 && population[rank[which][rankIndex]].ev[which] < population[rank[which][wtb-1]].ev[which])
			wtb --;
		if (wtb < rankIndex) {
			 // one big swap
			tmp = rank[which][rankIndex];
			while (rankIndex > wtb) {
				rank[which][rankIndex] = rank[which][rankIndex - 1];
				rankIndex --;
			}
			rank[which][wtb] = tmp;
		}
	}
*/
}

// get rank of a chromosome 
int get_wheel_selection()
{
	double area = wheel[0];
	int i;
	double rand_num = MYRANDF();
	for (i=0; i<pop_size-1; i++)
	{
		if (rand_num < area)
			return i;
		else area += wheel[i+1];
	}
	return (pop_size - 1);
}
// select two chromosomes as parents
// return the index of each one in population.
// we always use tournament selection: pick several candidates,
// then return a parent; then do it again to get the other
void selection(int *p1, int *p2)
{
	switch (strategy_selection) {
		case 'r': // rank-based
			rank_tournament(p1, p2, EV_FITV);
			break;
		case 'b': // binary selection
			binary_tournament(p1, p2);
			break;
		default:  // random selection
			*p1 = MYRANDI(pop_size);
			*p2 = MYRANDI(pop_size);
			if (*p1 == *p2) *p2 = (*p1 + 1 + MYRANDI(pop_size - 1)) % pop_size;
			break;
	}
}
void rank_tournament(int *p1, int *p2, int which)
{
	int ir1, ir2;
	// look up rank table, get rand index of the random number 
	ir1 = get_wheel_selection();
	ir2 = get_wheel_selection();
	// make sure they are different
	if (ir1 == ir2) ir2 = (ir1 + 1 + MYRANDI(pop_size - 1)) % pop_size;
	// return the chrom indices
	*p1 = rank[which][ir1];
	*p2 = rank[which][ir2]; 
}
void binary_tournament(int *p1, int *p2)
{
	int i1, i2;
	// randomly pick two chroms
	i1 = MYRANDI(pop_size);
	i2 = MYRANDI(pop_size);
	if (i1 == i2) i2 = (i1 + 1 + MYRANDI(pop_size - 1)) % pop_size;
	*p1 = (population[i1].ev[EV_FITV] < population[i2].ev[EV_FITV])?i1:i2;
	i1 = MYRANDI(pop_size);
	i2 = MYRANDI(pop_size);
	if (i1 == i2) i2 = (i1 + 1 + MYRANDI(pop_size - 1)) % pop_size;
	*p2 = (population[i1].ev[EV_FITV] < population[i2].ev[EV_FITV])?i1:i2;

}

void crossover(Chrom parent1, Chrom parent2, Chrom *child)
{
	// randomly select a cut-off point
	int cut = MYRANDI(n);
	double d = MYRANDF();
	int ff = (d>0.5)?1:0;
	int i;
	for (i=0; i<n; i++) {
		if (i < cut) 
			(child->solution)[i] = (ff==0)?parent1.solution[i]:parent2.solution[i];
		else
			(child->solution)[i] = (ff==0)?parent2.solution[i]:parent1.solution[i];
	}
	//TODO: crossover can produce two offspring, should we look at the 2nd?
	if (debug) eval_chrom(child);
}
void mutate(Chrom *child)
{
	int i1 = MYRANDI(n);
	int i2 = MYRANDI(n);
	if (i1 == i2) i2 = (i1 + 1 + MYRANDI(n - 1)) % n;
	int bin = child->solution[i1];
	child->solution[i2] = child->solution[i1];
	child->solution[i1] = bin;
	if (debug) eval_chrom(child);
}
// replacement
// return the index of child in population. -1: no replacement
int replacewith(Chrom child)
{
	int newRPos = -1;
	Chrom temp; 
	temp.solution = (int *)malloc(sizeof(int) * n);
	if (strategy_replacement == 'u' && population[rank[EV_UFITV][pop_size - 1]].ev[EV_UFITV] > 0) {
		// eliminate unfittest individual
		if (child.ev[EV_UFITV] < population[rank[EV_UFITV][pop_size - 1]].ev[EV_UFITV]) {
			chrom_copy(population[rank[EV_UFITV][pop_size - 1]], &temp);
			chrom_copy(child, &(population[rank[EV_UFITV][pop_size - 1]]));
			// we can do this earlier to save space, just keep logic clear
			eval_pop_update(child, temp);
			//newIndex = rank[EV_UFITV][pop_size - 1];
			newRPos = pop_size - 1;
		}
	} else { // 'w' or all feasible
		// eliminate worst individual
		if (child.ev[EV_FITV] < population[rank[EV_FITV][pop_size - 1]].ev[EV_FITV]) {
			chrom_copy(population[rank[EV_FITV][pop_size - 1]], &temp);
			chrom_copy(child, &(population[rank[EV_FITV][pop_size - 1]]));
			// we can do this earlier to save space, just keep logic clear
			eval_pop_update(child, temp);
			//newIndex = rank[EV_FITV][pop_size - 1];
			newRPos = pop_size - 1;
		}
	}
	free(temp.solution);
	return newRPos;
}
// utils
void roulette_init()
{
	int i;
	// initialize wheel. \cite{Wang-1997}
	float R = 1.0 + 1.0/pop_size; // ratio of rank
	double wheeltemp = 0.0; 
	// roulette wheel construction
	// each chrom's angle on the wheel depends on only pop_size,
	// not a chrom's fitness value. it only says higher ranked
	// chrom gets higher probability to be selected.
	for (i=0; i<pop_size; i++) {
		wheel[i] = pow(R, (pop_size-i-1)) * (R-1) / (pow(R, pop_size) - 1);
		wheeltemp+=wheel[i];
	}
	// debug 
	if (debug) {
		for (i=0;i<pop_size;i++) fprintf(myout, "%lf ", wheel[i]);
		fprintf(myout, "\nsum of wheel sectors: %lf\n", wheeltemp);
	}
}
double get_ga_time()
{
	struct timeval tsec;
	struct timezone tzone;
	gettimeofday(&tsec, &tzone);
	return (double)(tsec.tv_sec + tsec.tv_usec/1000000.0);
}
void chrom_copy(Chrom c1, Chrom *c2)
{
	int i;
	for (i=0; i<n; i++)
		c2->solution[i]=c1.solution[i];
	for (i=0; i<EV_COUNT; i++)
		c2->ev[i]=c1.ev[i];
}
void print_pop(int howmany)
{
	int i; char title[80];
	fprintf(myout, "********Population**********\n");
	if (howmany == 0) howmany=pop_size;
	for (i=0; i<howmany; i++) {
		sprintf(title, "chromosome%d", i);
		print_chrom(population[i], title);
	}
	fprintf(myout, "********End**********\n");
}
void print_chrom(Chrom chrom, char * title)
{
#ifndef NOIO
	int j, i; 
	int bins[m][n];
	int bi[m];
	long long cap[m];
	memset(cap, 0, sizeof(long long) * m);
	memset(bi, 0, sizeof(int) * m);
	for (j=0; j<n; j++) {
		cap[chrom.solution[j]] += M(w, (chrom.solution)[j], j);
		bins[chrom.solution[j]][bi[chrom.solution[j]]] = j;
		bi[chrom.solution[j]] ++;
	}
	// title='Z' means detailed output
	if (title!=NULL && *title=='Z') {
		fprintf(myout, "-----------%s-------------\n", title);
		for (i=0; i<m; i++) {
			fprintf(myout, "bin %d[%d, %lld] ", i, b[i], cap[i]);
			for (j=0; j<bi[i]; j++) {
				fprintf(myout, "%d ", bins[i][j]);
			}
			fprintf(myout, "\n");
		}
	}
	//fprintf(myout, "\n=====chromosome%s[objv,fitv,ufitv,iter,bestT] %lld %lld %lld %lld %lf =====\n", (title==NULL?"":title), chrom.ev[EV_OBJV], chrom.ev[EV_FITV], chrom.ev[EV_UFITV], iterations, (gastat.bestT-gastat.startT));
	fprintf(myout, "\n=====chromosome%s[] %lld %lld %lld %lld %lf =====\n", (title==NULL?"":title), chrom.ev[EV_OBJV], chrom.ev[EV_FITV], chrom.ev[EV_UFITV], iterations, (gastat.bestT-gastat.startT));
	
 /*
	for (j=0; j<n; j++) {
		if ((j % 10) == 0) printf("\n");
		fprintf(myout, "%d->%d; ", j, chrom.solution[j]);
	}
	fprintf(myout, "\n");	   
*/
	fflush(myout);
#endif
}
void print_stat(FILE * myout)
{
#ifdef PGAMODE
	//fprintf(myout, "=====stat[totalImpv,migImpv,bestT,T,commT, min,max,avg] %d %d %lf %lf %lf %lld %lld %lld\n", gastat.total_improve, gastat.mig_improve, gastat.bestT - gastat.startT, gastat.endT - gastat.startT, gastat.commT, gastat.min_fitv, gastat.max_fitv, gastat.avg_fitv);
	// totalImpv, migImprove, bestT, sendT, recvT
	fprintf(myout, "=====statZ%dZ %d %d %lf %lf %lf ", myrank, gastat.total_improve, gastat.mig_improve, gastat.bestT - gastat.startT, gastat.comm_sendT, gastat.comm_recvT);
#ifdef T_PROFILING
	fprintf(myout, "%lf %lf %lf %lf %lf %lf %lf %lf ", gastat.immigrateT, gastat.selectionT, gastat.injectT, gastat.crossmutT, gastat.feasibT, gastat.improveT, gastat.evalT, gastat.replaceT);
#endif
	fprintf(myout, "\n");
#else
	fprintf(myout, "=====stat[totalImpv,migImpv,T,bestT,min,max,avg] %d %d %lf %lf %lld %lld %lld\n", gastat.total_improve, 0, gastat.bestT - gastat.startT, gastat.endT - gastat.startT, gastat.min_fitv, gastat.max_fitv, gastat.avg_fitv);
#endif
	fflush(myout);
}
void print_rank(int ev_i)
{
	int ev_index[EV_COUNT], i, j, k;
	if (ev_i <0 || ev_i>EV_COUNT) return;
	memset(ev_index, 0, sizeof(int) * EV_COUNT);
	if (ev_i == EV_COUNT) {
		for (i=0; i<EV_COUNT; i++)
			ev_index[i] = 1;
	} else {
		ev_index[ev_i] = 1;
	}
	for (i=0; i<EV_COUNT; i++) {
		if (ev_index[i]) {
			fprintf(myout, "%s I ", ev_text[i]);
			k=0;
			for (j=0; j<pop_size; j++) {
				fprintf(myout, "%d\t", rank[i][j]);
				if (k==9) {
					k=0; // format control
					fprintf(myout, "\n");
					fprintf(myout, "%s I ", ev_text[i]);
				} else k++;
			}
			fprintf(myout, "\n");
			fprintf(myout, "%s V ", ev_text[i]);
			k=0;
			for (j=0; j<pop_size; j++) {
				fprintf(myout, "%lld\t", population[rank[i][j]].ev[i]);
				if (k==9) {
					k=0; // format control
					fprintf(myout, "\n");
					fprintf(myout, "%s V ", ev_text[i]);
				} else k++;
			} 
			fprintf(myout, "\n");
		}
	}
	fprintf(myout, "\n");
}
void output_result(int * solution, long long objv, double exec_time)
{
	char assignment[n*10+3], astr[12];
	int i;
	for (i=0; i<n; i++) {
		sprintf(astr,"%d ", solution[i]);
		strcat(assignment, astr);
	}
	glog("---Result: objv=%lld, exectime=%lf", objv, exec_time);
	glog("%s", assignment);
}
#ifdef PGAMODE
// put 2 selected chroms in emigration buffer
int emigrate(int * emi_buffer, int newimprove, int noimprove)
{
	int n_elite, n_rand, origin;
	// emigrant selection
	if (newimprove >= 5) {
		// if lots of improvements, send 2 elites
		n_elite = 2; n_rand = 0;
	} else if (newimprove > 0) {
		// has improvements, send 1 elite & 1 rand
		n_elite = 1; n_rand = 1;
	} else {
		// no improvements, we should not emigrate
		// but to inject diversity, we do it w/ 1/20 chance
		n_elite = 0; n_rand = (MYRANDI(100)>5)?0:2;
		// if there had been no improvement, why not inject
		// some noise to others? this could stir up new evolutions
		if (n_rand == 0 && noimprove < 2000)
			return 0; // don't emigrate
	}
	// fill emi_buffer
	int *emi_index = emi_buffer;
	// copy elites
	int i, j;
	for (i=0; i<n_elite; i++) {
		*emi_index = MIG_ELITE; emi_index++; // chrom type
		// history only records chroms with UFITV=0,so unfeasible chrom
		// will always return -1
		origin = imi_find_origin(population[rank[EV_FITV][i]].ev[EV_FITV]);
		*emi_index = (origin==-1)?myrank:origin; emi_index++; // origin rank
		*emi_index = send_seq; emi_index++; // send sequence number
		//TODO: choose 2 elites w/ larger diversity difference
		if (population[rank[EV_UFITV][i]].ev[EV_UFITV] > 0) {
			// no feasible solutions, fill least unfeasible ones
			for (j=0; j<n; j++) {
				*emi_index = population[rank[EV_UFITV][i]].solution[j];
				emi_index ++;
			}
		} else {
			// find best feasible sol'n
			for (j=0; j<n; j++) {
				*emi_index = population[rank[EV_FITV][i]].solution[j];
				emi_index ++;
			}
		}
	}
	// copy random
	int p1, p2;
	for (i=0; i<n_rand; i++) {
		binary_tournament(&p1, &p2);
		*emi_index = MIG_RANDOM; emi_index++; // chrom type
		*emi_index = myrank; emi_index++; // origin rank
		*emi_index = send_seq; emi_index++; // send sequence number
		for (j=0; j<n; j++) {
			*emi_index = population[p1].solution[j];
			emi_index ++;
		}
	}
	send_seq ++; // update sequence number
	return 1;
}
// return an immigrant from immigrant pool and update its processing state
// return also chrom type and its origin
int immigrate(Chrom * imi, int *origin)
{
	int t, orig;
	t = fetch_imi(imi->solution, &orig);
	*origin = orig;
	return t;
}
/*
void print_imi()
{
	int * imi_index = imi_buffer;
	int i,j;
	fprintf(myout, "+++imi_buffer");
	for (i=0; i<emi_size*neighbor_count; i++) {
		fprintf(myout, ".%d: %d", i, *imi_index);
		for (j=0; j<n; j++)
			fprintf(myout, " %d", *(imi_index+j+1));
		fprintf(myout, "\n");
		imi_index += mig_msglen;
	}
	fflush(myout);
}
*/
#endif
// ga search routine
void * search(void * args) 
{
	int i, j, i1, i2;
	//int converged = 0; // all chromosomes converge 
	int no_improve = 0;  // improved elite chromosome
#ifdef PGAMODE
	int new_improve = 0;  //count improved elite chromosome b/w migrations 
	int from_mig = 0; // mark migration-related iteration
	int imi_origin = 0;
	unsigned long long last_emig_iter = 0;
	unsigned long long last_imig_iter = 0;
	int imiType = 0;
	double commT1, commT2, commT3;
	int coord1[2], coord2[2], cart_dist; // MPI_Cart topo coordinates and dist
#endif
#ifdef T_PROFILING
	double T0, T1;
#endif
	Chrom *parent1, *parent2, child;
	Chrom elite_chrom;
	char buffer[256];
	// initialize random number generator
	if (debug) {
		global_seed = 0;
	} else {
#ifdef PGAMODE
		global_seed = 0; // required by SPRNG
#else
		if (global_seed < 0) // not set by command line
			global_seed = (int) (get_ga_time() + getpid());
#endif
	}
#ifdef PGAMODE
	// print configuration
	if (myrank == 0) 
		print_config();
#endif
#ifdef GSL_SPRNG
#ifdef PGAMODE
	//MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	//MPI_Comm_size( MPI_COMM_WORLD, &net_size);
	gsl_sprng_r = gsl_sprng_alloc(myrank, net_size, global_seed, NULL);
#else
	gsl_sprng_r = gsl_sprng_alloc(0, 1, global_seed, NULL);
#endif
#else
#ifdef SPRNG
	mysprng_init();
#else
	srand48(global_seed); 
#endif
#endif

	// initialize population arrays
	population = (Chrom *) malloc(pop_size * sizeof(Chrom));
	if (population == NULL) {
		fprintf(myout, "ERROR: could not malloc population memory of %ld bytes\n", pop_size * sizeof(Chrom));
		fflush(myout);
		exit (1);
	}
	for (j=0; j<pop_size; j++) {
		population[j].solution = (int *) malloc(n * sizeof(int));
	}
	child.solution = (int *) malloc(n * sizeof(int));
	elite_chrom.solution = (int *) malloc(n * sizeof(int));
#ifdef PGAMODE
	Chrom immigrant;
	immigrant.solution = (int *) malloc(n * sizeof(int));
#endif
	// rank's index is ranking, each element is the index of chrom in pop 
	for (i=0; i<EV_COUNT; i++) {
		rank[i] = (int *) malloc(pop_size * sizeof(int));
	}
	wheel = (double *) malloc(pop_size * sizeof(double));
	gastat_init();

	// GA: generate initial population 
	gen_init_pop(seeding);
	// GA: evaluate the whole population
	eval_pop();
	// GA: rank population
	for (i=0; i<EV_COUNT; i++) {
		rank_pop(i);
	}
	// init elite chrom. might not be feasible at beginning
	chrom_copy(population[rank[EV_FITV][0]], &elite_chrom);
	if (debug) {fprintf(myout, "rank_pop debug\n"); quicksort_verify();fprintf(myout, "rank_pop debug\n"); } 
//	if (debug)
//		print_rank(EV_COUNT);
	if (debug) print_pop(10);
	// GA: setup rank-based rouletee wheel
	roulette_init();   
	// GA: start evolution 
	iterations = 0;

#ifdef T_PROFILING
	double tickStart = get_ga_time();
	int tickCount =0;
	int tickInterval = 300; // output something every tickInterval seconds
#endif
	do
	{
#ifdef T_PROFILING
		T0 = get_ga_time();
#endif
#ifdef PGAMODE
		// if there are immigrants, take them as new children
		// but don't crossover and mutate
		imiType = immigrate(&immigrant, &imi_origin);
#endif
#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.immigrateT += T1 - T0;
		T0 = T1;
#endif
		// GA: selection
		selection(&i1, &i2);
#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.selectionT+= T1 - T0;
		T0 = T1;
#endif
		parent1 = population + i1;
#ifdef PGAMODE
		// if immigrant is randomly picked by neighbor, inject as a parent
		if (imiType == MIG_RANDOM) {
			improve_feasibility_chu(&immigrant);
			parent2 = &immigrant;
			from_mig = 1;
		} else {
#endif
		parent2 = population + i2;
#ifdef PGAMODE
		}
#endif
#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.injectT+= T1 - T0;
		T0 = T1;
#endif
		if (debug) print_chrom(*parent1, "parent1");	
		if (debug) print_chrom(*parent2, "parent2");
#ifdef PGAMODE
		// if immigrant is other island's elite, take as offspring
		if (imiType == MIG_ELITE) {
			chrom_copy(immigrant, &child);
			from_mig = 1;
			// performance study only
			eval_chrom(&child);
			commT3 = get_ga_time(); 
			// calc distance b/w self and chrom origin
			MPI_Cart_coords(topoComm, myrank, 2, coord1);
			MPI_Cart_coords(topoComm, imi_origin, 2, coord2);
			// Note: our topo is wired on boundaries, we only calc
			// the dist b/w two procs in one way, not necessarily the 
			// number of hops the chrom actually traveled.
			cart_dist = abs(coord1[0] - coord2[0]) + abs(coord1[1] - coord2[1]);
			if (pdebug) fprintf(myout, "+++imi(time,iter,bestSolnQ,distance): %lf %lld %lld %d (%d-%d,<%d,%d>-<%d,%d>)\n", commT3, iterations, child.ev[EV_FITV], cart_dist, imi_origin, myrank, coord2[0],coord2[1],coord1[0],coord1[1]);
		} else {
#endif
			// GA: crossover 
			crossover(*parent1, *parent2, &child);
			// GA: mutation 
			mutate(&child);
			if (debug) print_chrom(child, "after mating");   
#ifdef PGAMODE
		}
#endif
#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.crossmutT += T1 - T0;
		T0 = T1;
#endif
		// GA: improve feasibility
		improve_feasibility_chu(&child);
#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.feasibT += T1 - T0;
		T0 = T1;
#endif
		if (debug) print_chrom(child, "after improve feasi");	
		// GA: improve solution quality
		improve_quality_chu(&child);
#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.improveT += T1 - T0;
		T0 = T1;
#endif
		if (debug) print_chrom(child, "after improve quality");	
		// GA: evaluation
		eval_chrom(&child);
#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.evalT += T1 - T0;
		T0 = T1;
#endif
		// count improvement: feasible and improved fitness
		if (child.ev[EV_UFITV] == 0 && population[rank[EV_FITV][0]].ev[EV_FITV] > child.ev[EV_FITV]) {
			no_improve = 0;
#ifdef PGAMODE
			new_improve ++;
			if (from_mig) {
				gastat.mig_improve ++;
			}
#endif
			gastat.total_improve ++;
			gastat.bestT = get_ga_time();
			chrom_copy(child, &elite_chrom);
			sprintf(buffer, "New elite at iteration %lld", iterations);
			//print_chrom(child, buffer);
#ifdef PGAMODE
			if (from_mig) {
				print_chrom(child, "Mig");
				// put it in immigrant history buffer
				imi_hist_chrom[imi_hist_index] = child.ev[EV_FITV];
				imi_hist_origin[imi_hist_index] = imi_origin;
				imi_hist_index = (imi_hist_index + 1) % IMI_HIST_BUFFER_SIZE;
			} else {
#endif
			print_chrom(child, NULL); // for performance study only, skip solution printing
#ifdef PGAMODE
			}
#endif
#ifndef NOIO
			print_stat(myout);
#endif
			solution_verify(child);
		} else {
			no_improve++;
		}
		if (!is_duplicate(pop_size, child, 0, 0)) {
		// GA: replacement
			int oldRank = replacewith(child);
		// rank population according to fitness value
			if (oldRank >= 0) { 
				for (i=0; i<EV_COUNT; i++) {
					insertsort(i); // rank update
				}
				if (debug) quicksort_verify();
			}
		}

#ifdef T_PROFILING
		T1 = get_ga_time();
		gastat.replaceT += T1 - T0;
		T0 = T1;
#endif
#ifdef PGAMODE
		if (from_mig) from_mig = 0; // reset mig improve indicator
#endif		
		iterations++;

#ifdef PGAMODE
		// emigration and immigration might overlap in async mode	   
		if (iterations - last_emig_iter >= migrate_freq) {
			commT1 = get_ga_time(); 
#ifdef PGAMODE_SYNC
			MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef PGA_THREAD_MODE
			int ei;
			if (pthread_mutex_trylock(&emi_mutex) == 0) {
				if (!emi_ready) { // no need to check, just write?
					for (ei=0; ei<emi_size; ei++) {
						//chrom_copy(..., &(emi_buffer[ei]));
					}
					emi_ready = 1;
				}
				pthread_mutex_unlock(&emi_mutex);
			}
			if (pthread_mutex_trylock(&imi_mutex) == 0) {
				if (imi_ready) {
					// process immigrants

					imi_ready = 0;
				}
				pthread_mutex_unlock(&imi_mutex);
			}
			//sched_yield(); // yield on single cpu arch
#else
			// emigrate
			int * emi_buffer;
			if ((emi_buffer = send_req()) != NULL && emigrate(emi_buffer, new_improve, no_improve)) {
				if (new_improve)
					if (pdebug) fprintf(myout, "+++emi(time,iter,bestSolnQ,type): %lf %lld %lld %d new improves\n", commT1, iterations, elite_chrom.ev[EV_OBJV], new_improve);
#ifdef PGA_NONBLOCK_MODE
				//send
				send_emi(); // the actuall send
#else
				// block mode: use MPI_Gather
				MPI_Comm_rank(migComm, &myrank);
				MPI_Barrier(migComm);
				MPI_Gather(emi_buffer, emi_size*mig_msglen, MPI_INT, imi_buffer, emi_size * mig_msglen, MPI_INT, myrank, migComm);
				MPI_Barrier(migComm);
				//if (debug) print_imi();
#endif
			}
			commT2 = get_ga_time(); 
			gastat.comm_sendT += commT2 - commT1;
			gastat.commT += commT2 - commT1; 
			new_improve = 0; // reset
			last_emig_iter = iterations;
#ifdef PGA_NONBLOCK_MODE
		}
		if (iterations - last_imig_iter >= immigrate_freq) {
			commT1 = get_ga_time(); 
#ifdef PGAMODE_SYNC
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			recv_imi();
			commT2 = get_ga_time(); 
			gastat.comm_recvT += commT2 - commT1;
			gastat.commT += commT2 - commT1; 
			last_imig_iter = iterations;
#endif
		}
#endif
#endif
#ifdef T_PROFILING
		double tickT = get_ga_time();
		if (tickT - tickStart > tickInterval) {
			tickCount ++;
#ifdef PGAMODE
			fprintf(myout, "%d tick mark %d seconds\n", myrank, tickInterval * tickCount);
			fprintf(stderr, "%d tick mark %d seconds\n", myrank, tickInterval * tickCount);
#else
			fprintf(myout, "%d tick mark %d seconds\n", 0, tickInterval * tickCount);
			fprintf(stderr, "%d tick mark %d seconds\n", 0, tickInterval * tickCount);
#endif
			fflush(myout);
			tickStart = tickT;
		}
#endif
	// GA: stopping rules 
	} while ((strategy_stop=='f' && iterations<max_iterations) || (strategy_stop=='i' && no_improve < improve_threshold) || (strategy_stop=='q' && stopping_quality>0 && ((elite_chrom.ev[EV_FITV]>stopping_quality && elite_chrom.ev[EV_UFITV]==0) || (elite_chrom.ev[EV_UFITV]!=0))) || (strategy_stop=='t' && (get_ga_time() - gastat.startT <= walltime)));
		
	gastat.endT = get_ga_time();
	
#ifndef NOIO
	// output results (elite) 
	print_chrom(elite_chrom, "ZTheBestSolution");
	fprintf(myout, "Stat-iterations: %lld iterations\n", iterations);
	fprintf(myout, "Stat-improve_factor: %d\n", no_improve);
	fprintf(myout, "Stat-time: %lf\n", gastat.endT - gastat.startT);
	fprintf(myout, "Stat-best-solution-time: %lf\n", gastat.bestT - gastat.startT);
	// print_chrom(elite_chrom[pre_pop], task_num); 
	//output_result(population[rank[EV_OBJV][0]].solution, population[rank[EV_OBJV][0]].ev[EV_OBJV], gastat.endT - gastat.startT);
#endif

	// output to stdout too
#ifdef PGAMODE
	fprintf(stdout, "\n=====chromosome[Z%dZ] %lld %lld %lld %lld %lf %lf=====\n", myrank, elite_chrom.ev[EV_OBJV], elite_chrom.ev[EV_FITV], elite_chrom.ev[EV_UFITV], iterations, (gastat.bestT-gastat.startT), (gastat.endT - gastat.startT));
#else
	fprintf(stdout, "\n=====chromosome[ZZ] %lld %lld %lld %lld %lf %lf=====\n", elite_chrom.ev[EV_OBJV], elite_chrom.ev[EV_FITV], elite_chrom.ev[EV_UFITV], iterations, (gastat.bestT-gastat.startT), (gastat.endT - gastat.startT));
#endif
	fflush(stdout);
	print_stat(stdout);

#ifdef PGAMODE
	// wait until other processes are done if we care
	if (strategy_stop != 'q')	
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	// free resources
	for (j=0; j<pop_size; j++) {
		free(population[j].solution);
	}
	free(population);
	free(child.solution);
	free(elite_chrom.solution);
#ifdef PGAMODE
	free(immigrant.solution);
#endif
	free(wheel);
	for (i=0; i<EV_COUNT; i++) {
		free(rank[i]);
	}
#ifdef GSL_SPRNG
	gsl_sprng_free (gsl_sprng_r);
#endif
#ifdef SPRNG
	mysprng_free();
#endif
#ifdef PGAMODE
	ga_done = 1;
#endif
	return NULL;
}

int main(int argc, char ** argv)
{
#ifdef PGAMODE
	MPI_Init( &argc, &argv );
#endif
	strcpy(dataset, "");
	strcpy(myoutputdir, "");
	config(argc, argv);

	// set output file descriptor
	myout = stdout;
	char outputpath[1024];
	strcpy(outputpath, "");
#ifndef NOIO
#ifdef PGAMODE
	// must set
	if (strlen(myoutputdir)==0) strcpy(myoutputdir, ".");
	int myrank1;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank1);
	sprintf(outputpath, "%s/%s.%d.myout", myoutputdir, glogPrefix, myrank1);
#else
	if (strlen(myoutputdir)>0) {
		sprintf(outputpath, "%s/%s.%d.myout", myoutputdir, glogPrefix, getpid());
	} else {
		strcpy(outputpath, "");
	}
#endif 
#endif
	if (strlen(outputpath)>0) {
		if ((myout=fopen(outputpath, "w")) == NULL) {
			fprintf(stderr, "error open %s for writing\n", outputpath);
			myout = stdout;
		}
	}
	//glog("Start of GAP GA. pid: %d", (int)getpid());

	// read data
	readData(dataset);
	//glog("Dataset %s is read", argv[7]);

#ifdef PGAMODE	
	psearch();
#else
	search(NULL);
#endif

	if (myout != stdout)
		fclose(myout);

#ifdef PGAMODE
	// if we are not pursuing a bound, let others finish before quit
	if (strategy_stop != 'q')	
		MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	return 0;
}

// get command line parameters
void config(int argc, char **argv)
{
	glogId = globalId = 0;
	int i; 
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
		switch (argv[i][1]) {
			case 's':
				i++;
				if (i >= argc || (sscanf (argv[i], "%c", &strategy_selection) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
			case 'r':
				i++;
				if (i >= argc || (sscanf (argv[i], "%c", &strategy_replacement) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
			case 'i':
				i++;
				if (i >= argc || (sscanf (argv[i], "%c", &strategy_init_pop) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
			case 'p':
				i++;
				if (i >= argc || (sscanf (argv[i], "%d", &pop_size) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
			case 't':
				i++;
				if (i >= argc || (sscanf (argv[i], "%c", &strategy_stop) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				i++;
				long long number;
				if (i >= argc || (sscanf (argv[i], "%lld", &number) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				}
				switch (strategy_stop) {
					case 'f':
						max_iterations = (int)number;
						break;
					case 'q':
						stopping_quality = number;
						break;
					case 't':
						walltime = number;
						break;
					case 'i':
					default:
						improve_threshold = (int)number; 
						break;
				}
				break;
			case 'f':
				i++;
				if (i >= argc ) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				}
				strcpy(dataset, argv[i]); 
				break;
			case 'o':
				i++;
				if (i >= argc ) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				}
				strcpy(myoutputdir, argv[i]);
				break;
			case 'l':
				i++;
				if (i >= argc ) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				}
				sprintf(glogPrefix, "%s", argv[i]); 
				glogId = (int)getpid(); 
			break;
/*
			case 'c':
				i++;
				if (i >= argc || (sscanf (argv[i], "%f", &Pcross) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
			case 'm':
				i++;
				if (i >= argc || (sscanf (argv[i], "%f", &Pmut) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
*/
			case 'e':
				elitism = 1; 
				break;
			case 'd':
				debug = 1; 
				break;
/*
			case 'q':
				i++;
				if (i >= argc || (sscanf (argv[i], "%lld", &stopping_quality) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
*/
#ifdef PGAMODE
/*
			case 'B':
				comm_blocking = 1; 
				break;
			case 'S':
				comm_sync = 1; 
				break;
			case 'I':
				i++;
				if (i >= argc || (sscanf (argv[i], "%c", &imi_strategy) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
*/
			case 'R':
				i++;
				if (i >= argc || (sscanf (argv[i], "%ld", &global_seed) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				break;
			case 'M':
				i++;
				if (i >= argc || (sscanf (argv[i], "%d", &migrate_freq) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				i++;
				if (i >= argc || (sscanf (argv[i], "%d", &emi_size) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				}
				i++;
				if (i >= argc || (sscanf (argv[i], "%d", &immigrate_freq) != 1)) {
					//printf ("\"%s -?\" for help.\n", argv[0]);
					//exit (0);
					immigrate_freq = migrate_freq / 2; // default value
				} 
#ifdef HETERO
				// MIC processor: speed up export and import
				migrate_freq /= HETERO;
				immigrate_freq /= HETERO;
#endif
				break;
			case 'P':
				i++;
				if (i >= argc || (sscanf (argv[i], "%d", &snd_parallelism) != 1)) {
					printf ("\"%s -?\" for help.\n", argv[0]);
					exit (0);
				} 
				if (snd_parallelism <1) snd_parallelism = 1;
#ifdef HETERO_BUFFERCAPDIFF
				// MIC processor: hold less sends before mpi wait
				snd_parallelism /= HETERO_BUFFERCAPDIFF;
#endif
				break;
#endif
			case 'h':
			case '?':
			default:
				//printf ("Usage: %s [-s r|b] [-r u|n] [-i n|c|r] [-p pop_size] [-t i|f number] [-d data] [-c Pcross] [-m Pmute] [-e] [-B] [-S] [-M freq rate] \n", argv[0]);
				printf ("Usage: %s [-s r|b] [-r u|n] [-i n|c|r] [-p pop_size] [-t i|f number] [-f data] [-d] [-e] [-B] [-S] [-M freq rate] [-P snd_parallelism]\n", argv[0]);
/* >>> */
				printf ("  -?,-h : this output\n");
				printf ("  -d : debug mode\n");
				printf ("	  -e : turn on elitism strategy\n");
				printf ("  -f : specify dateset file w/ OR-LIB format\n");
				printf ("  -i : init pop creation&improve option. \n");
				printf ("	   (n)one*|(c)onstraint|(f)easibility \n");
				printf ("  -l : specify log prefix. \n");
				printf ("  -o : output dir\n");
				printf ("  -p : population size\n");
				printf ("  -r : replacement. u- unfitness; w*-worst\n");
				printf ("  -s : selection. r- rank-based; b*-binary\n");
				printf ("  -t : termination rule. improve or fixed iter\n");
				printf ("	   i N: quit if w/o improve w/ N chrom\n");
				printf ("	   f N: quit after N iterations\n");
				printf ("	   q bound: quit if bound is reached\n");
				printf ("	   t walltime: quit after walltime seconds\n");
/*
				printf ("  -c : prob of crossover [0.0 .. 1.0]\n");
				printf (" -m : prob of mutation [0.0 .. 1.0]\n");
				printf ("  -q : specify soln quality as stopping rule.\n");
				printf ("  -B : turn on blocking comm. PGA only\n");
				printf ("  -S : turn on sync comm. PGA only\n");
				printf ("  -I : immigration strategy: o-as offspring; p-as parent\n");
*/
				printf ("  -R : global seed for random num generator. PGA only\n");
				printf ("  -M : mig_freq mig_rate immig_freq. Specify migration freq and rate. PGA only\n");
				printf ("  -P : snd_parallelism. Specify the number of consecutive exports to go before sync, also determines the size of user-provided MPI sending buffer. PGA only\n");
				exit (0);
				break;
			}
		}
	}
}

void print_config() {
#ifdef PGAMODE
	if (myrank == 0) {
#endif
	printf("GA parameters:\n");
	printf("\tselection strategy:\t%c\n", strategy_selection);
	printf("\treplacement strategy:\t%c\n", strategy_replacement);
	printf("\tpopulation initilization:\t%c\n", strategy_init_pop);
	printf("\tpopulation size:\t%d\n", pop_size);
	printf("\tstopping rule:\t%c", strategy_stop);
	switch (strategy_stop) {
	case 'f':
		printf(" %d", max_iterations);
		break;
	case 'q':
		printf(" %lld", stopping_quality);
		break;
	case 't':
		printf(" %d", walltime);
		break;
	case 'i':
	default:
		printf(" %d", improve_threshold);
		break;
	}
	printf("\n");
	printf("\tdata:\t%s\n", dataset);
	printf("\toutput dir:\t%s\n", myoutputdir);
	printf("\telitism:\t%d\n", elitism);
	printf("PGA parameters:\n");
	printf("\tglobal seed:\t%ld\n", global_seed);
#ifdef PGAMODE
	printf("\texport interval:\t%d\n", migrate_freq);
	printf("\texport rate:\t%d\n", emi_size);
	printf("\timport interval:\t%d\n", immigrate_freq);
	printf("\tnumber of non-blocking exports:\t%d\n", snd_parallelism);
#endif
	fflush(stdout);
#ifdef PGAMODE
	}
#endif
}
void gastat_init() {
	// initialize stat
	gastat.total_improve = 0;
#ifdef PGAMODE
	gastat.mig_improve = 0;
#endif
	gastat.startT = get_ga_time();
	gastat.bestT = gastat.startT;
	gastat.endT = gastat.startT;
	gastat.min_fitv = -1;
	gastat.max_fitv = -1;
	gastat.avg_fitv = -1;
	gastat.commT = 0;
	gastat.comm_sendT = 0;
	gastat.comm_recvT = 0;
#ifdef T_PROFILING
	gastat.ioT = 0;
	gastat.immigrateT = 0;
	gastat.selectionT = 0;
	gastat.injectT = 0;
	gastat.crossmutT = 0;
	gastat.feasibT = 0;
	gastat.improveT = 0;
	gastat.evalT = 0;
	gastat.replaceT = 0;
#endif
}
#endif
