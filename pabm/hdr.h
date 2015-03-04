// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
//#define NDEBUG
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <unistd.h>

#include "include/util.h"
#include "include/sugar.h"

#define DEBUGLVL 2 

// DEFINITIONS // DEFINITIONS // DEFINITIONS // DEFINITIONS // DEFINITIONS // DEFINITIONS // DEFINITIONS // DEFINITIONS

// Enable or disable the load balancer
#define LB 
#define LB_FREQ 50 


#define NUMROWS 128

// Set number of columns equal to the number of rows
#define NUMCOLS NUMROWS

// Set block size for loadbalancer
#define BLOCK 4 

#define NUMAGENTS 1024
#define ITERATIONS 100

#define AGENTGROUPSIZE 128 

#ifdef LB
#define LOADBALANCE 1
#else
#define LOADBALANCE 0
#endif

#define RMULT 1 // Row multiplier
#define CMULT 3 // Col multiplier



#define INFLUENCEDIST 2


// OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS // OPERATORS 

typedef struct {
    int r,c;
} locstruct;


typedef struct {
    int (*update)();
    int (*init)();
    int (*del)();
    int (*move)();
    int (*info)(); // Function which prints the info of the agent
    int (*interact)(); // Interact agent1 with agent2 
} behaviorstruct;

typedef struct {
    void * (*allocgroup)(); // Alloc a group of agents of this type
    void * (*agentingroup)(); // Get an agent from a group
    int (*init)(); // Initialize the population params and setup simulation for this type
    int (*del)(); // Del any initialized population params for the simulation
    int (*alive)(); // Function which returns whether or not agent is alive 
    locstruct (*loc)(); // Function which returns location of agent 
    genericstate * (*gs)(); // Function which returns a pointer to the generic state 
} operatorstruct;


// AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS // AGENTS

enum AGENTTYPE {EXAMPLETYPE,LASTAGENT};

struct aagentstruct {
    void *a;
    struct aagentstruct *next;
};
typedef struct aagentstruct agentstruct;

struct ggroupstruct {
    void *g;
    int type;
    agentstruct a[AGENTGROUPSIZE];
    struct ggroupstruct *next;
};
typedef struct ggroupstruct  groupstruct;

typedef struct {
    // Population
    groupstruct *pop; // Population of agents
    int population;  // The actual population
    int type;      // The type of agents in this population
    int sizeofagent; // Size of an agent

    // A list of agents that will be migrating this iteration
    agentstruct *migratingagents;

    // Time saving structures (free agents and groups)
    groupstruct *freegroups;
    agentstruct *freeagents;

    // Operators and behaviors for this population
    operatorstruct ops;
    behaviorstruct behavior; 
} populationstruct;

typedef struct {
    populationstruct pop[LASTAGENT];
} agentsstruct;


// ENVIRONMENT // ENVIRONMENT // ENVIRONMENT // ENVIRONMENT // ENVIRONMENT // ENVIRONMENT // ENVIRONMENT // ENVIRONMENT 

typedef struct {
    int r,c; // row and column of tile
    int h,w; // height and width of tile
} tilestruct;

typedef struct {
    float val;
    float cap;
} cellstruct;
#define CELLSTRUCTELEMENTS 2+1

typedef float lbtype; // Could be an int

typedef struct {
    lbtype c; // cost for computation
    lbtype n,s,e,w; // cost for communication
} blockstruct;

typedef struct {
    groupstruct *b;
    MPI_Status s;   
    MPI_Request r;
    int proc;
} bufstruct;

typedef struct {
    int blocksize;  // block size for grain coarsening
    int hb,wb;      // h,w for block array
    int dir;        // direction for LB (N/S)
    int pairoffset; // determines the pairs for distributed LB calculation (odd always sends, po sets E/W direction)
    int lind,rind;  // save left and right partition index

    blockstruct **blocks; // blocks

    blockstruct *superblocks[2]; // superblock aggregations in each direction

    blockstruct *superblockbuf; // superblock buffer
    blockstruct *superblockbuf2; // superblock buffer for neighboring processor
    blockstruct **superblockarray; // superblock array
    int sbasize;
} lbstruct;


typedef struct {
    int r,c; // row and column in parallel environment
    int h,w; // height and width of this parallel abm
    int habs,wabs; // absolute height and width (includes buffer on all sides)
    int hmax,wmax; // maximum height and width (evenly distributed env * row/col multiplier)
    int hg,wg; // global height and width
    int buf; // buffer size (ghost zone)
    cellstruct **grid,**gridabs,**gridmax;
    lbstruct lb;
} envstruct;


// MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI // MPI

enum ENUMDIR {N,S,E,W,LASTDIR};

typedef struct {
    int rank,size; // rank and size of comm
    MPI_Comm comm; // communicator

} commstruct;

typedef struct {
    int send,recv; // send processor and recv processor
    tilestruct tilesend,tilerecv; // tile send and recv 
    cellstruct *envbuf;
    groupstruct *groupbuf[LASTAGENT];
    commstruct *comm;
} dirstruct;

typedef struct {
    int r,c; // row and column in processor grid
    int h,w; // maximum row and column in processor grid organization
    int rank,size; // rank and size in MPI_COMM_WORLD
    commstruct cart,row,col;
    dirstruct dir[LASTDIR];
    MPI_Datatype MPI_CELL,MPI_AGENT,MPI_GROUP;
} mpistruct;



// SIMULATION // SIMULATION // SIMULATION // SIMULATION // SIMULATION // SIMULATION // SIMULATION // SIMULATION // SIMULATION 
typedef struct {
    int as[LASTDIR],ar[LASTDIR]; // agents sent and recv'ed
    double iterstart,iterend;
} statstruct;

typedef struct {
    int iteration; // current iteration
    int maxiter; // number of iterations in this simulation
    int seed; // random seed
    int argc; // argc from command line
    char **argv; // argv from command line
} simstruct;


// MACROS // MACROS // MACROS // MACROS // MACROS // MACROS // MACROS // MACROS // MACROS // MACROS // MACROS // MACROS 

enum ENUMRC {R,C};
enum ENUMAP {AGENT,PROXY}; // proxy should be 1
enum ENUMLR {LOCAL,REMOTE};
#define ALIVE(agent) OPS(type).alive((agent))
#define max(a,b) (a>b)?a:b
#define min(a,b) (a<b)?a:b
#define distance(ar,ac,br,bc) sqrt((ar-br)*(ar-br) + (ac-bc)*(ac-bc))

//shortcuts to agents functions
#define POP(type) agents->pop[type]
#define BEHAVIOR(type) agents->pop[type].behavior
#define OPS(type) agents->pop[type].ops


#define DEBUG(lvl,funct) if(DEBUGLVL>=lvl) funct

// GLOBALS // GLOBALS // GLOBALS // GLOBALS // GLOBALS // GLOBALS // GLOBALS // GLOBALS // GLOBALS // GLOBALS // GLOBALS 

extern mpistruct *mpi;
extern envstruct *env;
extern agentsstruct *agents;
extern simstruct *sim;


// FUNCTIONS // FUNCTIONS // FUNCTIONS // FUNCTIONS // FUNCTIONS // FUNCTIONS // FUNCTIONS // FUNCTIONS // FUNCTIONS

//int simulationmosqinit(populationstruct *pop);
int simulationsuginit(populationstruct *pop);
groupstruct *groupalloc(int type);
agentstruct* agentgetnew(int type);
int envprint();
cellstruct *envgetrow();
tilestruct getnewtile(int r,int c);
int agentidentifymigrating(agentstruct *agent,int type);
int agentlink(agentstruct *agent,int type); 
agentstruct * agentgetintile(tilestruct t,int type);

lbtype partition2d(lbtype **ps,lbtype **comm,int height,int size,int cuts);
int probe(lbtype *ps,int size,lbtype wt,int cuts); 
int getpartitions(int *partitions,int numparts,lbtype **arr,int height,int size,lbtype wt); 
int search(lbtype *arr,int size,lbtype w,int s);
int search2d(lbtype **arr,int height,int size,lbtype wt,int s); 

lbtype part2d(lbtype **ps,lbtype **comm,int height,int size,int cuts);
int partition_find2d(lbtype **ps,int height,int size,int s,int cuts);
lbtype maxweight2d(lbtype **ps,int height,int size,int s,int e);
int probe2d(lbtype **arr,int height,int size,lbtype wt,int cuts);


int popmigrate(populationstruct *pop);

