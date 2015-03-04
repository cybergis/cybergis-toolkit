// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "hdr.h"

// Get a unique ID
int globalid;
#define GETID ++globalid
#define SETID(id) globalid=id;

// Delete a sugar agent
int sugdel(agentstruct *a) {
    xplstruct *m;
    m=(xplstruct *)a->a;
    
    m->g.id=-1; 
    m->g.alive=0; 
    return 0;
}

// Print information about agent
int suginfo(agentstruct *a) {
    xplstruct *m;
    m=(xplstruct *)a->a;

    printf("SUG id=%i astate value=%i ABS(%i,%i)\n",m->g.id,m->s.astate,env->r+m->g.r,env->c+m->g.c);
}

int weightedrand(int N) {
    int i;
    int sum,cur,rnd;

    sum=(1+N)*(N/2);

    rnd=rand()%sum;
    cur=0;
    for(i=0;i<N;i++) {
        cur+=i;
        if(cur>rnd)
                break;
    }
    return i; // more likely to be close to N than 0 (linear weight)

}

int sugsetrc(agentstruct *a,int r,int c) {
    xplstruct *m;
    m=(xplstruct *)a->a;
    m->g.r=r;
    m->g.c=c;
}

// Initialize agent
int suginit(agentstruct *a) {
    xplstruct *m;
    m=(xplstruct *)a->a;
    
    m->g.alive=1; // Is alive
    m->g.id=GETID; //  Get a unique ID
    m->s.astate=10; // Give an initial astate value for this example

    //Randomly assign a location
    m->g.r=(int)((float)drand48()*(float)(env->h-1));
    m->g.c=(int)((float)drand48()*(float)(env->w-1));

    // Sanity check
    assert(m->g.id>0);
    assert(m->g.r>=0&&m->g.r<env->h);
    assert(m->g.c>=0&&m->g.c<env->w);

    return 0;
}

int _moved(cellstruct **grid,int *ro,int *co) {
    int i,j;
    int r,c;
    int maxr,maxc;
    float maxsugar;

    r=*ro;
    c=*co;

    maxsugar=FLT_MIN;
    //Look over all the rows
    for(i=r-INFLUENCEDIST;i<=r+INFLUENCEDIST;i++)
        if(grid[i][c].val>maxsugar) {
            maxsugar=grid[i][c].val;
            maxr=i;
            maxc=c;
        }

    //Look over all the cols 
    for(j=c-INFLUENCEDIST;j<=c+INFLUENCEDIST;j++)
        if(grid[r][j].val>maxsugar) {
            maxsugar=grid[r][j].val;
            maxr=r;
            maxc=j;
        }

    assert(maxr==r||maxc==c); // only move in 1 direction or the other, but not both

    // Move to maximum value
    *ro=maxr;
    *co=maxc;

    return 0;
}

int sugmove(agentstruct *a) {
    int r,c;
    int buf;
    cellstruct **grid;
    int origr,origc;
    int maxr,maxc;
    int i,j;

    xplstruct *m;
    m=(xplstruct *)a->a;

    // Skip proxies
    if(m->g.id<0)
        return 0;

    // Original row and column
    origr=m->g.r;
    origc=m->g.c;

    //when working with agent locations in this area use absolute grid location (not regular grid)
    grid=env->gridabs;
    buf=env->buf;

    // Exact grid
    r=m->g.r+buf;
    c=m->g.c+buf;

    // Move to maximum value cell in nearby area
    _moved(grid,&r,&c);

    // Sanity check
    assert(abs((r-buf)-origr)<=INFLUENCEDIST);
    assert(abs((c-buf)-origc)<=INFLUENCEDIST);

    // Assign location
    m->g.r=r-buf;
    m->g.c=c-buf;

    return 0;
}

// Interact with another agent
int suginteract(agentstruct *a1,agentstruct *a2) {
    xplstruct *m2;
    m2=(xplstruct *)a2->a;
    m2->g.alive++;
}

// Update an agent
int sugupdate(agentstruct *a) {
    xplstruct *m;
    m=(xplstruct *)a->a;
    
    // In this example agents extract 10% of the environment resources
    m->s.astate+=env->grid[m->g.r][m->g.c].val*0.1;
    env->grid[m->g.r][m->g.c].val-=env->grid[m->g.r][m->g.c].val*0.1;

    return 0;
}

// Check if agent is alive
int sugalive(agentstruct *a) {
    xplstruct *m;
    m=(xplstruct *)a->a;
    return m->g.alive;
}

// Extract genericstate from agent
genericstate *suggs(agentstruct *a) {
    xplstruct *m;
    m=(xplstruct *)a->a;
    return &m->g;
}

// Extract location from agent
locstruct sugloc(agentstruct *a) {
    xplstruct *m;
    locstruct l;
    m=(xplstruct *)a->a;
    l.r=m->g.r;
    l.c=m->g.c;
    return l;
}

void *sugallocgroup() {
    void *g;

    g=calloc(AGENTGROUPSIZE,sizeof(xplstruct));
    return g;
}

void *sugagentingroup(void *g,int i) {
    xplstruct *group;
    group=(xplstruct *)g;
    return (void *)&group[i];
}

int simulationsugdel(populationstruct *pop) { 
    //free(pop->behavior);
}


// Initialize agent population
int simulationsuginit(populationstruct *pop) {
    int i,j;
    agentstruct *a;
    char fn[50];
 
    // This example assumes few agents
    SETID(mpi->rank*10000000);

    // Initialize operations for this population 
    pop->ops.del=&simulationsugdel;
    pop->ops.allocgroup=&sugallocgroup;
    pop->ops.agentingroup=&sugagentingroup;
    pop->ops.alive=&sugalive;
    pop->ops.loc=&sugloc;
    pop->ops.gs=&suggs;

    // Initialize population behavior 
    pop->behavior.init=&suginit;
    pop->behavior.del=&sugdel;
    pop->behavior.move=&sugmove;
    pop->behavior.update=&sugupdate;
    pop->behavior.info=&suginfo;
    pop->behavior.interact=&suginteract;

    // Initialize values
    pop->sizeofagent=sizeof(xplstruct);
    pop->population=NUMAGENTS + (NUMAGENTS*.8*mpi->c);

    // Initialize the agents in the initial population popsize
    for(i=0;i<pop->population;i++) {
        a=agentgetnew(EXAMPLETYPE);
        suginit(a);

    }

    return 0;
}

