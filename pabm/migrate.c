// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "hdr.h"

// Macros to get access to proxy and agent information
#define GETINFO(agent,type) BEHAVIOR(type).info(agent)
#define ISPROXY(isproxy) isproxy==PROXY?1:0

// Checks for buffers
int npbndcheck(int r,int c) { return NPBNDCHECK(r,c); }
int spbndcheck(int r,int c) { return SPBNDCHECK(r,c); }
int epbndcheck(int r,int c) { return EPBNDCHECK(r,c); }
int wpbndcheck(int r,int c) { return WPBNDCHECK(r,c); }

int nbndcheck(int r,int c) { return NBNDCHECK(r,c); }
int sbndcheck(int r,int c) { return SBNDCHECK(r,c); }
int ebndcheck(int r,int c) { return EBNDCHECK(r,c); }
int wbndcheck(int r,int c) { return WBNDCHECK(r,c); }


inline int migrateproxyboundarycheck(int r,int c) {
    return PBNDCHECK(r,c); 
}

inline int migrateboundarycheck(int r,int c) {
    return BNDCHECK(r,c); 
}

// This example maintains local coordinates for each core
int globaltolocalcoord(agentstruct *agent,int type,int isproxy) {

    genericstate *gs;
    gs=OPS(type).gs(agent);

    gs->r-=env->r;
    gs->c-=env->c;

    if(ISPROXY(isproxy)) {
        if(mpi->size>1) {
            if(gs->r<-env->buf) {
                assert(mpi->r==mpi->h-1);
                gs->r+=env->r+env->h;
            }
            if(gs->r>=env->habs) {
                assert(mpi->r==0);
                gs->r-=env->hg;
            }
            if(gs->c<-env->buf) {
                assert(mpi->c==mpi->w-1);
                gs->c+=env->c+env->w;
            }
            if(gs->c>=env->wabs) {
                assert(mpi->c==0);
                gs->c-=env->wg;
            }
        } else { // boundary case (size==1)
            if(gs->r<env->buf)
                gs->r+=env->h;
            else if(gs->r>=env->h-env->buf) 
                gs->r-=env->h;
            if(gs->c<env->buf)
                gs->c+=env->w;
            else if(gs->c>=env->w-env->buf) 
                gs->c-=env->w;
        }
    }

}

int localtoglobalcoord(agentstruct *agent,int type) {
    genericstate *gs;
    gs=OPS(type).gs(agent);

    // Change to global coordinates
    gs->r+=env->r+env->hg;
    gs->c+=env->c+env->wg;
    // Wrap the coordinates if overflow
    gs->r%=env->hg;
    gs->c%=env->wg;
    assert(gs->r>=0&&gs->c>=0);
}

int proxychange(agentstruct *agent,int type) {
    //change agent to be a proxy
    genericstate *gs;
    gs=OPS(type).gs(agent);
    // Sanity checks
    if(gs->id<=0)
        printf("WRONG ID=%i alive?=%i\n",gs->id,gs->alive);
    assert(gs->id>0);
    gs->id=-gs->id;
    assert(gs->id<0);
}

int acceptagent(agentstruct *agent,int type,int isproxy) {
    // Change global coordinates to local
    globaltolocalcoord(agent,type,isproxy);

    // Accept agent
    if(ISPROXY(isproxy))
        proxychange(agent,type);
}

int buffertopop(agentstruct *abuf,agentstruct *apop,int type) {
    memcpy(apop->a,abuf->a,POP(type).sizeofagent);
}

int poptobuffer(agentstruct *apop,agentstruct *abuf,int type) {
    memcpy(abuf->a,apop->a,POP(type).sizeofagent);
    BEHAVIOR(type).del(apop); 
}


int agentemptyrecvbuffer(groupstruct *group,int isproxy) {
    agentstruct *agent;
    genericstate *gs;
    int ind,type;

    type=group->type;

    for(ind=0;ind<AGENTGROUPSIZE;ind++) {
        gs=OPS(type).gs(&group->a[ind]); // pull agent out of group
        if(gs->alive) {
            agent=agentgetnew(type);
            buffertopop(&group->a[ind],agent,type);
            acceptagent(agent,type,isproxy);
        } else
            break; // ran out of agents
    }

    return ind;
}

int agentfillsendbuffer(groupstruct *group,int (*bndcheck)()) {
    agentstruct *agent,*prev;
    genericstate *gs;
    int ind,type;

    ind=0;
    type=group->type;
    prev=agent=POP(type).migratingagents;

    //zero the buffer first
    memset(group->g,0,POP(type).sizeofagent*AGENTGROUPSIZE);

    while(ind<AGENTGROUPSIZE&&agent!=NULL) {
        gs=OPS(type).gs(agent);

        //if(bndcheck(a->g.r,a->g.c)) { // if in the boundary we are interested in
        if(bndcheck(gs->r,gs->c)) { // if in the boundary we are interested in
            //Pull agent out of linked list
            if(agent==POP(type).migratingagents)
                POP(type).migratingagents=prev=agent->next;
            else
                prev->next=agent->next;

            localtoglobalcoord(agent,type); // Before putting into buffer change coordinates to global
            poptobuffer(agent,&group->a[ind],type);

            ind++;
        } else { // if its not in the boundaries
            prev=agent; // set as the previous agent

        }
        agent=agent->next;
    }

    return ind;
}

int mpisr(int snd,int rcv,int dir,int type) {
    MPI_Status status;
    // Using MPI_BYTE only works on homogeneous clusters 

    if(snd&&rcv) { //send and recv
        MPI_Sendrecv_replace(mpi->dir[dir].groupbuf[type]->g,AGENTGROUPSIZE*POP(type).sizeofagent,MPI_BYTE,
                             mpi->dir[dir].send,700+dir,
                             mpi->dir[dir].recv,700+dir,
                             mpi->dir[dir].comm->comm,&status);
    } else if(snd) { //send
        MPI_Send(mpi->dir[dir].groupbuf[type]->g,AGENTGROUPSIZE*POP(type).sizeofagent,MPI_BYTE,
                 mpi->dir[dir].send,700+dir,
                 mpi->dir[dir].comm->comm);
    } else { // recv
        MPI_Recv(mpi->dir[dir].groupbuf[type]->g,AGENTGROUPSIZE*POP(type).sizeofagent,MPI_BYTE,
                 mpi->dir[dir].recv,700+dir,
                 mpi->dir[dir].comm->comm,&status);
    }

}

int agentmigrate(populationstruct *pop,int isproxy) {
    int type;
    int i,snd,rcv;
    int rbs,sbs; // rcv/snd buf size
    ADMIN(3,"Agents migrating");

    type=pop->type;

    //Array of function pointers for boundary checking
    int (*dirbndcheck[LASTDIR])();

    if(ISPROXY(isproxy)) {
    dirbndcheck[N]=&npbndcheck;
    dirbndcheck[S]=&spbndcheck;
    dirbndcheck[E]=&epbndcheck;
    dirbndcheck[W]=&wpbndcheck;
    } else {
    dirbndcheck[N]=&nbndcheck;
    dirbndcheck[S]=&sbndcheck;
    dirbndcheck[E]=&ebndcheck;
    dirbndcheck[W]=&wbndcheck;
    }

    // For all directions, send and receive all messages.
    for(i=0;i<LASTDIR;i++) {
        snd=1;
        rcv=1;

        while(snd||rcv) {

            if(snd) { 
                //fill buffer
                sbs=agentfillsendbuffer(mpi->dir[i].groupbuf[type],dirbndcheck[i]);
            }

            mpisr(snd,rcv,i,type);

            if(rcv) { // recv
                //accept agents
                rbs=agentemptyrecvbuffer(mpi->dir[i].groupbuf[type],isproxy);
            }
            if(sbs<AGENTGROUPSIZE)  // If snd buffer is not full
                snd=0;              // No more agents to send
            if(rbs<AGENTGROUPSIZE)  // If rcv buffer is not full
                rcv=0;              // No more agents to recv
        }
    }
}

// Identify agent proxies that must be migrated
int agentidentifymigratingproxies(agentstruct *agent,int type) {
    genericstate *gs;
    int r,c;

    gs=OPS(type).gs(agent);
    if(gs->alive) {
        if(migrateproxyboundarycheck(gs->r,gs->c)) {
            r=gs->r;
            c=gs->c;
            assert(r>=-env->buf&&c>=-env->buf&&r<env->h+env->buf&&c<env->w+env->buf);

            agent->next=POP(type).migratingagents;
            POP(type).migratingagents=agent;
        }
    }

}

// See if alive agent is in local environment
int agentinenv(agentstruct *agent,int type) {
    genericstate *gs;

    gs=OPS(type).gs(agent);
    if(gs->alive) {
        if(!migrateboundarycheck(gs->r,gs->c)) {
            return 1;
        }
    }
    return 0;
}

int agentidentifymigrating(agentstruct *agent,int type) {
    genericstate *gs;
    gs=OPS(type).gs(agent);
    if(gs->alive) { // If agent is alive and outside boundaries
        if(migrateboundarycheck(gs->r,gs->c)) {
            agent->next=POP(type).migratingagents; // Migrate
            POP(type).migratingagents=agent;
            return 1;
        }
    }
    return 0;
}


int popmigrate(populationstruct *pop) {

    EXEC0(DEBUG(2,printf(" [ Agent transfer ]\n")));
    //Make sure no agents are migrating from another iteration
    assert(pop->migratingagents==NULL);

    // Find all agents that are migrating
    gaiter(pop->pop,&agentidentifymigrating);

    // Migrate the agents 
    agentmigrate(pop,AGENT); 

    assert(pop->migratingagents==NULL);

}

int popmigrateproxies(populationstruct *pop) {

    EXEC0(DEBUG(2,printf(" [ Agent proxy transfer ]\n")));

    //Make sure no agents are migrating from another iteration
    assert(pop->migratingagents==NULL);

    // Find all agent proxies that are migrating
    gaiter(pop->pop,&agentidentifymigratingproxies);

    //printf("migratingproxies=%i\n",agentcount(pop->migratingagents));

    agentmigrate(pop,PROXY); 

    assert(pop->migratingagents==NULL);

}


