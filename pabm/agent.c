// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "hdr.h"

// AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT AGENT //

#define ISLOCAL(islocal) islocal==LOCAL?1:0


int agentdelproxy(agentstruct *agent,int type) {
    genericstate *a;
    a=(genericstate *)agent->a;
    if(a->id<0)
        agents->pop[type].behavior.del(agent);
    return 0; 
}


// Call the specific population update functions
inline int agentupdate(agentstruct *a,int type) {
        
    if(ALIVE(a)) {
        BEHAVIOR(type).update(a);
        return 0;
    }
    return 1;
}

// Call the specific population move functions
inline int agentmove(agentstruct *a,int type) {
        
    if(ALIVE(a)) {
        BEHAVIOR(type).move(a);
        return 0;
    }
    return 1;
}

// Get a new agent from the population list (reuse or alloc)
agentstruct *agentgetnew(int type) {
    agentstruct *freeagent;

    if(POP(type).freeagents==NULL) {
        //alloc new space and link new free agents
        groupgetnew(type);
    }
    freeagent=POP(type).freeagents;
    POP(type).freeagents=freeagent->next;
    assert(freeagent!=NULL);
    return freeagent;
}   

// Apply function to the list of agents a
void agentiter(agentstruct *a,int (*function)()) {
    agentstruct *an;
    while(a!=NULL) {
        an=a->next;
        (*function)(a); 
        a=an;
    }
}

// Shift agent (genericstate) by r,c
int gsshiftloc(genericstate *gs,int r,int c) {
    gs->r+=r;
    gs->c+=c;
}

// Shift agent (agentstruct) by r,c
int agentshiftloc(agentstruct *agent,int type,int r,int c) {
    genericstate *gs;
    gs=OPS(type).gs(agent);
    gsshiftloc(gs,r,c);
}

// Count agents following linked list
int agentcount(agentstruct *a) {
    int sum=0;
    while(a!=NULL) {
        sum++;
        a=a->next;
    }
    return sum;
}

// If agent is in tile, add to list
agentstruct * agentgetintile(tilestruct t,int type) {
    int i;
    agentstruct *a;
    genericstate *gs;
    groupstruct *g;
    int r,c;


    g=POP(type).pop;

    a=NULL;

    while(g!=NULL) {
        for(i=0;i<AGENTGROUPSIZE;i++) {
            if(OPS(type).alive(&g->a[i])) {
                gs=OPS(type).gs(&g->a[i]);
                r=gs->r;
                c=gs->c;
                if(r>=t.r&&r<t.r+t.h&&c>=t.c&&c<t.c+t.w) {
                    g->a[i].next=a;
                    a=&g->a[i];
                }
            }
        }
        g=g->next;
    }
    return a;
}

int agentcellaaint(agentstruct *agent,cellstruct *cell,int type) {
/*
    // find all agents in cell and call aa interaction with passed in agent
    agentstruct *aint;
    aint=cell->a;
    while(aint!=NULL) {
        BEHAVIOR(type).interact(agent,aint);
        aint=aint->next;
    }
*/
}

/*
int agentaaint(agentstruct *agent,int type,int islocal) {
    int i,r,c;
    genericagentstruct *g;

    g=(genericagentstruct *)agent->a;
    r=g->g.r;
    c=g->g.c;
    if(ISLOCAL(islocal)) {
        // only interact with local agents
        for(i=0;i<INFLUENCEDIST;i++) {
            if(!BNDCHECK(r,c+i))
                agentcellaaint(agent,&env->grid[r][c+i],type);
            if(!BNDCHECK(r,c-i))
                agentcellaaint(agent,&env->grid[r][c-i],type);
            if(!BNDCHECK(r+i,c))
                agentcellaaint(agent,&env->grid[r+i][c],type);
            if(!BNDCHECK(r-i,c))
                agentcellaaint(agent,&env->grid[r-i][c],type);
        }
    } else { // is remote
        // only interact with remote agents
        for(i=0;i<INFLUENCEDIST;i++) {
            if(BNDCHECK(r,c+i))
                agentcellaaint(agent,&env->grid[r][c+i],type);
            if(BNDCHECK(r,c-i))
                agentcellaaint(agent,&env->grid[r][c-i],type);
            if(BNDCHECK(r+i,c))
                agentcellaaint(agent,&env->grid[r+i][c],type);
            if(BNDCHECK(r-i,c))
                agentcellaaint(agent,&env->grid[r-i][c],type);
        }
    }
}
*/

int cellaaint(cellstruct *cell,int type,int islocal) {
/*
    agentstruct *a;
    // apply agentaaint to all agents in this cell
    a=cell->a;
    while(a!=NULL) {
        agentaaint(a,type,islocal); 
        a=a->next;
    }
*/
}


// GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP GROUP //

// Allocate a generic group
groupstruct *groupgenalloc(int type) {
    groupstruct *group;

    group=(groupstruct *)calloc(1,sizeof(groupstruct));
    group->type=type;

    return group;
}

// Handle extract allocations for this population
void *groupspecificalloc(int type) {
    void *g;

    g=NULL;
    g=OPS(type).allocgroup();
    if(g==NULL)
        DIE("problem with allocgroupspecific");

    return g;
}

// Allocate a new group structure
groupstruct *groupalloc(int type) {
    agentstruct *agent;
    groupstruct *group;
    int i;

    // Allocate group structures
    group=groupgenalloc(type);
    group->g=groupspecificalloc(type); 
    assert(group!=NULL);
    assert(group->g!=NULL); 
    // Setup agents in group as available
    for(i=0;i<AGENTGROUPSIZE;i++) {
        agent=&group->a[i];
        agent->a=OPS(type).agentingroup(group->g,i);
        BEHAVIOR(type).del(agent);
    }
    return group;
}

// Allocate a new group and add to population list 
int groupgetnew(int type) {
    groupstruct *group;
    agentstruct *agent;
    int i;

    group=groupalloc(type);

    group->next=POP(type).pop; 
    POP(type).pop=group; 
    for(i=0;i<AGENTGROUPSIZE;i++) {
        agent=&group->a[i];
        agent->next=POP(type).freeagents;
        POP(type).freeagents=agent;
        BEHAVIOR(type).del(agent);
    }
}

int groupspecificfree(void *g) {
    free(g);
    return 0;
}

int groupfree(groupstruct *group) {
    groupspecificfree(group->g);
    free(group);
    return 0;
}

int groupdel(groupstruct *g) {
    groupfree(g); // Free the group
    return 0;
}

// Apply function to group
void groupiter(groupstruct *g,int (*function)()) {
    groupstruct *gn;
    while(g!=NULL) {
        gn=g->next;
        (*function)(g);
        g=gn;
    }
}

// Apply function to agents in group
void gaiter(groupstruct *g,int (*function)()) {
    int i,type;
    type=g->type;
    while(g!=NULL) {
        for(i=0;i<AGENTGROUPSIZE;i++) {
            (*function)(&g->a[i],type);
        }
        g=g->next;
    }
}

// POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP POP //

// Shift agents in population by r,c
int popshiftloc(populationstruct *pop,int r,int c) {
    groupstruct *g;
    int i;
    int type;

    type=pop->type;
    g=pop->pop;
    while(g!=NULL) {
        for(i=0;i<AGENTGROUPSIZE;i++) {
            if(ALIVE((&g->a[i])))
                agentshiftloc(&g->a[i],type,r,c);
        }
        g=g->next;
    }

}

int agentaaint(agentstruct *agent,populationstruct *pop) {
    agentstruct *aint;
    groupstruct *g;
    genericstate *gs;
    int i;
    int type;

    g=pop->pop;
    while(g!=NULL) {
        for(i=0;i<AGENTGROUPSIZE;i++) {
            if(ALIVE((&g->a[i]))) {
                aint=&g->a[i];
                gs=(genericstate*)aint;
                if(aint==agent) continue; // skip self
                // Agent-agent interaction has been removed in this example
                // Agent proxies are available for agent-agent interactions
                // Proxies are indicated with a negative ID, which corresponds
                // to the positive ID on the remote core, which is identified by
                // the ghost zone region on which the proxy resides.
                if(gs->id<0) { // It is a proxy
                    // Add interaction message to send queue
                    // Once the message is sent to the appropriate core
                    // A receive operation will initiate this interaction behavior
                    // for each interaction message. Since it will be interacting
                    // with a local agent, it will enter the else.  If a response
                    // is generated it will follow this logic again by responding
                    // to the proxy, which will again invoke an interaction message
                    // to be added to the send queue, sent and interacted on the other core
                } else {
                    // local agent
                    // In this example all agents interact with all other local agents
                    // Models are likley to restrict interaction to nearby distance
                    // or shared characteristics.
                    BEHAVIOR(type).interact(agent,aint);
                }
            }
        }
        g=g->next;
    }
}

int popaaint(populationstruct *pop) {
    int i,j,type;
    agentstruct *a;
    groupstruct *g;
    genericstate *gs;

    type=pop->type;

    EXEC0(DEBUG(2,printf(" [ Agent-agent interaction ]\n")));


    g=pop->pop;
    while(g!=NULL) {
        for(i=0;i<AGENTGROUPSIZE;i++) {
            // For all agents that are alive and not a proxy
            if(ALIVE((&g->a[i]))) {
                gs=(genericstate *)g->a[i].a;
                if(gs->id>0)
                    agentaaint(&g->a[i],pop); // Initiate agent-agent interaction
            }
        }
        g=g->next;
    }

}

// Count the population
int popcount(populationstruct *pop) {
    int i,type;
    groupstruct *group;
    int count;

    count=0;
    type=pop->type;
    group=pop->pop;
    while(group!=NULL) {
        for(i=0;i<AGENTGROUPSIZE;i++) {
            if(OPS(type).alive(&group->a[i]))
                count++;
        }
        group=group->next;
    }
    return count;
}

// Update the population
int popupdate(populationstruct *pop) {

    EXEC0(DEBUG(2,printf(" [ Agent move ]\n")));
    gaiter(pop->pop,&agentmove);        // move

    popmigrate(pop);                    // migrate

    EXEC0(DEBUG(2,printf(" [ Agent update ]\n")));
    gaiter(pop->pop,&agentupdate);      // update

    // To enable agent-agent interactions
    // proxies must be migrated
    // then agent-agent interaction is complete
    // where proxies are triggered to communicate
    // finally the proxies should be deleted
    popmigrateproxies(pop);             // migrate proxies 

    popaaint(pop);                      // agent-agent interactions

    gaiter(pop->pop,agentdelproxy);     // delete proxies
}

// Initialize population variables
int popinit(populationstruct *pop,int type) {
    pop->pop=NULL;
    pop->migratingagents=NULL;
    pop->freeagents=NULL;
    pop->freegroups=NULL;
    pop->population=0;
    pop->type=type;
    pop->ops.init(pop);
}

// Clear the population
int popdel(populationstruct *pop) {

    //Specific population delete instruction
    pop->ops.del(); 

    groupiter(pop->pop,&groupdel);
    groupiter(pop->freegroups,&groupdel);
}



// AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS AGENTS //

int agentsshiftloc(int r,int c) {
    int i;
    //printf("agentsshiftloc()\n");
    for(i=0;i<LASTAGENT;i++) {
        popshiftloc(&agents->pop[i],r,c);
    }
}

int agentsiter(int (*funct)()) {
    int i;
    for(i=0;i<LASTAGENT;i++)
        gaiter(POP(i).pop,funct);
}

int agentscount() {
    int i,count;

    count=0;
    for(i=0;i<LASTAGENT;i++) {
        count+=popcount(&agents->pop[i]);
    }
    return count;
}

int agentscountglobal() {
    int count,globalcount;
    count=agentscount();
    MPI_Allreduce(&count,&globalcount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    return globalcount;
}

int agentsupdate() {
    int i;
    ADMIN(1,"Agents update");

    // Iterate over all populations
    for(i=0;i<LASTAGENT;i++) {
        popupdate(&agents->pop[i]);
    }
}

// Initialize all agents
int agentsinit() {
    int i,j;
    int gcount;

    DEBUG(2,printf("Agent init\n"););

    // Agent population initialization function
    OPS(0).init=&simulationsuginit;

    // Initialize all agent populations (this example only has one population)
    for(i=0;i<LASTAGENT;i++) {
        popinit(&agents->pop[i],i);
    }

    // Allocate group buffers for all agent populations
    for(i=0;i<LASTDIR;i++) 
        for(j=0;j<LASTAGENT;j++) 
            mpi->dir[i].groupbuf[j]=groupalloc(j); 

    // Print agent statistics
    EXEC0(printf("Agents are unevenly distributed amongst processor cores for demonstration\n"));

    gcount=agentscountglobal();
    EXEC0(printf(" Agents (global): %i\n",gcount));
    printf("  Core[%i] Number of agents : %i\n",mpi->rank,agentscount());
}

// Clear all agents
int agentsdel() {
    int i,j;
    ADMIN(1,"Agents del");

    // MPI related
    // These are the agent buffers for MPI, but they couldn't be created, because the agent functions hadn't been initialized yet
    for(i=0;i<LASTDIR;i++)
        for(j=0;j<LASTAGENT;j++) 
            groupfree(mpi->dir[i].groupbuf[j]);


    for(i=0;i<LASTAGENT;i++) {
        popdel(&agents->pop[i]);
    }
}

