// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "hdr.h"

// In this example agent and environment weights for communication and computation have not been calculated
// The computation and communication weights are dependent on a given model (e.g. interactions, movements, etc)
// Changing the relative weights will influence the load-balancing strategy

#define AGENTCOMPWEIGHT 20 // Weight for computing agent
#define AGENTCOMMWEIGHT 20 // Weight for communicating agent (account for ghost zones/proxies)
#define ENVCOMPWEIGHT    1 // Weight for computing environment
#define ENVCOMMWEIGHT    1 // Weight for communication env  (account for ghost zones/proxies)

// Set the tile values
tilestruct tileset(int r,int c,int h,int w) {
    tilestruct t;
    t.r=r;
    t.c=c;
    t.h=h;
    t.w=w;
    return t;
}

// Shift the environment right
int shiftenvright(int start,int size,int offset) { // shift all rows (start-start+size-1) right by offset
    int i,j;

    assert(start+size+offset<env->wmax-env->buf*2); // don't overflow
    //Note: must work backwards, otherwise will overwrite data which is bad

    for(i=0+env->buf;i<env->h;i++) { // all rows (add buf b/c in maxolute grid)
        for(j=start+size-1+env->buf;j>=start;j--) { // work right to left (add buf b/c in maxolute grid)
            env->gridmax[i][j+offset]=env->gridmax[i][j]; // move the cell at i,j to i,j+offset (i.e. shift right by offset)
        }
    }
}

// Shift the environment left
int shiftenvleft(int start,int size,int offset) { // shift all rows (start-start+size-1) left by offset
    int i,j;

    assert(start-0>=0); // don't overflow

    for(i=0;i<env->h;i++) { // all rows
        for(j=start;j<size;j++) { // work from left to right
            env->grid[i][j-offset]=env->grid[i][j]; // move the cell at i,j to i,j-offset (i.e. shift left by offset)
        }
    }
}

int sendagtdir(groupstruct *buf,agentstruct *list,int type,int dir) { // list is a linked list of agents to be sent
    agentstruct *a,*n;
    int ind;

    //zero out buffer
    memset(buf->g,0,POP(type).sizeofagent*AGENTGROUPSIZE);

    ind=0;
    a=list;
    while(a!=NULL) { // while agents in the list
        n=a->next;
        localtoglobalcoord(a,type); // switch to global coords
        poptobuffer(a,&buf->a[ind],type); // move to buffer
        ind++;
        if(ind==AGENTGROUPSIZE) { // send full group
            // Send buf
            MPI_Send(buf->g,AGENTGROUPSIZE*POP(type).sizeofagent,MPI_BYTE,
                     mpi->dir[dir].send,900,mpi->dir[dir].comm->comm);
            ind=0;
            memset(buf->g,0,POP(type).sizeofagent*AGENTGROUPSIZE); //zero out buffer

        }
        a=n;
    }
    //Last send to empty out partial group
    MPI_Send(buf->g,AGENTGROUPSIZE*POP(type).sizeofagent,MPI_BYTE,
             mpi->dir[dir].send,900,mpi->dir[dir].comm->comm);
    
}

int sendagtleft(groupstruct *buf,int width) { 
    int dir;
    agentstruct *list;
    int type=0;

    tilestruct tile;
    dir=W;

    tile.r=0;
    tile.h=env->h;
    tile.c=0;
    tile.w=width;

    list=agentgetintile(tile,type);

    sendagtdir(buf,list,type,dir);
    sendagtdir(buf,list,type,dir);
}

int sendagtright(groupstruct *buf,int width) { // buf is array of agents size AGENTGROUPSIZE 
    int dir;
    agentstruct *list;
    int type=0;

    tilestruct tile;
    dir=E;

    tile.r=0;
    tile.h=env->h;
    tile.c=env->w-width;
    tile.w=width;

    list=agentgetintile(tile,type);

    sendagtdir(buf,list,type,dir);
}


int recvagtdir(groupstruct *buf,int dir) { // buf is array of agents size AGENTGROUPSIZE
    int moreagents,count;
    MPI_Status status;
    int type=0;

    moreagents=1;

    // Loop over receiving agents, until no more are sent.
    while(moreagents) {
        // Recv buf
        MPI_Recv(buf->g,AGENTGROUPSIZE*POP(type).sizeofagent,MPI_BYTE,
                 mpi->dir[dir].recv,900,mpi->dir[dir].comm->comm,&status);

        count=agentemptyrecvbuffer(buf,0);
        if(count<AGENTGROUPSIZE)
            moreagents=0;
    }
}

// Receive agents from left
int recvagtleft(groupstruct *buf) {
    recvagtdir(buf,E);
}

// Receive agents from right 
int recvagtright(groupstruct *buf) {
    recvagtdir(buf,W);
}


int sendenvleft(cellstruct *buf,int size,int width) {
    int i,j;
    int bs;
    int iter;
    tilestruct tile;

    bs=env->lb.blocksize;

    for(iter=0;iter<size/width;iter++) {
        // Fill buf
        tile.r=0;
        tile.h=env->h;
        tile.c=iter*width;
        tile.w=width;
        blocktobuffer(env->grid,buf,&tile);

        // Send buf
        MPI_Send(buf,env->h*bs*width,mpi->MPI_CELL,
                 mpi->dir[W].send,500,mpi->dir[W].comm->comm);
    }
}

int sendenvright(cellstruct *buf,int size,int width) {
    int i,j;
    int bs;
    int iter;
    tilestruct tile;

    bs=env->lb.blocksize;

    for(iter=0;iter<size/width;iter++) {
        // Fill buf
        tile.r=0;
        tile.h=env->h;
        tile.c=env->w-size+iter*width;
        tile.w=width;
        blocktobuffer(env->grid,buf,&tile);

        // Send buf
        MPI_Send(buf,env->h*bs*width,mpi->MPI_CELL,
                 mpi->dir[E].send,500,mpi->dir[E].comm->comm);
    }
}

int recvenvleft(cellstruct *buf,int size,int width) {
    int i,j;
    int bs;
    int iter;
    MPI_Status status;
    tilestruct tile;

    bs=env->lb.blocksize;

    for(iter=0;iter<size/width;iter++) {
        // Recv buf
        MPI_Recv(buf,env->h*bs*width,mpi->MPI_CELL,
                 mpi->dir[E].recv,500,mpi->dir[E].comm->comm,&status);

        // Empty buf
        tile.r=0;
        tile.h=env->h;
        tile.c=iter*width;
        tile.w=width;
        buffertoblock(buf,env->grid,&tile);
    }

}

int recvenvright(cellstruct *buf,int size,int width) {
    int i,j;
    int bs;
    int iter;
    MPI_Status status;
    tilestruct tile;

    bs=env->lb.blocksize;

    for(iter=0;iter<size/width;iter++) {
        // Recv buf
        MPI_Recv(buf,env->h*bs*width,mpi->MPI_CELL,
                 mpi->dir[W].recv,500,mpi->dir[W].comm->comm,&status);

        // Empty buf
        tile.r=0;
        tile.h=env->h;
        tile.c=env->w+iter*width;
        tile.w=width;
        buffertoblock(buf,env->grid,&tile);
    }

}


int shiftenvandagents() {
    int bs,width;
    cellstruct *buf;
    groupstruct *grpbuf;
    int size;

    bs=env->lb.blocksize;
    width=1;
    buf=(cellstruct *)calloc(env->h*bs*width,sizeof(cellstruct));
    grpbuf=groupalloc(0);

    if(!edgeprocessorinexchange()) {
        if(mpi->c%2) { // odd

            if(env->c/bs==env->lb.lind) { // then no change
            } else if(env->c/bs<env->lb.lind) {
                size=env->lb.lind*bs-env->c;
                sendenvleft(buf,size,width); // send buffer left
                shiftenvleft(size,env->w-size,size); // shift everything left
            } else {
                size=env->c-env->lb.lind*bs;
                shiftenvright(0,env->w,size); // shift everything right
                recvenvleft(buf,size,width); // the receive from the right
            }                

        } else { // even
            if((env->c+env->w)/bs==env->lb.rind) { // then no change
            } else if((env->c+env->w)/bs<env->lb.rind) {
                recvenvright(buf,env->lb.rind*bs-(env->c+env->w),width); // recv env from right side
            } else {
                sendenvright(buf,(env->c+env->w)-env->lb.rind*bs,width); // send env to right
            }                
        }
    }


    if(env->lb.lind*bs!=env->c) {
        agentsshiftloc(0,env->c-env->lb.lind*bs); // shift agents by the new left index (or new env->c) current-new should shift the agents correctly
    }

    groupfree(grpbuf);
    free(buf);
}

void printblock(int i,int j,blockstruct *b) {
    printf(" (%i,%i) %f | N %f S %f E %f W %f\n",i,j,b->c,b->n,b->s,b->e,b->w);
}

int printsuperblockbuf() {
    int i,j;
    if(mpi->r!=0)
        return 0;

    printf("printsuperblockbuf\n");
    for(i=0;i<mpi->h;i++)
        for(j=0;j<env->w/env->lb.blocksize;j++) {
            printblock(i,j,&env->lb.superblockbuf[i*env->lb.wb+j]);
        }

}

int printsuperblocks() {
    int i,j;
    blockstruct *sb;

    printf("printsuperblocks\n");
    if(env->lb.dir==R)
        for(i=0;i<env->h/env->lb.blocksize;i++)
            printblock(i,0,&env->lb.superblocks[R][i]);
    else
        for(j=0;j<env->w/env->lb.blocksize;j++) 
            printblock(0,j,&env->lb.superblocks[C][j]);

}

int printblocks() {
    int i,j;
    blockstruct *b;

    if(mpi->rank!=0)
        return 0;

    printf("printblocks\n");
    for(i=0;i<env->h/env->lb.blocksize;i++)
        for(j=0;j<env->w/env->lb.blocksize;j++) {
//            b=&env->lb.blocks[i][j];
            printblock(i,j,&env->lb.blocks[i][j]);
//            printf(" (%i,%i) %f | N %f S %f E %f W %f\n",i,j,b->c,b->n,b->s,b->e,b->w);
        }
}

int edgeprocessorinexchange() {
    // if sending R, then proc 0,w-1 shouldn't wrap the environment so drop them out
    if(env->lb.pairoffset) { // odd send R (E)
        if(mpi->c==0||mpi->c==mpi->w-1) {
            return 1;
        }
    }
    return 0;
}

int sharepartition(int *optindex) {
    MPI_Status status;
    int dir; // decides the direction to exchange superblockinfo

    if(env->lb.pairoffset==0) // IMPORTANT direction is backwards, because we are sending backwards
        dir=E;
    else
        dir=W;

    // share optimal partition
    if(mpi->c%2==0) {  // send
        MPI_Send(optindex,1,MPI_INT,
                 mpi->dir[dir].send,1,mpi->dir[dir].comm->comm);
    } else { // recv (odd always receives)
        MPI_Recv(optindex,1,MPI_INT,
                 mpi->dir[dir].recv,1,mpi->dir[dir].comm->comm,&status);
    }

    // share sba portion between old partition and new partition

    //printf("shared optindex=%i\n",*optindex);
}

float calculaterowworkload(blockstruct *rb,int s,int e) {
    int i;
    float workload;

    workload=0;

    // add up E,W
    workload+=rb[s].w;
    workload+=rb[e].e;

    // sum n,s for all s-e
    for(i=s;i<=e;i++) {
        workload+=rb[i].c; // comp
        workload+=rb[i].n; // n
        workload+=rb[i].s; // s
    }

    return workload;
}

float calculatemaxworkload(blockstruct **sba,int index,int sbasize) {
    float maxl,maxr,max;
    int row;

    max=-1.0;
    for(row=0;row<mpi->h;row++) { // look at all processors in the column 
        maxl=calculaterowworkload(sba[row],0,index); // calculate l
        maxr=calculaterowworkload(sba[row],index+1,sbasize-1); // calculate r

        // Record maximum
        max=max(max,maxl);
        max=max(max,maxr);
    }

    return max;
}

// Calculate the partitions for superblock array
int calculatepartition() {
    float min,curmax;
    int index,optindex;
    int sbasize;
    blockstruct **sba;
    int indstart,indstop;

    sba=env->lb.superblockarray;

    sbasize=env->lb.sbasize;

    // Calculate the optimal partition index between the pair of processors (information in superblock array)

    optindex=-1; // initialize to invalid index
    min=FLT_MAX;

    indstart=1; // by default start with index between 0 | 1
    indstop=sbasize-1; // by default end with index between n-1 | n

    if(sbasize*env->lb.blocksize>=env->wmax-env->buf*2) { // potential conflict with partition being too big!
        // so take the overflow (sbasize-env->wmax) and add to the beginning and subtract from the end, so that no conflict happens
        indstart=1+(sbasize-((env->wmax-env->buf*2)/env->lb.blocksize));
        indstop=sbasize-1-(sbasize-((env->wmax-env->buf*2)/env->lb.blocksize))-1;
    }

    for(index=indstart;index<indstop;index++) {  //consider all partition indices
        curmax=calculatemaxworkload(sba,index,sbasize);
        if(curmax<min) { // found a new optimal minimum
            min=curmax;
            optindex=index;
        }
    }

    return optindex+1; // Note: must add 1 to optimal index, because C assumes that you are 1 above the actual cut so that you can perform < operations
}

// Copy buffer to superblock array
int buftosuperblockarray(blockstruct *buf,int width,blockstruct **sba,int offset) {
    int i,j;
    int ind;

    for(i=0;i<mpi->h;i++) {
        for(j=0;j<width;j++) {
            ind=i*env->lb.wb+j;
            sba[i][offset+j]=buf[ind];
        }
    }
}


int pairsuperblockinfo() {
    int width; // width of superblock array from remote processor
    int dir; // decides the direction to exchange superblockinfo
    int offset;

    if(env->lb.pairoffset==0)
        dir=W;
    else
        dir=E;

    // Pick the pair of processors to exchange superblock info with
    // Send/recv information
    // Make the superblockarray

    if(mpi->c%2==0) { // recv
        MPI_Status status;

        // Get width of pair processor
        MPI_Recv(&width,1,MPI_INT,
                 mpi->dir[dir].recv,1,mpi->dir[dir].comm->comm,&status);

        if(env->lb.pairoffset==0)
            buftosuperblockarray(env->lb.superblockbuf,env->w/env->lb.blocksize,env->lb.superblockarray,0);
        else
            buftosuperblockarray(env->lb.superblockbuf,env->w/env->lb.blocksize,env->lb.superblockarray,width);
            offset=width;

        // Receive new superblocks (overwrite own superblocks)
        MPI_Recv(env->lb.superblockbuf2,mpi->h*env->lb.wb*sizeof(blockstruct),MPI_BYTE,
                 mpi->dir[dir].recv,2,mpi->dir[dir].comm->comm,&status);

        if(env->lb.pairoffset==0)
            buftosuperblockarray(env->lb.superblockbuf2,width,env->lb.superblockarray,env->w/env->lb.blocksize);
        else
            buftosuperblockarray(env->lb.superblockbuf2,width,env->lb.superblockarray,0);

        env->lb.sbasize=env->w/env->lb.blocksize+width; // combine self + size

    } else { // send
        width=env->w/env->lb.blocksize;

        // Sending blocks to sender
        MPI_Send(&width,1,MPI_INT,
                 mpi->dir[dir].send,1,mpi->dir[dir].comm->comm);
        MPI_Send(env->lb.superblockbuf,mpi->h*env->lb.wb*sizeof(blockstruct),MPI_BYTE,
                 mpi->dir[dir].send,2,mpi->dir[dir].comm->comm);
        // Then don't do anything, because recv processor is doing everything ...
    }

}

// Send superblocks to selected processors (rank==0 in their column or row)
int sendsuperblocks() {
    MPI_Gather(env->lb.superblocks[env->lb.dir],env->lb.wb*sizeof(blockstruct),MPI_BYTE,
               env->lb.superblockbuf,           env->lb.wb*sizeof(blockstruct),MPI_BYTE,
               0,mpi->row.comm);
}

int sendpartitions() {
    // Send partitions to remaining processors from  (rank==0 in their column or row) to rest

    // Set left and right index to all processors in column or row
    MPI_Bcast(&env->lb.lind,1,MPI_INT,
              0,mpi->row.comm);
    MPI_Bcast(&env->lb.rind,1,MPI_INT,
              0,mpi->row.comm);
    // Sanity check
    assert(env->lb.lind<env->lb.rind);

}

// Calculate superblocks for this example
int calculatesuperblocks() {
    int i,j;
    blockstruct *sb,*b;

    env->lb.dir=C;
    memset(&env->lb.superblocks[C][0],0,sizeof(blockstruct)*env->w/env->lb.blocksize);
    if(env->lb.dir==C) {
        for(j=0;j<env->w/env->lb.blocksize;j++) { // foreach column of blocks
            sb=&env->lb.superblocks[C][j]; // get superblocks for the column
            for(i=0;i<env->h/env->lb.blocksize;i++) {
                b=&env->lb.blocks[i][j];
                sb->c+=b->c; // add all computation for each block in the superblocks
                sb->e+=b->e; // same for E and W communication
                sb->w+=b->w;
            }
            sb->n=env->lb.blocks[0][j].n; // northern comm value of northern most proc
            sb->s=env->lb.blocks[env->h/env->lb.blocksize-1][j].s; // southern comm value from southern most proc
        }
    } 
}


// Calculate block computation and communication costs
// This example only shows a very basic method for calculating weights
// Realistic ABM should have more complex calculations that consider spatial arrangements and interactions amongst entities
int calculateblocks() {
    int i,j;
    int type;
    blockstruct *b;
    float comm;
    groupstruct *g;
    int bs;
    genericstate *gs;

    bs=env->lb.blocksize;

    // ENVIRONMENT
    // Set computation and communication values for environment cells
    // This example take advantage of the fact that the blocks are the same size, but other models may not
    for(i=0;i<env->lb.hb;i++)
        for(j=0;j<env->lb.wb;j++) {
            b=&env->lb.blocks[i][j];
            b->c=ENVCOMPWEIGHT*bs*bs;
            b->n=b->s=b->e=b->w=ENVCOMMWEIGHT*bs;
        }

    // AGENTS
    // Set computation and communication values for agents
    // In this example we assign communication weights for all agents, but more realistic methods may consider distance-to-edge
    // or agent-agent / agent-environment interactions 
    for(type=0;type<LASTAGENT;type++) {
        g=POP(type).pop;
        while(g!=NULL) {
            for(i=0;i<AGENTGROUPSIZE;i++) {
                gs=agents->pop[type].ops.gs(&g->a[i]);
                if(gs->alive) {
                    b=&env->lb.blocks[gs->r/bs][gs->c/bs];
                    b->c+=AGENTCOMPWEIGHT;
                    if((gs->r%bs)>bs/2) // if agent is more south
                        b->s+=AGENTCOMMWEIGHT; // add to south comm weight
                    else
                        b->n+=AGENTCOMMWEIGHT; // else add to north
                    if((gs->c%bs)>bs/2)
                        b->e+=AGENTCOMMWEIGHT; // if more east, add to east comm weight
                    else
                        b->w+=AGENTCOMMWEIGHT; // else add to west
                }
            }
            g=g->next;
        }
    }
}

// Partition the superblocks
int partitionsuperblocks() {
    int arr[100];
    int i;
    int tste;
    int optindex;

    // Test if edge 
    if(mpi->r==0)
        tste=1;
    else
        tste=0;


    if(tste) { // If an edge

        MPI_Gather(&env->w,1,MPI_INT,
                   arr,1,MPI_INT,
                    0,mpi->col.comm);
        if(mpi->rank==0) {
            printf("Old partitions : ");
            for(i=0;i<mpi->w;i++)
                printf("%i ",arr[i]);
            printf("\n");
        }

        // odd always sends, even always receives
        // when pairoffset==0 odd sends left (west) when pairoffset==1 odd send right (east)
        // pairoffset changes the direction of the odd sender (E/W)

        env->lb.lind=env->c/env->lb.blocksize; 
        env->lb.rind=env->lb.lind+env->w/env->lb.blocksize; 

        if(!edgeprocessorinexchange()) {
            // Calculate and share partitions
            optindex=partitioncalculationandsharing();

            if(mpi->c%2) { // odd
                if(env->lb.pairoffset==0)   // 0 | 1 
                    env->lb.lind=optindex;
                else                        // 1 | 0
                    env->lb.rind=optindex;
            } else { // even
                if(env->lb.pairoffset==0)   // 0 | 1
                    env->lb.rind=optindex;
                else                        // 1 | 0
                    env->lb.lind=optindex;
            }
        }
    }
}

int partitioncalculationandsharing() {
    int optindex;

    // Select column pair and exchange superblock information, combine into superblockarray
    pairsuperblockinfo();

    // Calculate partition from superblock array
    if(mpi->c%2==0) { // Even processsor sends
        optindex=calculatepartition();

        if(env->lb.pairoffset==0) {   // 0 | 1
            optindex=env->c/env->lb.blocksize+optindex;
        } else  {                       // 1 | 0
            optindex=(env->c+env->w)/env->lb.blocksize-env->lb.sbasize+optindex; // get abs (c+w)/bs - sbasize to get far left, then add optindex to get index
        }
    }

    // Share information with other processor in pair
    // Share sba portion too, based on optindex
    sharepartition(&optindex);

    return optindex;
}

// Decide whether or not to load balance using a computationally efficient method
// This example only considers load balancing if a single core has a set number of more agents than average
// More complex decisions should be made for other models, this is only an illustration

int decide() {
    int cnt,gcnt;
    int min,max;
    int lb;
    int threshold;

    // Calculate local,maximum, and global counts
    cnt=agentscount();
    MPI_Reduce(&cnt,&max,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&cnt,&gcnt,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    // Set threshold to 1.25 more than average to trigger load balancing
    // Other methods consider min/min, spatial distribution of entities, etc.
    // This is a very basic example to illustrate the process of making a decision 
    threshold=(gcnt/mpi->size)*1.25;

    if(max>threshold)
        lb=1;
    else
        lb=0;

    // Share decision
    MPI_Bcast(&lb,1,MPI_INT,0,MPI_COMM_WORLD);

    return lb;
}

// Load balancer that was simplified for example code
int loadbalance() {
    int optindex,i,arr[100];

    // Phase 1: Decide to load balance or not 
    if(decide()) {
        ADMIN(0,"Load balancing ...");
        ADMIN(0,"Invoking communication-aware load-balancing strategy");
    } else {
        ADMIN(0,"Decided to skip load balancing, because the load imbalance is small");
        return 0;
    }

    // Phase 2: calculate blocks for localenvironment
    calculateblocks();

    // Phase 3: Create superblocks
    calculatesuperblocks();

    // Phase 4: send superblocks to selected processors
    sendsuperblocks();

    // Phase 5: partition superblocks
    partitionsuperblocks();

    // Phase 6: Share partitions with other processors
    sendpartitions();

    // Phase 7: Shift environment and agent entities around to balance the workload
    shiftenvandagents();

    // Reset env variables based on new indices
    env->c=env->lb.lind*env->lb.blocksize;
    env->w=(env->lb.rind-env->lb.lind)*env->lb.blocksize;
    env->wabs=env->w+env->buf*2;

    // after shifting the environment for the agents - migrate the corresponding population
    popmigrate(&POP(0));

    // Sanity checks for migration
    #ifndef NDEBUG
    gaiter(POP(0).pop,&agentidentifymigrating);
    agentiter(POP(0).migratingagents,POP(0).behavior.info);
    assert(POP(0).migratingagents==NULL);
    #endif
    
    // Reset ghost zones based on new subdomains
    envdelghostzones();
    envinitghostzones();

    // Repopulate ghost zones
    envupdateghostzones();

    // Alternate directions sending/receiving (even/odd)
    env->lb.pairoffset=(env->lb.pairoffset+1)%2; 

    MPI_Gather(&env->w,1,MPI_INT,
               arr,1,MPI_INT,
                0,mpi->col.comm);
    if(mpi->rank==0) {
        printf("New partitions : ");
        for(i=0;i<mpi->w;i++)
            printf("%i ",arr[i]);
        printf("\n");
    }
    envprint();
    printf("  Core[%i] Number of agents : %i; Number of environment cells : %i\n",mpi->rank,agentscount(),env->h*env->w);
}


// Initialize load-balancing structures
int lbinit() {
    int i;
    int h,w;

    ADMIN(2,"lb init");

    // Assign blocksize
    env->lb.blocksize=BLOCK;

    // Calculate height/width of block arrays
    env->lb.hb=env->hmax/env->lb.blocksize;
    env->lb.wb=env->wmax/env->lb.blocksize;

    // Alloc the structures
    env->lb.blocks=(blockstruct **)calloc(env->lb.hb,sizeof(blockstruct *));
    for(i=0;i<env->lb.hb;i++)
        env->lb.blocks[i]=(blockstruct *)calloc(env->lb.wb,sizeof(blockstruct));

    // Alloc structures for supeblocks
    env->lb.superblocks[R]=(blockstruct *)calloc(env->lb.hb,sizeof(blockstruct));
    env->lb.superblocks[C]=(blockstruct *)calloc(env->lb.wb,sizeof(blockstruct));

    // Setup superblock array on edge processors
    if(mpi->r==0||mpi->c==0) { 
        env->lb.superblockbuf=(blockstruct *)calloc(mpi->h*env->lb.wb,sizeof(blockstruct));
        env->lb.superblockbuf2=(blockstruct *)calloc(mpi->h*env->lb.wb,sizeof(blockstruct));
        env->lb.superblockarray=(blockstruct **)calloc(mpi->h,sizeof(blockstruct *));
        for(i=0;i<mpi->h;i++)
            env->lb.superblockarray[i]=(blockstruct *)calloc(env->lb.wb*2,sizeof(blockstruct));  
    }

    env->lb.pairoffset=0; // recv=even,send=odd
}

// Clear load balancing data structures
int lbdel() {
    int i;
    for(i=0;i<env->lb.hb;i++)
        free(env->lb.blocks[i]);
    free(env->lb.blocks);

    free(env->lb.superblocks[R]);
    free(env->lb.superblocks[C]);

    if(mpi->r==0||mpi->c==0) {
        free(env->lb.superblockbuf);
        free(env->lb.superblockbuf2);
        for(i=0;i<mpi->h;i++)
            free(env->lb.superblockarray[i]);
        free(env->lb.superblockarray);
    }

}

