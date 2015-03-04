// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "hdr.h"

int cellinit(cellstruct *cell,int i,int j);
int envinit();
int envdel();
int envupdatelocal();
int envupdateghostzone();

// REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE REMOTE 

// Calculate global height and width
int envgetglobal() {
    MPI_Allreduce(&env->h,&env->hg,1,MPI_INT,MPI_SUM,mpi->row.comm);
    MPI_Allreduce(&env->w,&env->wg,1,MPI_INT,MPI_SUM,mpi->col.comm);
}

// Update the environment values in ghost zones
int envupdateghostzones() {
    int i;
    tilestruct *tile; 
    MPI_Status status;
    ADMIN(2,"Env update ghost zone");

    // Update the buffers for each direction
    for(i=0;i<LASTDIR;i++) {
        tile=&mpi->dir[i].tilesend;
        // Group environment cells for a remote ghost zone region
        blocktobuffer(env->grid,mpi->dir[i].envbuf,tile);

        // Send group using NS,EW strategy
        MPI_Sendrecv_replace(mpi->dir[i].envbuf,tile->h*tile->w,mpi->MPI_CELL,
                            mpi->dir[i].send,200+i,
                            mpi->dir[i].recv,200+i,
                            mpi->dir[i].comm->comm,&status);
        tile=&mpi->dir[i].tilerecv;

        // Copy group of environments cells to appropraite ghost zone region
        buffertoblock(mpi->dir[i].envbuf,env->gridabs,tile);
    }
}

// Initialize ghost zones
int envinitghostzones() {
    int i;
    tilestruct *ts,*tr;
    int buf;

    DEBUG(3,printf("envinitghostzones()\n"));

    buf=env->buf;

    // This section sets the environment buffers in the mpi struct
    for(i=0;i<LASTDIR;i++) {
        ts=&mpi->dir[i].tilesend; // Tile send is in relation send grid
        ts->r=0;ts->c=0; // default r,c to 0
        ts->h=buf;ts->w=buf; // default h,w to buf
        switch(i) {
            case N:
                    ts->w=env->w;
                    break;
            case S:
                    ts->r=env->h-buf;
                    ts->w=env->w;
                    break;
            case E:
                    ts->c=env->w-buf;
                    ts->h=env->h;
                    break;
            case W:
                    ts->h=env->h;
                    break;
        }

        tr=&mpi->dir[i].tilerecv; // Tile recv is in relation send grid absolute
        tr->r=buf;tr->c=buf; // default r,c to buf (remember in relation to abs grid)
        tr->h=buf;tr->w=buf; // default h,w to buf

        switch(i) {
            case N: // recv S
                    tr->r=env->habs-buf;
                    tr->w=env->w;
                    break;
            case S: // recv N
                    tr->r=0;
                    tr->w=env->w;
                    break;
            case E: // recv W
                    tr->c=0;
                    tr->h=env->h;
                    break;
            case W: // recv E
                    tr->c=env->wabs-buf;
                    tr->h=env->h;
                    break;
        }

        // Allocate buffers for the structures
        mpi->dir[i].envbuf=(cellstruct *)malloc(ts->h*ts->w*sizeof(cellstruct));
        assert(mpi->dir[i].envbuf!=NULL);
    }

}

// Clear the ghost zones
int envdelghostzones() {
    int i;
    for(i=0;i<LASTDIR;i++) {
        free(mpi->dir[i].envbuf);
    }
}


// Copy a 2D block of cells to a 1D buffer
int blocktobuffer(cellstruct **block,cellstruct *buf,tilestruct *tile) {
    //Copy a rectangle size (hXw) at (rXc) to a 1D buffer with the size h*w
    //r=h=i and c=w=j
    int i,j;
    for(i=0;i<tile->h;i++)
        for(j=0;j<tile->w;j++) {
            buf[tile->w*i+j]=block[i+tile->r][j+tile->c];
        }
}

// Copy a 1D buffer to a 2D block of cells
int buffertoblock(cellstruct *buf,cellstruct **block,tilestruct *tile) {
    //Copy a 1D buffer with size h*w to a rectangle size (hXw) at (rXc) in grid 
    //r=h=i and c=w=j
    //printf("buffertoblock (%i,%i) by [%i,%i]\n",r,c,h,w);
    int i,j;

    for(i=0;i<tile->h;i++)
        for(j=0;j<tile->w;j++) {
            block[i+tile->r][j+tile->c]=buf[tile->w*i+j];
        }
}

// Copy a 2D block to another 2D block
int blocktoblock(cellstruct **blockrecv,cellstruct **blocksend,tilestruct *tilerecv,tilestruct *tilesend) {
    // Copy a rectangle size (hXw) from grid at (frXfc) to grid at (trXtc)
    //r=h=i and c=w=j
    //printf("blocktoblock copying from(%i,%i) to (%i,%i) by [%i,%i]\n",fr,fc,tr,tc,h,w);
    int i,j;
    assert(tilesend->h==tilerecv->h);
    assert(tilesend->w==tilerecv->w);
    for(i=0;i<tilesend->h;i++)
        for(j=0;j<tilesend->w;j++)
            blocksend[i+tilesend->r][j+tilesend->c]=blockrecv[i+tilerecv->r][j+tilerecv->c];
}


// LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL LOCAL 

// Get a new row for environment
cellstruct *envgetrow() {
    cellstruct *c;
    c=(cellstruct *)calloc(env->wmax,sizeof(cellstruct));
    return c;
}

int envprint() {
    printf("  Core[%i] Environment configuration (%i,%i) [%i,%i]\n",mpi->rank,env->r,env->c,env->h,env->w);
}

int tileprint(tilestruct *t) {
    printf(" (%i,%i) [%i,%i]\n",t->r,t->c,t->h,t->w);
}

int envprintgrid() {
    int i,j;
    printf("envprintgrid()\n");
    for(i=0;i<env->h;i++) {
        for(j=0;j<env->w;j++) {
            printf(" %3.0f",env->grid[i][j].val);
        }
        printf("\n");
    }

}
int envprintgridabs() {
    int i,j;
    printf("envprintgridabs()\n");
    for(i=0;i<env->habs;i++) {
        for(j=0;j<env->wabs;j++) {
            printf(" %4.0f ",env->gridabs[i][j].val);
        }
        printf("\n");
    }

}

// Update a cell
inline int cellupdate(cellstruct *cell,int r,int c) {
    // Add a fixed value, capped by a value
    cell->val=fmin(cell->val+0.2,cell->cap);
}

// Set the cell with set values
inline int cellset(cellstruct *cell,float v) {
    cell->val=v;
    cell->cap=v;
}

inline int cellinit(cellstruct *cell,int r,int c) {
    // Randomly assign a value and capacity from 0-100 
    float cap=drand48()*(float)100.0;
    cellset(cell,cap);

}

int celldel() {} // Does not have any function currently 

int envupdatelocal() {
    int i,j;
    float v,cap;

    ADMIN(2,"Env update local");

    for(i=0;i<env->h;i++)
        for(j=0;j<env->w;j++) 
            cellupdate(&env->grid[i][j],i,j);
}

// Setup pointers for gridabs and grid based on location in gridmax
int envupdategrids() {
    int i;

    for(i=0;i<env->habs;i++)
        env->gridabs[i]=&env->gridmax[i][0];
    for(i=0;i<env->h;i++)
        env->grid[i]=&env->gridabs[env->buf+i][env->buf];

}

int envprintenv() {
    printf(" Local env (%i,%i) [%i,%i]\n",env->r,env->c,env->h,env->w);
}

// ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV // ENV 

int envupdate() {
    ADMIN(1,"Env update");

    // First update the local environment
    envupdatelocal();

    // Then update the ghost zones
    envupdateghostzones();
}

int envinit() {
    int i,j;

    DEBUG(2,printf("Env init\n"););

    // Each subdomain is sized NUMROWSxNUMCOLS
    env->h=NUMROWS;
    env->w=NUMCOLS;

    // Determine first row and column index
    env->r=env->h*mpi->row.rank;
    env->c=env->w*mpi->col.rank;

    // Calculate ghost zone buffer size (set as fixed variable in this example)
    env->buf=INFLUENCEDIST;

    // Current absolute height and width (consider height/width of subdomain plus buffers
    env->habs=env->h+env->buf*2;
    env->wabs=env->w+env->buf*2;

    // Maximum possible height/width allowed by load balancing and decomposition 
    // Based on current subdomain size in this example
    env->hmax=env->h*RMULT;
    env->wmax=env->w*CMULT;

    // Ensure maximum possible is at least as big as current size
    if(env->hmax<env->habs)
        env->hmax=env->habs;
    if(env->wmax<env->wabs)
        env->wmax=env->wabs;

    // Setup max grid
    env->gridmax=(cellstruct **) malloc((env->hmax)*sizeof(cellstruct *));
    assert(env->gridmax!=NULL);
    for(i=0;i<env->hmax;i++)
        env->gridmax[i]=envgetrow();

    // Setup absolute grid
    env->gridabs=(cellstruct **) malloc((env->habs)*sizeof(cellstruct *));
    assert(env->gridabs!=NULL);

    // Setup grid inside absolute grid (subtract buffers)
    env->grid=(cellstruct **) malloc((env->h)*sizeof(cellstruct *));
    assert(env->grid!=NULL);

    // Set environment variables (local and global)
    envupdategrids(); // set gridabs and grid
    envgetglobal();   // calculate global dimensions of environment

    // Initialize values in grid (ghost zones will be updated soon) 
    for(i=0;i<env->h;i++)
        for(j=0;j<env->w;j++) 
            cellinit(&env->grid[i][j],i,j);

    // Setup ghost zone regions
    envinitghostzones();
    envupdateghostzones();

/*
envprintgrid();
printf("\n##############\n\n");
printf("##############\n\n");
envprintgridabs();
*/
    //mpi->MPI_CELL is created in mpi.c (mpiinitcelltype)

    //EXEC0(printf(" Environment (local) [%i,%i] \n",env->h,env->w));
    EXEC0(printf(" Environment is rectilinearly decomposed (spatial extent) [%i, %i] \n",env->hg,env->wg));
    envprint();

    // Initialize the load balancer
    if(LOADBALANCE)
        lbinit();

}

// Cleanup environment
int envdel() {
    int i;

    ADMIN(1,"Env del");

    if(LOADBALANCE)
    lbdel();

    //mpi->MPI_CELL is freed in mpi.c

    envdelghostzones();

    free(env->grid);
    free(env->gridabs);
    for(i=0;i<env->hmax;i++)
        free(env->gridmax[i]);
    free(env->gridmax);
}




