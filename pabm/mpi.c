// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include <math.h>
#include "hdr.h"

int mpidelcelltype() {
    MPI_Type_free(&mpi->MPI_CELL);
}

int mpiinitcelltype() {
    int i;
    // This section sets the MPI_CELL in mpi so that we do not have to pack and unpack the environment datatype
    cellstruct cell[1];

    // This code is somewhat complex to read, but defines an equivalent to a struct in MPI
    // This will need to be changed as environment cell structure changes
    // Creating an MPI_Type oftentimes results in slightly better performance, but most importantly
    // enables heterogeneous machines to send/receive data when using MPI.
    // An alternative to this approach is to define everything as MPI_BYTE and use sizeof() operator to send/recv
    MPI_Datatype cell_types[CELLSTRUCTELEMENTS]={MPI_FLOAT,MPI_FLOAT,MPI_UB};
    int cell_bl[CELLSTRUCTELEMENTS]=            {        1,        1,     1};
    MPI_Aint cell_ind[CELLSTRUCTELEMENTS],base;
    MPI_Address(&cell[0],      &base);
    MPI_Address(&cell[0].val,  &cell_ind[0]);
    MPI_Address(&cell[0].cap,  &cell_ind[1]);
    cell_ind[CELLSTRUCTELEMENTS-1]=sizeof(cellstruct);

    // Must subtract base address
    for(i=0;i<CELLSTRUCTELEMENTS-1;i++)
        cell_ind[i]-=base;

    // Define CELL as an MPI struct
    MPI_Type_struct(CELLSTRUCTELEMENTS,cell_bl,cell_ind,cell_types,&mpi->MPI_CELL);
    MPI_Type_commit(&mpi->MPI_CELL);
}

int mpiupdate() { } // Currently nothing here

int mpiinit() {
    int i,j;
    int wrap[2],dim[2],coords[2];
    int reorder;

    // Initialize global MPI environment
    MPI_Init(&sim->argc,&sim->argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi->rank);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi->size);

    // Find a working grid, even if 1xC
    mpi->h=floor(sqrt(mpi->size));
    mpi->w=mpi->size/mpi->h;
    while(mpi->size>mpi->w*mpi->h) { 
        mpi->h--;
        mpi->w=mpi->size/mpi->h;
    }

    // Sanity check, for this simple case
    if(mpi->size>mpi->h*mpi->w) {
        EXEC0(printf(" [ ERROR ] Processors cannot be made into a semi-square grid of (%ix%i) with %i processors\n",mpi->h,mpi->w,mpi->size));
        MPI_Finalize();
        exit(0);
    }

    // Prepare environment for MPI 
    wrap[0]=wrap[1]=1;
    dim[0]=mpi->h; // Set height
    dim[1]=mpi->w; // Set width
    reorder=0; // Set reorder to 0, because decomposition is explicitly managed 

    //Setup cartesian grid
    MPI_Cart_create(MPI_COMM_WORLD,2,dim,wrap,reorder,&(mpi->cart.comm));
    MPI_Comm_rank(mpi->cart.comm,&mpi->cart.rank);
    MPI_Comm_size(mpi->cart.comm,&mpi->cart.size);
    MPI_Cart_coords(mpi->cart.comm,mpi->rank,2,coords);
    mpi->r=coords[0];
    mpi->c=coords[1];

    //Setup row
    dim[0]=1;
    dim[1]=0;
    MPI_Cart_sub(mpi->cart.comm,dim,&mpi->row.comm);
    MPI_Comm_rank(mpi->row.comm,&mpi->row.rank);
    MPI_Comm_size(mpi->row.comm,&mpi->row.size);
    mpi->r=mpi->row.rank;
    mpi->h=mpi->row.size;

    //Setup col
    dim[0]=0;
    dim[1]=1;
    MPI_Cart_sub(mpi->cart.comm,dim,&mpi->col.comm);
    MPI_Comm_rank(mpi->col.comm,&mpi->col.rank);
    MPI_Comm_size(mpi->col.comm,&mpi->col.size);
    mpi->c=mpi->col.rank;
    mpi->w=mpi->col.size;

    // Initialize all the communicators 
    // Note: Environment bufs and ghost zones are setup in env.c
    for(i=0;i<LASTDIR;i++) {
        if(i==N||i==S) {
            mpi->dir[i].comm=&mpi->row;
        } else if(i==E||i==W) {
            mpi->dir[i].comm=&mpi->col;
        }

        // Determine which cores to send/recv to based on locations in cartesian grid
        switch(i) {
            case N:
                mpi->dir[i].send=(mpi->r-1+mpi->h)%mpi->h;
                mpi->dir[i].recv=(mpi->r+1)%mpi->h;
                break;
            case S:
                mpi->dir[i].send=(mpi->r+1)%mpi->h;
                mpi->dir[i].recv=(mpi->r-1+mpi->h)%mpi->h;
                break;
            case E: 
                mpi->dir[i].send=(mpi->c+1)%mpi->w;
                mpi->dir[i].recv=(mpi->c-1+mpi->w)%mpi->w;
                break;
            case W:
                mpi->dir[i].send=(mpi->c-1+mpi->w)%mpi->w;
                mpi->dir[i].recv=(mpi->c+1)%mpi->w;
                break;
            default:
                printf("[ ERROR ] Case statement problem: %i\n",i);
        }
    }

    // Set seed to be system time, in seconds
    sim->seed=time(NULL);

    // Print simulation information now that rank==0 is established
    EXEC0(printf("Simulation iterations: %i\n",sim->maxiter)); 
    EXEC0(printf("Random seed: %i (Note: seed is multiplied by mpi->rank)\n",sim->seed)); 
    EXEC0(printf("Processor arrangement: [%i,%i]\n",mpi->h,mpi->w));
    #ifndef NDEBUG 
        EXEC0(printf("ASSERT statements on\n"));
    #endif
    if(DEBUGLVL>0)
        EXEC0(printf("DEBUG set to %i\n",DEBUGLVL));


    // Broadcast seed to all cores (from core 0)
    MPI_Bcast(&sim->seed,1,MPI_INT,0,MPI_COMM_WORLD);

    // Offset seed by rank (to ensure all cores don't follow exact same random sequence)
    sim->seed+=sim->seed*mpi->rank;

    // Seed the random number generators
    srand(sim->seed);
    srand48(sim->seed);

    // Define an environment cell as an MPI data type
    mpiinitcelltype();
}

int mpidel() {
    int i,j;

    ADMIN(1,"MPI del");

    // Barrier to sync all processors before teardown (without this you get segfaults)
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm_free(&mpi->col.comm);
    MPI_Comm_free(&mpi->row.comm);
    MPI_Comm_free(&mpi->cart.comm);

    mpidelcelltype();

    MPI_Finalize();
}



