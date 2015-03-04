// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "hdr.h"

// Global variables 
mpistruct *mpi;
envstruct *env;
agentsstruct *agents;
simstruct *sim;
simstruct simstr;

int main(int argc, char **argv) {
    mpistruct mpistr;
    envstruct envstr;
    agentsstruct agentsstr;
    double iterstart,iterstop;

    // Initialize global variables
    mpi=&mpistr;
    env=&envstr;
    agents=&agentsstr;
    sim=&simstr;

    // Initialize sim,agents,env,mpi
    init(argc,argv);

    // Run the simulation
    iterstart=MPI_Wtime();
    for(sim->iteration=0;sim->iteration<sim->maxiter;sim->iteration++) {
        MPI_Pcontrol(1,"iteration");

        MPI_Pcontrol(1,"sim");
        simupdate();
        MPI_Pcontrol(-1,"sim");

        MPI_Pcontrol(1,"env");
        envupdate();
        MPI_Pcontrol(-1,"env");

        MPI_Pcontrol(1,"agent");
        agentsupdate();
        MPI_Pcontrol(-1,"agent");

        MPI_Pcontrol(-1,"iteration");
    }

    // Make sure all cores are finished
    MPI_Barrier(MPI_COMM_WORLD);
    iterstop=MPI_Wtime();
    EXEC0(printf("Iteration time=%f\n",iterstop-iterstart));

    // Cleanup
    del(); 

    return 0;
}

void EXIT() {
    if(sim->iteration<sim->maxiter) {
        printf("EXITing\n");
        MPI_Abort(MPI_COMM_WORLD,-1);
        del();
    } else {
        EXEC0(printf("Simulation exiting ...\n"));
        EXEC0(printf("IPM captures simulation performance characteristics, please refer to documentation for compilation instructions\n"));
    }
}

int init(int argc, char **argv) {

    siminit(argc,argv);
    mpiinit();
    envinit();
    agentsinit();

    atexit(&EXIT); // exit at the end 
}

int del() {
    agentsdel();
    envdel();
    mpidel();
    simdel();
}

int DIE(char *str) {
    printf(" [ Died : %s ]\n",str);
    del();
    return 13;
}


