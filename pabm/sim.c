// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "hdr.h"

// Update simulation each iteration
int simupdate() {

    // Print iteration information
    if(DEBUGLVL>0) {
        MPI_Barrier(MPI_COMM_WORLD); // This barrier is included to maintain output printing consistency 
        EXEC0(printf("\t\t\t\t\t\t                    [ Iteration: %i ]\n",sim->iteration));
    } else 
    if(!(sim->iteration%5)) {
        EXEC0(printf(" [ Iteration: %i ]\n",sim->iteration));
    }
    
    // If load balancing is enabled
    if(LOADBALANCE) {
        MPI_Pcontrol(1,"lb");
        // In this example, load balancing checks happen by a given frequency
        if((sim->iteration%(LB_FREQ))==0) {
            loadbalance();
        }
        MPI_Pcontrol(-1,"lb");
    }
}

// Initialize a simulation
int siminit(int argc, char **argv) {
    sim->argc=argc;
    sim->argv=argv;
    sim->iteration=0;
    sim->maxiter=ITERATIONS;
}

int simdel() { 
    // Nothing to do
} 


