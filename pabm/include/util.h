// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.

// Utility macros

#define EXEC0(X) if(mpi->rank==0) X

#define ADMIN(dbg,msg) if(DEBUGLVL>=dbg) EXEC0(printf(" [ %s ]\n",msg);)
