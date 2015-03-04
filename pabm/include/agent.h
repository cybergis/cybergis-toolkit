// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#ifndef DEFAGENTHEADER
#define DEFAGENTHEADER

// Agent boundary check macros 
#define NPBNDCHECK(r,c) (r<env->buf)
#define SPBNDCHECK(r,c) (r>=env->h-env->buf)
#define EPBNDCHECK(r,c) (c>=env->w-env->buf)
#define WPBNDCHECK(r,c) (c<env->buf)
#define PBNDCHECK(r,c) NPBNDCHECK(r,c)||SPBNDCHECK(r,c)||EPBNDCHECK(r,c)||WPBNDCHECK(r,c)

#define NBNDCHECK(r,c) (r<0)
#define SBNDCHECK(r,c) (r>=env->h)
#define EBNDCHECK(r,c) (c>=env->w)
#define WBNDCHECK(r,c) (c<0)
#define BNDCHECK(r,c)  NBNDCHECK(r,c)||SBNDCHECK(r,c)||EBNDCHECK(r,c)||WBNDCHECK(r,c)

// Generic structures for agents
typedef struct {
    int id;
    int r,c;
    int alive;
} genericstate;

#endif

