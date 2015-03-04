// © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
#include "agent.h"

// Behaviors for sugar agents
typedef struct {
    int (*init)();
    int (*move)();
    int (*update)();
    int (*del)();
} behaviorxplstr;

typedef struct {
    int astate;
} xplstate;

typedef struct {
    genericstate g;
    xplstate s;
} xplstruct;



