#include <cstdlib>
#include <cstdio>

#include "SimInputUnv.hh"
#include "unv_dat.hh"

char *UNV_FILE=NULL;

point *NODES=NULL;
int N_NODES=0;
double NODES_X_MAX, NODES_X_MIN, NODES_Y_MAX, NODES_Y_MIN;

element_data* ELEMENTS=NULL;
int N_ELEMENTS=0;

// introduced by Simon Triebenbacher
int DIM=0;

set_55* SETS55=NULL;
set_56* SETS56=NULL;
int N_SETS55=0;
int N_SETS56=0;
set_55* CURRENT=NULL;
set_56* ELEMENT_CURRENT=NULL;

int* DISP=NULL;
int N_DISP=0;

int N_EDGES=0;
int* EDGE_DATA=NULL;
int* EDGE_COLOR=NULL;

int* NODE_LABELS = NULL;
int  MAX_NODE_LABEL = 0;

workData work;


#ifdef NO_CAPAPOST
//
// Dummy-Routinen zum Linken ohne X11
//
void wprint(workData w, char* s)
{
  std::cerr << s << std::endl;
}
int XtParent(int i) {return i; }
void ErrorDialog(workData w, const char* s) { EXCEPTION(s); }
void TimeoutCursors(int i, int j) {; }
int CheckForInterrupt(void) { return 0; }

#endif

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
