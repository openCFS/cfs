/* headerfile for capapost */
#ifndef UNV_DAT_H
#define UNV_DAT_H


#include "unv.hh"
#include "General/defs.hh"

#define MAX_NUM_NODES 1000
#define MAX_NUM_ELEMENTS 1000


#ifdef NO_CAPAPOST

#define Boolean int
#define Widget int
#define XtMalloc malloc
#define XtCalloc calloc
#define XtRealloc realloc
#define XtFree free
#define workData int    // dummy-Definition fuer wprintf
int XtParent(workData);
#define FALSE 0         // wir machen es mal so, wie in C ueblich
#define TRUE 1          // um mit Akuflow kompatibel zu sein

#endif

/*   structures           */

typedef struct {
  Double x1,x2,x3;
} point;

struct element_data {
  int label;
  int color;
  int n_nodes;
  int fe_type;
  int *points;
  int *edges;
};

struct node_dat {
  Double data[6];
};

struct set_55 {
  header_55 header;
  node_dat *dat;
};

struct elem_dat {
  Double data[12];
};

struct set_56 {
  header_56 header;
  elem_dat *dat;
};

struct set_58 {
  int capaDataType;
  header_58 header;
  data_58 *data;
};

struct history_node {
  int num_datasets;
  set_58 *sets;
};


extern char* UNV_FILE;

//---- Node data ----
extern point* NODES;
extern point *DEF_NODES;
extern int N_NODES;
extern Double NODES_X_MAX, NODES_X_MIN, NODES_Y_MAX, NODES_Y_MIN;

//---- Element data -----
extern element_data* ELEMENTS;
extern int N_ELEMENTS;

// introduced by Simon Triebenbacher
//---- Dimension of the mesh ----
extern int DIM;

// --- Data structures -----
extern set_55* SETS55;
extern set_56* SETS56;
extern history_node *NODES58;
extern set_56* ELEMENT_CURRENT;
extern int N_SETS55;
extern int N_SETS56;
extern set_55* CURRENT;

extern int* DISP;
extern int N_DISP;

extern int N_EDGES;
extern int* EDGE_DATA;
extern int* EDGE_COLOR;

extern int* NODE_LABELS;
extern int  MAX_NODE_LABEL;

extern workData work;

extern Boolean CheckForInterrupt(void);
extern void TimeoutCursors(Boolean, Boolean);
extern void wprint(workData work, char *buffer);
extern void ErrorDialog(Widget w, const char* error);
extern int ReadUniversalFile(Widget w);

#endif

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
