/* headerfile for capapost */
#ifndef UNV_DAT_H
#define UNV_DAT_H


#include "unv.h"

#include "def_capapost.h"

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
    double x1,x2,x3;
} point;

struct element_data {
    int label;
    int color;
    int n_nodes;
    int *points;
    int *edges;
};

struct node_dat {
    double data[6];
};

struct set_55 {
    header_55 header;
    node_dat *dat;
};

struct elem_dat {
    double data[12];
};

struct set_56 {
    header_56 header;
    elem_dat *dat;
};


extern char* UNV_FILE;

//---- Node data ----
extern point* NODES;
extern point *DEF_NODES;
extern int N_NODES;
extern double NODES_X_MAX, NODES_X_MIN, NODES_Y_MAX, NODES_Y_MIN;

//---- Element data -----
extern element_data* ELEMENTS;
extern int N_ELEMENTS;

// --- Data structures -----
extern set_55* SETS55;
extern set_56* SETS56;
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
extern void ErrorDialog(Widget w, char* error);
extern int ReadUniversalFile(Widget w);

#endif












