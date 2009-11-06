/* C interface for snopt */

#ifndef SNOPTINTERFACE_HH_
#define SNOPTINTERFACE_HH_

#include <stdio.h>
#include "f2c.h"

/*open*/
typedef struct
{ 
  short oerr;
  short ounit;
  char *ofnm;
  int ofnmlen;
  char *osta;
  char *oacc;
  char *ofm;
  short orl;
  char *oblnk;
} olist;

/*close*/
typedef struct
{
  short cerr;
  short cunit;
  char *csta;
} cllist;

extern "C" {

typedef int (*U_vfp)();
typedef int (*U_fp)(int*, int*, double*, int*, int*, double*, int*, int*, double*, char*, int*, int*, int*, double*, int*);

int snopta_(
    int *start, int *nef, int *n, 
    int *nxname, int *nfname, double *objadd, const int *objrow, char *prob,
    U_fp usrfun, U_vfp snlog, int *iafun, int *javar, 
    int *lena, int *nea, double *a, int *igfun, int *jgvar, int *leng, int *neg, 
    double *xlow, double *xupp, char *xnames, double *flow, double *fupp, char *fnames, 
    double *x, int *xstate, double *xmul, double *f, 
    int *fstate, double *fmul, int *inform__, int *mincw, 
    int *miniw, int *minrw, int *ns, int *ninf, 
    double *sinf, char *cu, int *lencu, int *iu,
    int *leniu, double *ru, int *lenru, char *cw, int *lencw, 
    int *iw, int *leniw, double *rw, int *lenrw,
    int prob_len, int xnames_len, int fnames_len, int cu_len, int cw_len
);

int sninit_(
    int *iPrint, int *iSumm,
    char *cw, int *lencw, int *iw, int *leniw,
    double *rw, int *lenrw, int cw_len
);

int sngeti_(
    const char *buffer, int *ivalue, int *inform__,
    char *cw, int *lencw, int *iw, 
    int *leniw, double *rw, int *lenrw, 
    int buffer_len, int cw_len
);

int snset_(
    const char *buffer, int *iprint, int *isumm, 
    int *inform__, char *cw, int *lencw, 
    int *iw, int *leniw, 
    double *rw, int *lenrw, 
    int buffer_len, int cw_len
);

int snseti_(
    const char *buffer, int *ivalue, int *iprint, 
    int *isumm, int *inform__, char *cw, 
    int *lencw, int *iw, int *leniw, 
    double *rw, int *lenrw, int buffer_len,
    int cw_len
);

int snsetr_(
    const char *buffer, double *rvalue, int * iprint, 
    int *isumm, int *inform__, char *cw, 
    int *lencw, int *iw, int *leniw, 
    double *rw, int *lenrw, int buffer_len, 
    int cw_len
);

int snspec_(
    int *ispecs, int *inform__, char *cw, 
    int *lencw, int *iw, int *leniw, 
    double *rw, int *lenrw, int cw_len
);

int snmema_(
    int *nef, int *n, int *nxname, 
    int *nfname, int *nea, int *neg,
    int *mincw, int *miniw, int *minrw,
    char *cw, int *lencw, int *iw, 
    int *leniw, double *rw, int *lenrw, 
    int cw_len
);

int snjac_(
    int *inform__, int *iprint, int *
    isumm, int *nef, int *n, U_fp userfg, 
    int *iafun, int *javar, int *lena, 
    int *nea, double *a, int *igfun, 
    int *jgvar, int *leng, int *neg, 
    double *x, double *xlow, double *xupp,
    int *mincw, int *miniw, 
    int *minrw, char *cu, int *lencu, 
    int *iu, int *leniu, double *ru, 
    int *lenru, char *cw, int *lencw, int *iw, 
    int *leniw, double *rw, int *lenrw, 
    int cu_len, int cw_len
);

int snlog_();

int snoptclose_(int *iunit);

int snoptopen_(int *iunit);

#undef real
} // extern "C"

#endif
