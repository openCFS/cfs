/* C interface for snopt */

#ifndef SNOPTINTERFACE_HH_
#define SNOPTINTERFACE_HH_

#include <stdio.h>
#include "f2c.h"

typedef double doublereal;

extern "C" {

typedef int (*U_fp)(int*, int*, double*, int*, int*, double*, int*, int*, double*, char*, int*, int*, int*, double*, int*);

 int snopta_
( integer *start, integer *nef, integer *n,
  integer *nxname, integer *nfname, doublereal *objadd, const integer *objrow,
  char *prob, U_fp usrfun, integer *iafun, integer *javar,
  integer *lena, integer *nea, doublereal *a, integer *igfun,
  integer *jgvar, integer *leng, integer *neg,
  doublereal *xlow, doublereal *xupp,
  char *xnames, doublereal *flow, doublereal *fupp, char *fnames,
  doublereal *x, integer *xstate, doublereal *xmul, doublereal *f,
  integer *fstate, doublereal *fmul, integer *inform, integer *mincw,
  integer *miniw, integer *minrw, integer *ns, integer *ninf,
  doublereal *sinf, char *cu, integer *lencu, integer *iu, integer *leniu,
  doublereal *ru, integer *lenru, char *cw, integer *lencw,
  integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen
  prob_len, ftnlen xnames_len, ftnlen fnames_len, ftnlen cu_len, ftnlen cw_len );

 int sninit_
( integer *iPrint, integer *iSumm, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  double *rw, integer *lenrw, ftnlen cw_len );


 int sngeti_
( char *buffer, integer *ivalue, integer *inform,
  char *cw, integer *lencw, integer *iw,
  integer *leniw, double *rw, integer *lenrw,
  ftnlen buffer_len, ftnlen cw_len);

 int snset_
( const char *buffer, integer *iprint, integer *isumm,
  integer *inform, char *cw, integer *lencw,
  integer *iw, integer *leniw,
  double *rw, integer *lenrw,
  ftnlen buffer_len, ftnlen cw_len );

 int snseti_
( const char *buffer, integer *ivalue, integer *iprint,
  integer *isumm, integer *inform, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  double *rw, integer *lenrw, ftnlen buffer_len,
  ftnlen cw_len );

 int snsetr_
( const char *buffer, double *rvalue, integer * iprint,
  integer *isumm, integer *inform, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  double *rw, integer *lenrw, ftnlen buffer_len,
  ftnlen cw_len );

 int snspec_
( integer *ispecs, integer *inform, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  double *rw, integer *lenrw, ftnlen cw_len);

 int snmema_
( integer *iExit, integer *nef, integer *n,
  integer *nxname, integer *nfname,
  integer *nea, integer *neg,
  integer *mincw, integer *miniw, integer *minrw,
  char *cw, integer *lencw,
  integer *iw, integer *leniw,
  doublereal *rw, integer *lenrw,
  ftnlen cw_len);

 int snjac_
( integer *iExit, integer *nef, integer *n, U_fp userfg,
  integer *iafun, integer *javar, integer *lena, integer *nea, doublereal *a,
  integer *igfun, integer *jgvar, integer *leng, integer *neg,
  doublereal *x, doublereal *xlow, doublereal *xupp,
  integer *mincw, integer *miniw, integer *minrw,
  char *cu, integer *lencu, integer *iu, integer *leniu, double *ru, integer *lenru,
  char *cw, integer *lencw, integer *iw, integer *leniw, double *rw, integer *lenrw,
  ftnlen cu_len, ftnlen cw_len);

 // from filewrapper.h
 int snopenappend_
 ( integer *iunit, char *name, integer *inform, ftnlen name_len);

 int snfilewrapper_
 ( char *name__, integer *ispec, integer *
   inform__, char *cw, integer *lencw, integer *iw,
   integer *leniw, doublereal *rw, integer *lenrw,
   ftnlen name_len, ftnlen cw_len);

 int snclose_
 ( integer *iunit);

 int snopen_
 ( integer *iunit);

#undef real
} // extern "C"

#endif
