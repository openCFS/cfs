/* C interface for snopt */

#ifndef SNOPTINTERFACE_HH_
#define SNOPTINTERFACE_HH_


#include <stdio.h>
#include <stdint.h>

#include <def_cfs_fortran_interface.hh>

extern "C" {

  // callback for snopt to request function and gradient evaluation
  typedef int (*My_fp)( int32_t *Status, int32_t *n,
       double x[],     int32_t *needF, int32_t *neF,  double F[],
       int32_t    *needG,  int32_t *neG,  double G[],
       char       *cu,     int32_t *lencu,
       int32_t    iu[],    int32_t *leniu,
       double ru[],    int32_t *lenru);

  // patched callback for snopt to inform about increase of the major iteration
  typedef int (*My_usrmjr)(int32_t* nMajor );

  void snopta
     ( int32_t *start, int32_t *nf, int32_t *n,
       int32_t *nxname, int32_t *nfname, double *objadd, const int32_t *objrow,
       char *prob, My_fp usrfun, My_usrmjr usrmjr, int32_t *iafun, int32_t *javar,
       int32_t *lena, int32_t *nea, double *a, int32_t *igfun,
       int32_t *jgvar, int32_t *leng, int32_t *neg, double *xlow,
       double *xupp, char *xnames, double *flow, double *fupp,
       char *fnames, double *x, int32_t *xstate, double *xmul,
       double *f, int32_t *fstate, double *fmul, int32_t *inform__,
       int32_t *mincw, int32_t *miniw, int32_t *minrw, int32_t *ns,
       int32_t *ninf, double *sinf, char *cu, int32_t *lencu, int32_t *iu,
       int32_t *leniu, double *ru, int32_t *lenru, char *cw, int32_t *lencw,
       int32_t *iw, int32_t *leniw, double *rw, int32_t *lenrw,
       intptr_t prob_len, intptr_t xnames_len, intptr_t fnames_len, intptr_t cu_len,
       intptr_t cw_len);

  // it is unclear, where the intptr_t cw_len argument comes from.
  // https://web.stanford.edu/group/SOL/guides/sndoc7.pdf documents 8 parameters ?!
  void sninit
     ( int32_t *iPrint, int32_t *iSumm, char *cw,
       int32_t *lencw, int32_t *iw, int32_t *leniw,
       double *rw, int32_t *lenrw, intptr_t cw_len );

  void sngeti
     ( char *buffer, int32_t *ivalue, int32_t *inform__,
       char *cw, int32_t *lencw, int32_t *iw,
       int32_t *leniw, double *rw, int32_t *lenrw,
       intptr_t buffer_len, intptr_t cw_len);


  void sngetr
     ( char *buffer, double *ivalue, int32_t *inform__,
       char *cw, int32_t *lencw, int32_t *iw,
       int32_t *leniw, double *rw, int32_t *lenrw,
       intptr_t buffer_len, intptr_t cw_len);

  void snset
     ( const char *buffer, int32_t *iprint, int32_t *isumm,
       int32_t *inform__, char *cw, int32_t *lencw,
       int32_t *iw, int32_t *leniw,
       double *rw, int32_t *lenrw,
       intptr_t buffer_len, intptr_t cw_len);

  void sngetc
     ( char *buffer, char *ivalue, int32_t *inform__,
       char *cw, int32_t *lencw, int32_t *iw,
       int32_t *leniw, double *rw, int32_t *lenrw,
       intptr_t buffer_len, intptr_t ivalue_len, intptr_t cw_len);

  void snseti
     ( const char *buffer, int32_t *ivalue, int32_t *iprint,
       int32_t *isumm, int32_t *inform__, char *cw,
       int32_t *lencw, int32_t *iw, int32_t *leniw,
       double *rw, int32_t *lenrw, intptr_t buffer_len,
       intptr_t cw_len);

  void snsetr
     ( const char *buffer, double *rvalue, int32_t * iprint,
       int32_t *isumm, int32_t *inform__, char *cw,
       int32_t *lencw, int32_t *iw, int32_t *leniw,
       double *rw, int32_t *lenrw, intptr_t buffer_len,
       intptr_t cw_len);

  void snspec
     ( int32_t *ispecs, int32_t *inform__, char *cw,
       int32_t *lencw, int32_t *iw, int32_t *leniw,
       double *rw, int32_t *lenrw, intptr_t cw_len);

  void snmema
     ( int32_t *iexit, int32_t *nf, int32_t *n, int32_t *nxname,
       int32_t *nfname, int32_t *nea, int32_t *neg,
       int32_t *mincw, int32_t *miniw,
       int32_t *minrw, char *cw, int32_t *lencw, int32_t *iw,
       int32_t *leniw, double *rw, int32_t *lenrw,
       intptr_t cw_len);


  void snjac
     ( int32_t *inform__, int32_t *nf, int32_t *n, My_fp userfg,
       int32_t *iafun, int32_t *javar, int32_t *lena,
       int32_t *nea, double *a, int32_t *igfun,
       int32_t *jgvar, int32_t *leng, int32_t *neg,
       double *x, double *xlow, double *xupp,
       int32_t *mincw, int32_t *miniw,
       int32_t *minrw, char *cu, int32_t *lencu,
       int32_t *iu, int32_t *leniu, double *ru,
       int32_t *lenru, char *cw, int32_t *lencw, int32_t *iw,
       int32_t *leniw, double *rw, int32_t *lenrw,
       intptr_t cu_len, intptr_t cw_len );

  // from filewrapper.h
  // function for opening snopt output files
  int snopenappend (int32_t *iunit, char *name, int32_t *inform, intptr_t name_len);

  int snclose (int32_t *iunit);
} // extern "C"

#endif
