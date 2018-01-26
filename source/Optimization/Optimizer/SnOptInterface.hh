/* C interface for snopt */

#ifndef SNOPTINTERFACE_HH_
#define SNOPTINTERFACE_HH_


#include <stdio.h>
#include <stdint.h>

typedef intptr_t ftnlen;
typedef int32_t integer;
typedef double doublereal;

#ifdef USE_ADOLC

  typedef int (*U_fp)(...);
  typedef int (*My_fp)( integer *Status, integer *n,
       Double x[],     integer *needF, integer *neF,  Double F[],
       integer    *needG,  integer *neG,  Double G[],
       char       *cu,     integer *lencu,
       integer    iu[],    integer *leniu,
       Double ru[],    integer *lenru )

  void snopta_
     ( integer *start, integer *nf, integer *n,
       integer *nxname, integer *nfname, Double *objadd, const integer *objrow,
       char *prob, My_fp usrfun, integer *iafun, integer *javar,
       integer *lena, integer *nea, Double *a, integer *igfun,
       integer *jgvar, integer *leng, integer *neg, Double *xlow,
       Double *xupp, char *xnames, Double *flow, Double *fupp,
       char *fnames, Double *x, integer *xstate, Double *xmul,
       Double *f, integer *fstate, Double *fmul, integer *inform__,
       integer *mincw, integer *miniw, integer *minrw, integer *ns,
       integer *ninf, Double *sinf, char *cu, integer *lencu, integer *iu,
       integer *leniu, Double *ru, integer *lenru, char *cw, integer *lencw,
       integer *iw, integer *leniw, Double *rw, integer *lenrw,
       ftnlen prob_len, ftnlen xnames_len, ftnlen fnames_len, ftnlen cu_len,
       ftnlen cw_len){EXCEPTION("Not implemented.");}

  void sninit_
     ( integer *iPrint, integer *iSumm, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       Double *rw, integer *lenrw, ftnlen cw_len ){EXCEPTION("Not implemented.");}

  void sngeti_
     ( char *buffer, integer *ivalue, integer *inform__,
       char *cw, integer *lencw, integer *iw,
       integer *leniw, Double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen cw_len){EXCEPTION("Not implemented.");}


  void sngetr_
     ( char *buffer, Double *ivalue, integer *inform__,
       char *cw, integer *lencw, integer *iw,
       integer *leniw, Double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen cw_len){EXCEPTION("Not implemented.");}

  void snset_
     ( const char *buffer, integer *iprint, integer *isumm,
       integer *inform__, char *cw, integer *lencw,
       integer *iw, integer *leniw,
       Double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen cw_len){EXCEPTION("Not implemented.");}

  void sngetc_
     ( char *buffer, char *ivalue, integer *inform__,
       char *cw, integer *lencw, integer *iw,
       integer *leniw, Double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen ivalue_len, ftnlen cw_len){EXCEPTION("Not implemented.");}

  void snseti_
     ( const char *buffer, integer *ivalue, integer *iprint,
       integer *isumm, integer *inform__, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       Double *rw, integer *lenrw, ftnlen buffer_len,
       ftnlen cw_len){EXCEPTION("Not implemented.");}

  void snsetr_
     ( const char *buffer, Double *rvalue, integer * iprint,
       integer *isumm, integer *inform__, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       Double *rw, integer *lenrw, ftnlen buffer_len,
       ftnlen cw_len){EXCEPTION("Not implemented.");}

  void snspec_
     ( integer *ispecs, integer *inform__, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       Double *rw, integer *lenrw, ftnlen cw_len){EXCEPTION("Not implemented.");}

  void snmema_
     ( integer *iexit, integer *nf, integer *n, integer *nxname,
       integer *nfname, integer *nea, integer *neg,
       integer *mincw, integer *miniw,
       integer *minrw, char *cw, integer *lencw, integer *iw,
       integer *leniw, Double *rw, integer *lenrw,
       ftnlen cw_len){EXCEPTION("Not implemented.");}


  void snjac_
     ( integer *inform__, integer *nf, integer *n, My_fp userfg,
       integer *iafun, integer *javar, integer *lena,
       integer *nea, Double *a, integer *igfun,
       integer *jgvar, integer *leng, integer *neg,
       Double *x, Double *xlow, Double *xupp,
       integer *mincw, integer *miniw,
       integer *minrw, char *cu, integer *lencu,
       integer *iu, integer *leniu, Double *ru,
       integer *lenru, char *cw, integer *lencw, integer *iw,
       integer *leniw, Double *rw, integer *lenrw,
       ftnlen cu_len, ftnlen cw_len ){EXCEPTION("Not implemented.");}

  // from filewrapper.h
  // function for opening snopt output files
  int snopenappend_ 
  (integer *iunit, char *name, integer *inform, ftnlen name_len);

  int snclose_
  (integer *iunit);


#else

extern "C" {

  typedef int (*U_fp)(...);
  typedef int (*My_fp)( integer *Status, integer *n,
       double x[],     integer *needF, integer *neF,  double F[],
       integer    *needG,  integer *neG,  double G[],
       char       *cu,     integer *lencu,
       integer    iu[],    integer *leniu,
       double ru[],    integer *lenru );

  void snopta_
     ( integer *start, integer *nf, integer *n,
       integer *nxname, integer *nfname, double *objadd, const integer *objrow,
       char *prob, My_fp usrfun, integer *iafun, integer *javar,
       integer *lena, integer *nea, double *a, integer *igfun,
       integer *jgvar, integer *leng, integer *neg, double *xlow,
       double *xupp, char *xnames, double *flow, double *fupp,
       char *fnames, double *x, integer *xstate, double *xmul,
       double *f, integer *fstate, double *fmul, integer *inform__,
       integer *mincw, integer *miniw, integer *minrw, integer *ns,
       integer *ninf, double *sinf, char *cu, integer *lencu, integer *iu,
       integer *leniu, double *ru, integer *lenru, char *cw, integer *lencw,
       integer *iw, integer *leniw, double *rw, integer *lenrw,
       ftnlen prob_len, ftnlen xnames_len, ftnlen fnames_len, ftnlen cu_len,
       ftnlen cw_len);

  void sninit_
     ( integer *iPrint, integer *iSumm, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       double *rw, integer *lenrw, ftnlen cw_len );

  void sngeti_
     ( char *buffer, integer *ivalue, integer *inform__,
       char *cw, integer *lencw, integer *iw,
       integer *leniw, double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen cw_len);


  void sngetr_
     ( char *buffer, double *ivalue, integer *inform__,
       char *cw, integer *lencw, integer *iw,
       integer *leniw, double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen cw_len);

  void snset_
     ( const char *buffer, integer *iprint, integer *isumm,
       integer *inform__, char *cw, integer *lencw,
       integer *iw, integer *leniw,
       double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen cw_len);

  void sngetc_
     ( char *buffer, char *ivalue, integer *inform__,
       char *cw, integer *lencw, integer *iw,
       integer *leniw, double *rw, integer *lenrw,
       ftnlen buffer_len, ftnlen ivalue_len, ftnlen cw_len);

  void snseti_
     ( const char *buffer, integer *ivalue, integer *iprint,
       integer *isumm, integer *inform__, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       double *rw, integer *lenrw, ftnlen buffer_len,
       ftnlen cw_len);

  void snsetr_
     ( const char *buffer, double *rvalue, integer * iprint,
       integer *isumm, integer *inform__, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       double *rw, integer *lenrw, ftnlen buffer_len,
       ftnlen cw_len);

  void snspec_
     ( integer *ispecs, integer *inform__, char *cw,
       integer *lencw, integer *iw, integer *leniw,
       double *rw, integer *lenrw, ftnlen cw_len);

  void snmema_
     ( integer *iexit, integer *nf, integer *n, integer *nxname,
       integer *nfname, integer *nea, integer *neg,
       integer *mincw, integer *miniw,
       integer *minrw, char *cw, integer *lencw, integer *iw,
       integer *leniw, double *rw, integer *lenrw,
       ftnlen cw_len);


  void snjac_
     ( integer *inform__, integer *nf, integer *n, My_fp userfg,
       integer *iafun, integer *javar, integer *lena,
       integer *nea, double *a, integer *igfun,
       integer *jgvar, integer *leng, integer *neg,
       double *x, double *xlow, double *xupp,
       integer *mincw, integer *miniw,
       integer *minrw, char *cu, integer *lencu,
       integer *iu, integer *leniu, double *ru,
       integer *lenru, char *cw, integer *lencw, integer *iw,
       integer *leniw, double *rw, integer *lenrw,
       ftnlen cu_len, ftnlen cw_len );

  // from filewrapper.h
  // function for opening snopt output files
  int snopenappend_ 
  (integer *iunit, char *name, integer *inform, ftnlen name_len);

  int snclose_
  (integer *iunit);

#undef real
} // extern "C"

#endif
#endif
