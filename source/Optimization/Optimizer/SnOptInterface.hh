/* C interface for snopt */

#ifndef SNOPTINTERFACE_HH_
#define SNOPTINTERFACE_HH_


#include <stdio.h>
#include <stdint.h>

#include <def_cfs_fortran_interface.hh>

extern "C" {

  /** we use here a now really documented callback-mechanism by Snopt. Opus found it in the Fortran sources,
      it is not in the manuals. */

  // callback for snopt to request function and gradient evaluation
  typedef int (*My_fp)( int32_t *Status, int32_t *n,
       double x[],     int32_t *needF, int32_t *neF,  double F[],
       int32_t    *needG,  int32_t *neG,  double G[],
       char       *cu,     int32_t *lencu,
       int32_t    iu[],    int32_t *leniu,
       double ru[],    int32_t *lenru);

  // callback which snopt calls once per major iteration, including the initial (nMajor == 0)
  // and the final one. Corresponds to isnSTOP of the official snopt.h. Set iAbort != 0 to
  // terminate the run. Only usable via snkera, snopta always uses snopt's dummy version.
  typedef void (*My_snSTOP)
     ( int32_t *iAbort, int32_t KTcond[], int32_t *mjrPrtlvl, int32_t *minimize,
       int32_t *m, int32_t *maxS, int32_t *n, int32_t *nb,
       int32_t *nnCon0, int32_t *nnCon, int32_t *nnObj0, int32_t *nnObj, int32_t *nS,
       int32_t *itn, int32_t *nMajor, int32_t *nMinor, int32_t *nSwap,
       double *condZHZ, int32_t *iObj, double *scaleObj, double *objAdd,
       double *fObj, double *fMerit, double penParm[], double *step,
       double *primalInf, double *dualInf, double *maxVi, double *maxViRel, int32_t hs[],
       int32_t *neJ, int32_t *nlocJ, int32_t locJ[], int32_t indJ[], double Jcol[],
       int32_t *negCon, double scales[], double bl[], double bu[],
       double Fx[], double fCon[], double gCon[], double gObj[],
       double yCon[], double pi[], double rc[], double rg[], double x[],
       char cu[], int32_t *lencu, int32_t iu[], int32_t *leniu, double ru[], int32_t *lenru,
       char cw[], int32_t *lencw, int32_t iw[], int32_t *leniw, double rw[], int32_t *lenrw);

  // snopt's own log routines. We never call them, we only hand them over to snkera, hence an
  // empty parameter list is sufficient.
  typedef void (*My_snLogDummy)();
  void snlog();
  void snlog2();
  void sqlog();

  void snopta
     ( int32_t *start, int32_t *nf, int32_t *n,
       int32_t *nxname, int32_t *nfname, double *objadd, const int32_t *objrow,
       char *prob, My_fp usrfun, int32_t *iafun, int32_t *javar,
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

  // the kernel of snopta. Identical apart from the four callbacks, which snopta fills with
  // snopt's default versions. This is the documented way to supply an own snSTOP.
  void snkera
     ( int32_t *start, int32_t *nf, int32_t *n,
       int32_t *nxname, int32_t *nfname, double *objadd, const int32_t *objrow,
       char *prob, My_fp usrfun, My_snLogDummy snlog, My_snLogDummy snlog2,
       My_snLogDummy sqlog, My_snSTOP snstop, int32_t *iafun, int32_t *javar,
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

  // the trailing cw_len is not a snopt argument but the hidden length which Fortran appends for
  // every character argument. cw is declared as character cw(lencw)*8, hence cw_len = 8*lencw.
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
