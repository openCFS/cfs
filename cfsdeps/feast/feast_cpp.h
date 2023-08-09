#ifndef FEAST_RCI_CPP_HEADER
#define FEAST_RCI_CPP_HEADER

#ifdef __cplusplus

#include <complex>

#ifdef WIN32  // Windows name mangling
# define feastinit    FEASTINIT
// RCI
# define dfeast_srci  DFEAST_SRCI
# define dfeast_srcix  DFEAST_SRCIX
# define dfeast_grci  DFEAST_GRCI
# define dfeast_grcix  DFEAST_GRCIX
# define zfeast_srci  ZFEAST_SRCI
# define zfeast_grci  ZFEAST_GRCI
# define zfeast_grcix  ZFEAST_GRCIX
# define zfeast_hrci  ZFEAST_HRCI
# define zfeast_hrcix  ZFEAST_HRCIX
// contours (returns integration nodes and weights)
# define zfeast_customcontour ZFEAST_CUSTOMCONTOUR
# define zfeast_contour ZFEAST_CONTOUR
# define zfeast_gcontour ZFEAST_GCONTOUR
// solvers (normal mode, replaced by expert mode!)
# define dfeast_scsrev DFEAST_SCSREV
# define dfeast_scsrgv DFEAST_SCSRGV
# define dfeast_gcsrev DFEAST_GCSREV
# define dfeast_gcsrgv DFEAST_GCSRGV
# define zfeast_scsrev DFEAST_SCSREV
# define zfeast_scsrgv ZFEAST_SCSRGV
# define zfeast_gcsrev ZFEAST_GCSREV
# define zfeast_gcsrgv ZFEAST_GCSRGV
# define zfeast_hcsrev ZFEAST_HCSREV
# define zfeast_hcsrgv ZFEAST_HCSRGV
# define dfeast_scsrpev DFEAST_SCSRPEV
# define zfeast_scsrpev ZFEAST_SCSRPEV
# define dfeast_gcsrpev DFEAST_GCSRPEV
# define zfeast_gcsrpev ZFEAST_GCSRPEV
// solvers (expert mode)
# define dfeast_scsrevx DFEAST_SCSREVX
# define dfeast_scsrgvx DFEAST_SCSRGVX
# define dfeast_gcsrevx DFEAST_GCSREVX
# define dfeast_gcsrgvx DFEAST_GCSRGVX
# define zfeast_scsrevx DFEAST_SCSREVX
# define zfeast_scsrgvx ZFEAST_SCSRGVX
# define zfeast_gcsrevx ZFEAST_GCSREVX
# define zfeast_gcsrgvx ZFEAST_GCSRGVX
# define zfeast_hcsrevx ZFEAST_HCSREVX
# define zfeast_hcsrgvx ZFEAST_HCSRGVX
# define dfeast_scsrpevx DFEAST_SCSRPEVX
# define zfeast_scsrpevx ZFEAST_SCSRPEVX
# define dfeast_gcsrpevx DFEAST_GCSRPEVX
# define zfeast_gcsrpevx ZFEAST_GCSRPEVX
#else  // UNIX name mangling
# define feastinit    feastinit_
// RCI
# define dfeast_srci  dfeast_srci_
# define dfeast_srcix  dfeast_srcix_
# define dfeast_grci  dfeast_grci_
# define dfeast_grcix  dfeast_grcix_
# define zfeast_srci  zfeast_srci_
# define zfeast_srcix  zfeast_srcix_
# define zfeast_grci  zfeast_grci_
# define zfeast_grcix  zfeast_grcix_
# define zfeast_hrci  zfeast_hrci_
# define zfeast_hrcix  zfeast_hrcix_
// contours (returns integration nodes and weights)
# define zfeast_customcontour zfeast_customcontour_
# define zfeast_contour zfeast_contour_
# define zfeast_gcontour zfeast_gcontour_
// solvers (normal mode, replaced by expert mode!)
# define dfeast_gcsrev dfeast_gcsrev_
# define dfeast_gcsrgv dfeast_gcsrgv_
# define zfeast_scsrev zfeast_scsrev_
# define zfeast_scsrgv zfeast_scsrgv_
# define zfeast_gcsrev zfeast_gcsrev_
# define zfeast_gcsrgv zfeast_gcsrgv_
# define zfeast_hcsrev zfeast_hcsrev_
# define zfeast_hcsrgv zfeast_hcsrgv_
# define dfeast_scsrev dfeast_scsrev_
# define dfeast_scsrgv dfeast_scsrgv_
# define dfeast_scsrpev dfeast_scsrpev_
# define zfeast_hcsrpev zfeast_hcsrpev_
# define zfeast_scsrpev zfeast_scsrpev_
# define dfeast_gcsrpev dfeast_gcsrpev_
# define zfeast_gcsrpev zfeast_gcsrpev_
// solvers (expert mode)
# define dfeast_gcsrevx dfeast_gcsrevx_
# define dfeast_gcsrgvx dfeast_gcsrgvx_
# define zfeast_scsrevx zfeast_scsrevx_
# define zfeast_scsrgvx zfeast_scsrgvx_
# define zfeast_gcsrevx zfeast_gcsrevx_
# define zfeast_gcsrgvx zfeast_gcsrgvx_
# define zfeast_hcsrevx zfeast_hcsrevx_
# define zfeast_hcsrgvx zfeast_hcsrgvx_
# define dfeast_scsrevx dfeast_scsrevx_
# define dfeast_scsrgvx dfeast_scsrgvx_
# define dfeast_scsrpevx dfeast_scsrpevx_
# define zfeast_hcsrpevx zfeast_hcsrpevx_
# define zfeast_scsrpevx zfeast_scsrpevx_
# define dfeast_gcsrpevx dfeast_gcsrpevx_
# define zfeast_gcsrpevx zfeast_gcsrpevx_

#endif

extern "C" {
  // Initialization
  void feastinit(int *feastparam);

  // search contour
  void zfeast_contour(double *Emin, double *Emax, int *fpm2, int *fpm16, int *fpm18, double *Zne, double *Wne);

  void zfeast_gcontour(double *Emid, double *r, int *fpm8, int *fpm16, int *fpm18, int *fpm19, double *Zne, double *Wne);

  void zfeast_customcontour(int *numNC, int *numP, int *Nedge, int *Tedge, double *Zedge, double *Zne, double *Wne);

  // Double precision reverse communication interface
//  void dfeast_srci( int *ijob, int *N, std::complex<double> *Ze,                double *work1, std::complex<double> *zwork2, double *Aq,  double *Sq,  int *fpm, double *epsout, int *loop,               double *Emin, double *Emax, int *M0,               double *E, double *X,               int *M, double *res, int *info);
//  void zfeast_hrci( int *ijob, int *N, std::complex<double> *Ze, std::complex<double> *zwork1, std::complex<double> *zwork2, double *zAq, double *zSq, int *fpm, double *epsout, int *loop,               double *Emin, double *Emax, int *M0,               double *E, std::complex<double> *X, int *M, double *res, int *info);
//  void zfeast_srci( int *ijob, int *N, std::complex<double> *Ze, std::complex<double> *zwork1, std::complex<double> *zwork2, double *zAq, double *zSq, int *fpm, double *epsout, int *loop, std::complex<double> *Emid, double *r,    int *M0, std::complex<double> *E, std::complex<double> *X, int *M, double *res, int *info);
//  void dfeast_grci( int *ijob, int *N, std::complex<double> *Ze, std::complex<double> *zwork1, std::complex<double> *zwork2, double *zAq, double *zSq, int *fpm, double *epsout, int *loop, std::complex<double> *Emid, double *r,    int *M0, std::complex<double> *E, std::complex<double> *X, int *M, double *res, int *info);
//  void zfeast_grci( int *ijob, int *N, std::complex<double> *Ze, std::complex<double> *zwork1, std::complex<double> *zwork2, double *zAq, double *zSq, int *fpm, double *epsout, int *loop, std::complex<double> *Emid, double *r,    int *M0, std::complex<double> *E, std::complex<double> *X, int *M, double *res, int *info);
//  void zfeast_grcix(int *ijob, int *N, std::complex<double> *Ze, std::complex<double> *zwork1, std::complex<double> *zwork2, double *zAq, double *zSq, int *fpm, double *epsout, int *loop, std::complex<double> *Emid, double *r,    int *M0, std::complex<double> *E, std::complex<double> *X, int *M, double *res, int *info, std::complex<double> *Zne, std::complex<double> *Wne);
//  void zfeast_customcontour(int *numNC, int *numP, int *Nedge, int *Tedge, std::complex<double> *Zedge, std::complex<double> *Zne, std::complex<double> *Wne);
  void dfeast_srci(char *ijob,int *N,double *Ze,double *work,double   *workc,double   *Aq,double   *Sq,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void zfeast_hrci(char *ijob,int *N,double   *Ze,double   *work,double   *workc,double   *zAq,double   *zSq,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double   *q,int *mode,double *res,int *info);

  void dfeast_srcix(char *ijob,int *N,double *Ze,double *work,double   *workc,double   *Aq,double   *Sq,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);

  void zfeast_hrcix(char *ijob,int *N,double   *Ze,double   *work,double   *workc,double   *zAq,double   *zSq,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double   *q,int *mode,double *res,int *info,double *Zne, double *Wne);

  void zfeast_srci(char *ijob,int *N,double *Ze,double *workr, double *workl ,double *workc,double *Aq,double *Sq,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,int *mode,double *resr,int *info);

  void dfeast_grci(char *ijob,int *N,double *Ze,double *workr, double *workl ,double *workc,double *Aq,double *Sq,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr, double *ql,int *mode,double *resr,double *resl, int *info);

  void zfeast_grci(char *ijob,int *N,double *Ze,double *workr, double *workl ,double *workc,double *Aq,double *Sq,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,double *ql,int *mode,double *resr,double *resl,int *info);

  void zfeast_srcix(char *ijob,int *N,double *Ze,double *workr, double *workl ,double *workc,double *Aq,double *Sq,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,int *mode,double *resr,int *info,double *Zne, double *Wne);

  void dfeast_grcix(char *ijob,int *N,double *Ze,double *workr, double *workl ,double *workc,double *Aq,double *Sq,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,double *ql,int *mode,double *resr,double *resl,int *info, double *Zne, double *Wne);

  void zfeast_grcix(char *ijob,int *N,double *Ze,double *workr, double *workl ,double *workc,double *Aq,double *Sq,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,double *ql,int *mode,double *resr,double *resl,int *info, double *Zne, double *Wne);

//////////////////////////////////////////////////////////////

  void dfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void dfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void dfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void dfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void zfeast_hcsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void zfeast_hcsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void zfeast_hcsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void zfeast_hcsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void zfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void zfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void zfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void zfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void zfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void zfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void zfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void zfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void dfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void dfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

  void dfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);

  void dfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
  
  void dfeast_scsrpevx(char *UPLO,int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info,double *Zne,double *Wne );
  
  void dfeast_scsrpev(char *UPLO,int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info);

  void zfeast_hcsrpevx(char *UPLO,int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info, double *Zne,double *Wne );

  void zfeast_hcsrpev(char *UPLO,int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info);

  void zfeast_scsrpevx(char *UPLO,int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info, double *Zne,double *Wne );

  void zfeast_scsrpev(char *UPLO,int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info);
   
  void dfeast_gcsrpevx(int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info,double *Zne,double *Wne );

  void dfeast_gcsrpev(int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info);

  void zfeast_gcsrpevx(int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info,double *Zne,double *Wne );

  void zfeast_gcsrpev(int *dmax,int* N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int* loop,double* Emid,double* r,int *M0,double *lambda, double* q,int *mode,double *res,int* info);
} 

#endif

#endif
