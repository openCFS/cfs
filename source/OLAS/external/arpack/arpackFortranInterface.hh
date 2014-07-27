#ifndef OLAS_ARPACK_FORTRAN_INTERFACE_HH
#define OLAS_ARPACK_FORTRAN_INTERFACE_HH

#include "General/defs.hh"

#include <def_cfs_fortran_interface.hh>

using CoupledField::Integer;
using CoupledField::Double;
using CoupledField::Complex;

extern "C" {
  //! arpack debug common block
  extern struct {
    Integer logfil, ndigit, mgetv0;
    Integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    Integer mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
    Integer mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
  } debug;

  //! arpack driver routine controling calculation
  void dsaupd(Integer *ido, char *bmat,
              Integer *n, char *which,
              Integer *nev, Double *tol,
              Double *resid, Integer *ncv,
              Double *V, Integer *ldv,
              Integer *iparam, Integer *ipntr,
              Double *workd, Double *workl,
              Integer *lworkl, Integer *info);

  //! arpack driver routine controling eigenvalue and eigenvector extraction
  void dseupd(bool *rvec, char *howMny,
              Double *select, Double *d,
              Double *z, Integer *ldz,
              Double *shift, char *bmat,
              Integer* size, char *which,
              Integer *nev, Double *tol,
              Double *resid, Integer *ncv,
              Double *V, Integer *ldv,
              Integer *iparam, Integer *ipntr,
              Double *workd, Double *workl,
              Integer *lworkl, Integer *info);
  
  
  void znaupd(Integer *ido, char *bmat,
              Integer *n, char *which,
              Integer *nev, Double *tol,
              Complex *resid, Integer *ncv,
              Complex *V, Integer *ldv,
              Integer *iparam, Integer *ipntr,
              Complex *workd, Complex *workl,
              Integer *lworkl, 
              Double *workDbleD, Integer *info);

  void zneupd(bool *rvec, char *howMny,
              Double *select, Complex *d,
              Complex *z, Integer *ldz,
              Complex *shift, Complex *zwork, 
              char *bmat,
              Integer* size, char *which,
              Integer *nev, Double *tol,
              Complex *resid, Integer *ncv,
              Complex *V, Integer *ldv,
              Integer *iparam, Integer *ipntr,
              Complex *workd, Complex *workl,
              Integer *lworkl, 
              Double  *workDbleD,Integer *info);

}
#endif
