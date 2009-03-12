// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_ARPACK_SETTINGS_HH
#define OLAS_ARPACK_SETTINGS_HH

#include "General/defs.hh"

using CoupledField::Integer;
using CoupledField::Double;

//! arpack driver routine controling calculation
#define ARPACK_DSAUPD dsaupd_
//! arpack driver routine controling eigenvalue and eigenvector extraction
#define ARPACK_DSEUPD dseupd_
//! arpack debug common block
#define ARPACK_DEBUG_CMN debug_

extern "C" {
  void ARPACK_DSAUPD(Integer *ido, char *bmat, Integer *n, char *which,
                   Integer *nev, Double *tol, Double *resid,
                   Integer *ncv, Double *V, Integer *ldv,
                   Integer *iparam, Integer *ipntr, Double *workd,
                   Double *workl, Integer *lworkl, Integer *info);

  void ARPACK_DSEUPD(bool *rvec, char *howMny, Double *select, Double *d,
                   Double *z, Integer *ldz, Double *shift, char *bmat,
                   Integer* size, char *which, Integer *nev, Double *tol, Double *resid,
                   Integer *ncv, Double *V, Integer *ldv,
                   Integer *iparam, Integer *ipntr, Double *workd,
                   Double *workl, Integer *lworkl, Integer *info);

  // arpack debug "common" statement.

  extern struct {
    Integer logfil, ndigit, mgetv0;
    Integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    Integer mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
    Integer mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
  } ARPACK_DEBUG_CMN;

}
#endif
