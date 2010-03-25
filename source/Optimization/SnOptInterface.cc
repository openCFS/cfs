#include "Optimization/SnOptInterface.hh"

extern "C"
{
  
typedef short flag;
typedef short ftnint;

/*open*/
typedef struct
{ flag oerr;
  ftnint ounit;
  char *ofnm;
  ftnlen ofnmlen;
  char *osta;
  char *oacc;
  char *ofm;
  ftnint orl;
  char *oblnk;
} olist;

/*close*/
typedef struct
{ flag cerr;
  ftnint cunit;
  char *csta;
} cllist;

/* Subroutine */ int snfilewrapper_(char *name__, integer *ispec, integer *
  inform__, char *cw, integer *lencw, integer *iw, integer *leniw,
  double *rw, integer *lenrw, ftnlen name_len, ftnlen cw_len)
{
    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int snspec_(integer *, integer *, char *, integer
      *, integer *, integer *, double *, integer *, ftnlen);
    static integer iostat;

/*     ================================================================== */
/*     Read options for snopt from the file named name. inform .eq 0 if */
/*     successful. */

/*     09 Jan 2000: First version of snFilewrapper. */
/*     ================================================================== */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    o__1.oerr = 1;
    o__1.ounit = *ispec;
    o__1.ofnmlen = name_len;
    o__1.ofnm = name__;
    o__1.orl = 0;
    o__1.osta = const_cast<char*>("old");
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    iostat = f_open(&o__1);
    if (iostat != 0) {
  *inform__ = iostat + 2;
    } else {
  snspec_(ispec, inform__, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
    (ftnlen)8);
  cl__1.cerr = 0;
  cl__1.cunit = *ispec;
  cl__1.csta = 0;
  f_clos(&cl__1);
    }
    return 0;
} /* snfilewrapper_ */

/* Subroutine */ int snopenappend_(integer *iunit, char *name__, integer *
  inform__, ftnlen name_len)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);

/*     ================================================================== */
/*     Open file named name to Fortran unit iunit. inform .eq. 0 if */
/*     sucessful.  Although opening for appending is not in the f77 */
/*     standard, it is understood by f2c. */

/*     09 Jan 2000: First version of snOpenappend */
/*     ================================================================== */
    o__1.oerr = 1;
    o__1.ounit = *iunit;
    o__1.ofnmlen = name_len;
    o__1.ofnm = name__;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = const_cast<char*>("append");
    o__1.ofm = 0;
    o__1.oblnk = 0;
    *inform__ = f_open(&o__1);
    return 0;
} /* snopenappend_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snOpenappend */
/* Subroutine */ int snclose_(integer *iunit)
{
    /* System generated locals */
    cllist cl__1;

    /* Builtin functions */
    integer f_clos(cllist *);

/*     ================================================================== */
/*     Close unit iunit. */

/*     09 Jan 2000: First version of snClose */
/*     ================================================================== */
    cl__1.cerr = 0;
    cl__1.cunit = *iunit;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* snclose_ */


/* Subroutine */ int snopen_(integer *iunit)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *);

    /* Local variables */
    static char lfile[20];

/*     ================================================================= */
/*     Open file named name to Fortran unit iunit.  inform .eq. 0 */
/*     if successful. */
/*     ================================================================= */
    s_copy(lfile, const_cast<char*>("testing.out"), (ftnlen)20, (ftnlen)11);
    o__1.oerr = 0;
    o__1.ounit = *iunit;
    o__1.ofnmlen = 20;
    o__1.ofnm = lfile;
    o__1.orl = 0;
    o__1.osta = const_cast<char*>("UNKNOWN");
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    return 0;
} /* snopen_ */

}
