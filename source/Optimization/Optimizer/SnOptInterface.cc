#include "Optimization/Optimizer/SnOptInterface.hh"

extern "C"
{

#if defined(__x86_64__) || defined(__ia64__) || defined(__sparc64__) || defined(__alpha)
  typedef int flag;
  typedef int ftnlen;
  typedef int ftnint; 
#else
  typedef short flag;
  typedef short ftnlen;
  typedef short ftnint;
#endif

/*open*/
typedef struct
{ 
  flag oerr;
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
{ 
  flag cerr;
  ftnint cunit;
  char *csta;
} cllist;


/*     ================================================================== */
/*     Open file named name to Fortran unit iunit. inform .eq. 0 if */
/*     successful.  Although opening for appending is not in the f77 */
/*     standard, it is understood by f2c. */

/*     09 Jan 2000: First version of snOpenappend */
/*     ================================================================== */
int snopenappend_(int *iunit, char *name, int *inform, int name_len)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    int f_open(olist *); 

    char action[] = "append";
    o__1.oerr = 1;
    o__1.ounit = *iunit;
    o__1.ofnmlen = name_len;
    o__1.ofnm = name;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = action;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    *inform = f_open(&o__1);
    return 0;
} /* snopenappend_ */


/*     ================================================================== */
/*     Close unit iunit. */

/*     09 Jan 2000: First version of snClose */
/*     ================================================================== */
int snclose_(int *iunit)
{
    /* System generated locals */
    cllist cl__1;

    /* Builtin functions */
    int f_clos(cllist *);

    cl__1.cerr = 0;
    cl__1.cunit = *iunit;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* snclose_ */


}
