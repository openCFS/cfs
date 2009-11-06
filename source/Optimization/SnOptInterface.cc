#include "Optimization/SnOptInterface.hh"

extern "C"
{
  
int snoptclose_(int *iunit)
{
  /* System generated locals */
  cllist cl__1;

  /* Builtin functions */
  int f_clos(cllist *);

  /*     ================================================================== */
  /*     Close unit iunit. */
  /*     09 Jan 2000: First version of snoptclose */
  /*     ================================================================== */
  cl__1.cerr = 0;
  cl__1.cunit = *iunit;
  cl__1.csta = 0;
  f_clos(&cl__1);

  return 0;
}

int snoptopen_(int *iunit)
{
  /* System generated locals */
  olist o__1;

  /* Builtin functions */
  int s_copy(char *, char *, int, int);
  int f_open(olist *);

  /* Local variables */
  static char lfile[20];

  /*     ================================================================= */
  /*     Open file named name to Fortran unit iunit.  inform .eq. 0        */
  /*     if successful. */
  /*     ================================================================= */
  s_copy(lfile, const_cast<char*>("testing.out"), (int)20, (int)11);
  o__1.oerr = 0;
  o__1.ounit = *iunit;
  o__1.ofnmlen = 20;
  o__1.ofnm = lfile;
  o__1.orl = 0;
  o__1.osta = const_cast<char*>("UNKNOWN_");
  o__1.oacc = 0;
  o__1.ofm = 0;
  o__1.oblnk = 0;
  f_open(&o__1);

  return 0;
}

#undef Complex

} // extern C
