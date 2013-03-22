
#if defined (__cplusplus)
extern "C" {
#endif

extern long long ___chkstk ();

long long
__wrap___chkstk ()
{
//  printf ("chkstk called\n");
  return ___chkstk ();
}

#if defined (__cplusplus)
}
#endif
