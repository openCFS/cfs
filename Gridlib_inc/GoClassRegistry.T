// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2000/07/21 10:40:49  prkipfer
| dropped g++ support
|
| Revision 1.2  2000/07/20 13:32:42  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.1  2000/07/03 16:15:16  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
|
+---------------------------------------------------------------------*/

//#include <ace/Synch_T.h>

template<class CREATOR>
GoClassRegistry<CREATOR>::GoClassRegistry()
//  : mtx_()
  : creators_()
{
}

template<class CREATOR>
void 
GoClassRegistry<CREATOR>::registerCreator(const GbString& className, CREATOR loc)
{
//  ACE_Guard<ACE_Thread_Mutex> guard(mtx_);

  debugmsg("registerCreator("<<className<<") called");
  
  creators_[className]=loc;
}

