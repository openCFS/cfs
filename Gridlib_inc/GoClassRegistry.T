// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
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

