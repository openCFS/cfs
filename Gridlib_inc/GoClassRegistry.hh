/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOCLASSREGISTRY_HH
#define GOCLASSREGISTRY_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

//#include <ace/OS.h>
//#include <ace/Synch.h>

#include <map>

#include "GbTypes.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

// a class registry template
// reused (!) from the Vision rendering framework
// 1997 Tom Peuker
// 1998 Peter Kipfer

template<class CREATOR>
class GoClassRegistry
{
public:
  virtual ~GoClassRegistry() {}

  void   registerCreator(const GbString& className, CREATOR);
  INLINE GbBool containsCreator(const GbString& className);
  INLINE CREATOR   queryCreator(const GbString& className);

protected:
  typedef std::map<GbString, CREATOR> ClassMap;

  ClassMap creators_;
//  ACE_Thread_Mutex mtx_;

  //! Constructor is protected to prevent direct instantiation
  GoClassRegistry();
};

#ifdef __GNUG__

#include "GoClassRegistry.in"
#include "GoClassRegistry.T"

#else

#ifndef OUTLINE
#include "GoClassRegistry.in"
#endif

#endif // g++

#endif // GOCLASSREGISTRY_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2000/07/21 10:40:49  prkipfer
| dropped g++ support
|
| Revision 1.2  2000/07/20 13:32:42  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.1  2000/07/03 16:15:17  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
