/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOFILEREGISTRY_HH
#define GOFILEREGISTRY_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

//#include <ace/Singleton.h>
#include <map>

#include "GbTypes.hh"
#include "GoClassRegistry.hh"

class GoFile;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

//! prototype of class constructor
typedef GoFile* (*GoFileCreator)(const GbString& filename);


class GoFileRegistry 
  : public GoClassRegistry<GoFileCreator>
{
public:
  GoFileRegistry();
  virtual ~GoFileRegistry() {}
  
  void   registerExtension(const GbString& ex, const GbString& className);
  INLINE GbBool containsExtension(const GbString& ex);
  INLINE GbString  queryExtension(const GbString& ex);
  INLINE GbString  listExtensions();
  
  void   registerDescription(const GbString& className, const GbString& dsc);
  INLINE GbString  queryDescription(const GbString& className);
  INLINE GbString  listDescriptions();

  INLINE int numEntries() const;
  
protected:
  typedef std::map<GbString, GbString> NameMap;

  NameMap names_;
  NameMap descriptions_;
};

#ifndef __GNUG__
#pragma instantiate GoClassRegistry<GoFileCreator>
#endif

//typedef ACE_Singleton<GoFileRegistry, ACE_Recursive_Thread_Mutex> GO_FILE_REGISTRY;

#ifndef OUTLINE
#include "GoFileRegistry.in"
#endif

#endif // GOFILEREGISTRY_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.4  2000/07/25 09:10:39  prkipfer
| new list function
|
| Revision 1.3  2000/07/21 10:40:50  prkipfer
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
