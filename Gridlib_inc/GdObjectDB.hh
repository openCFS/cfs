/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GDOBJECTDB_HH
#define GDOBJECTDB_HH

#ifndef INLINE
#ifdef OUTLINE
#define INLINE
#else
#define INLINE inline
#endif
#endif

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algo.h>

using namespace std;

#include "GdDumpable.hh"
#include "GdDumpablePtr.hh"

/*----------------------------------------------------------------------
|       definitions
+---------------------------------------------------------------------*/

#if defined(DEBUG)
#define REGISTER_OBJECT(CLASS) \
	GdObjectDB::instance()->registerObject(new GdDumpableAdapter<CLASS>(this));
#define REMOVE_OBJECT \
	GdObjectDB::instance()->removeObject((void *)this);
#else
#define REGISTER_OBJECT(CLASS)
#define REMOVE_OBJECT
#endif

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GdObjectDB
{
public:
  //! Iterate through the entire set of registered objects
  //! and dump their state
  void dumpObjects();

  //! Add the tuple <dumper,this_> to the list of registered objects
  void registerObject(const GdDumpable *dumper);

  //! Use this_ to locate and remove the associated dumper from the list
  //! of registered objects
  void removeObject(const void *this_ptr);

  //! Factory method to get the singleton database
  static GdObjectDB *GdObjectDB::instance();

private:
  //! Singleton instance of this class
  static GdObjectDB *instance_;

  //! Ensure we have a singleton (nobody can create instances but
  //! this class)
  GdObjectDB();

  struct Tuple {
//    Tuple(const void *t, const GdDumpablePtr p)
//      : this_(t) , dumper_(p) {}

    //! Pointer to the registered C++ object
    const void *this_;
    
    //! Smart pointer to the GdDumpable object associated with this_
    const GdDumpablePtr dumper_;
  };

  typedef vector<Tuple> TupleVector;

  //! Holds all registered C++ objects
  TupleVector objectTable_;

  struct DumpTuple {
    bool operator()(const Tuple &t) {
      t.dumper_->dump();
    }
  };

  // STL predicate function object
  struct ThisMemberEqual
    : public binary_function<Tuple,Tuple,bool>
  {
    bool operator()(const Tuple &t1, const Tuple &t2) const {
      return t1.this_ == t2.this_;
    }
  };

};

//#ifndef OUTLINE
//#include "GdObjectDB.in"
//#endif

#endif // GDOBJECTDB_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1.1.1  2000/06/08 16:24:43  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
