/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBIMAGEELEMENT_HH
#define  GBIMAGEELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <typeinfo>

#include "GbTypes.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

#define GBDECLAREELEMENT(type) \
class GbE##type \
{ \
public: \
    GbE##type (type tValue = 0) { value_ = tValue; } \
    \
    GbE##type& operator= (GbE##type kElement) \
    { \
        value_ = kElement.value_; \
        return *this; \
    } \
    \
    operator type () { return value_; } \
    \
    static int getRTTI () { return rttiId_; } \
    \
protected: \
    type value_; \
    static const int rttiId_; \
}

// wrappers for native types
const unsigned int GB_ELEMENT_QUANTITY = 12;

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short rgb5;
typedef unsigned int rgb8;

GBDECLAREELEMENT(char);
GBDECLAREELEMENT(uchar);
GBDECLAREELEMENT(short);
GBDECLAREELEMENT(ushort);
GBDECLAREELEMENT(int);
GBDECLAREELEMENT(uint);
GBDECLAREELEMENT(long);
GBDECLAREELEMENT(ulong);
GBDECLAREELEMENT(float);
GBDECLAREELEMENT(double);
GBDECLAREELEMENT(rgb5);
GBDECLAREELEMENT(rgb8);

#undef GBDECLAREELEMENT

#ifndef OUTLINE
#include "GbImageElement.in"
#endif

#endif  // GBIMAGEELEMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:25:26  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:15  prkipfer
| introduced image handling tool
|
|
+---------------------------------------------------------------------*/
