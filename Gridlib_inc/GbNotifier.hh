/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBNOTIFIER_H
#define GBNOTIFIER_H

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GbNotifier 
{ 
public:
  GbNotifier();
  virtual ~GbNotifier();

  virtual void name(const char* s);
  virtual void process(const char* s);
  virtual void message(const char* s);
  virtual void status(float f);
  virtual void completion(float f);

protected:
  GbString myname_;
  int complete_;
};

//  #ifndef OUTLINE
//  #include "GbNotifier.in"
//  #endif

#endif // GBNOTIFIER_H
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/04/17 12:16:36  prkipfer
| introduced user notifier
|
|
+---------------------------------------------------------------------*/
