/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOGRCFILE_HH
#define GOGRCFILE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iomanip>

#include "GoFile.hh"
#include "GrCamera.hh"

#ifdef IRIX
// I know that read/write signatures don't match...
#pragma set woff 1681
#endif

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  This class implements the gridlib renderer config file format. This
  does not contain geometry information rather than a reference to the
  mesh file. Therefore, we do not register this format. You have to construct
  this class explicitly.
*/

class GoGrcFile
  : public GoFile
{
public:
  //! Constructor
  GoGrcFile(const GbString&);
  virtual ~GoGrcFile();

  virtual void write(const GrCamera& cam) const;
  virtual void read(GrCamera& cam) const;

  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoGrcFile&);

private:
  //! hide copy constructor
  GoGrcFile(const GoGrcFile&);
};

//#ifndef OUTLINE
//#include "GoGrcFile.in"
//#endif

#endif // GOGRCFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/19 16:30:54  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.1  2001/06/18 11:00:14  prkipfer
| added gridlib renderer config file format
|
|
+---------------------------------------------------------------------*/
