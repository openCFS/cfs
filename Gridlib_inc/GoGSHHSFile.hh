/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOGSHHSFILE_HH
#define GOGSHHSFILE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iomanip>

#include "GoFile.hh"
#include "GoFileRegistry.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GoGSHHSFileRegistration;

class GoGSHHSFile
  : public GoFile
{
public:
  //! Constructor
  GoGSHHSFile(const GbString&);
  virtual ~GoGSHHSFile();

//  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;
  virtual void read(std::vector<GoMesh*>& mv) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoOffFile&);

private:
  //! hide copy constructor
  GoGSHHSFile(const GoGSHHSFile&);

  struct GSHHS {	/* Global Self-consistant Hierarchical High-resolution Shorelines */
    int id;				/* Unique polygon id number, starting at 0 */
    int n;				/* Number of points in this polygon */
    int level;			/* 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake */
    int west, east, south, north;	/* min/max extent in micro-degrees */
    int area;			/* Area of polygon in 1/10 km^2 */
    short int greenwich;		/* Greenwich is 1 if Greenwich is crossed */
    short int source;		/* 0 = CIA WDBII, 1 = WVS */
  };

  struct POINT {
    int	x;
    int	y;
  };

  static const GoGSHHSFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoGSHHSFileRegistration
{
public:
  GoGSHHSFileRegistration() {
//    debugmsg("registering "<<GoGSHHSFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoGSHHSFile::className_,
				     GoGSHHSFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("b"),
				       GoGSHHSFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoGSHHSFile::className_,
					 GbString("Shoreline binary file (*.b)"));
  }
};


//#ifndef OUTLINE
//#include "GoGSHHSFile.in"
//#endif

#endif // GOGSHHSFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/02/13 11:34:22  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.1  2001/01/02 14:26:04  prkipfer
| added GSHHS file format
|
|
+---------------------------------------------------------------------*/
