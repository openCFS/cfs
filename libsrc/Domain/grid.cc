#include <vector>
#include <string>
#include <fstream>

#include "elements_header.hh"
#include "grid.hh"
#include "conffile.hh"

namespace CoupledField
{

Grid::Grid(FileType * aptFileType)
{
#ifdef TRACE
 (*trace) << " entering Grid::Grid \n";
#endif

  ptFileType=aptFileType;
  ptQ=new Quad1();
  ptTr=new Triangle1();
  ptTet=new Tetrahedral1();

  conf->getsubdom(listSD_);

}

Grid::~Grid()
{
 if (ptQ) delete ptQ;
 if (ptTr) delete ptTr;
 if (ptTet) delete ptTet;
}

}
