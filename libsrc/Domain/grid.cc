#include "elements_header.hh"
#include "grid.hh"

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
}

Grid::~Grid()
{
 if (ptQ) delete ptQ;
 if (ptTr) delete ptTr;
 if (ptTet) delete ptTet;
}

}
