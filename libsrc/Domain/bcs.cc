#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "bcs.hh"

namespace CoupledField
{

Boolean NodeRestraint::operator < (const NodeRestraint & t)
{
 if (nodalnum < t.nodalnum) return 1;
 else 0;
}

BCs :: BCs(FileType * const aInFile)
{
#ifdef TRACE
  (*trace) << "entering BCs::BCs" << std::endl;
#endif
 InFile_     = aInFile; 
}

BCs :: ~BCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::~BCs" << std::endl;
#endif
;
}

void BCs :: ReadBCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::ReadBCs" << std::endl;
#endif

 InFile_->ReadBoundRestr(restr_[0], numrestr_[0]);

}

}
