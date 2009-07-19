// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_REFELEMGEN_2009
#define FILE_REFELEMGEN_2009

#include <vector>
#include <map>

#include "General/environment.hh"
#include "Domain/elem.hh"

namespace CoupledField
{
  class ReferenceElementGenerator 
  {
  public:
    ReferenceElementGenerator();
    virtual ~ReferenceElementGenerator() {};

    void GenGrid(std::vector<double>& coords,
                 std::vector<UInt>& connect,
                 std::vector<UInt>& elemTypes,
                 UInt& maxNumElemNodes,
                 std::map<std::string, std::vector<UInt> >& regionElems,
                 std::map<std::string, std::vector<UInt> >& nodeGroups,
                 Elem::FEType elemType);

  private:
  };
    
}

#endif
