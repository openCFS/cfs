// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_RECTLINUNIFORMGEN_2008
#define FILE_RECTLINUNIFORMGEN_2008

#include <vector>
#include <map>

#include "General/Environment.hh"

namespace CoupledField
{
  class RectilinearUniformGridGenerator 
  {
  public:
    RectilinearUniformGridGenerator();
    virtual ~RectilinearUniformGridGenerator() {};

    enum HexaType 
    {
      LINEAR,
      SERENDIPITY,
      LAGRANGE
    };
    
    void GenUniformGrid(std::vector<double>& coords,
                        std::vector<UInt>& connect,
                        std::vector<UInt>& elemTypes,
                        UInt& maxNumElemNodes,
                        std::map<std::string, std::vector<UInt> >& regionElems,
                        std::map<std::string, std::vector<UInt> >& nodeGroups,
                        HexaType hexaType);

  private:
    void CreateMesh(const std::vector<Double>& xCoords,
                    const std::vector<Double>& yCoords,
                    const std::vector<Double>& zCoords,
                    const std::vector<UInt>& xNElems,
                    const std::vector<UInt>& yNElems,
                    const std::vector<UInt>& zNElems,
                    HexaType hexaType);
  private:
    UInt maxNumElemNodes_;
    std::vector<std::string> regionNames_;
    std::vector<std::string> regionBaseNames_;
    std::vector<Double> coords_;
    std::vector<UInt> connect_;
    std::vector<UInt> elemTypes_;
    std::vector<UInt> regions_;

    std::map<std::string, std::vector<UInt> > nodeGroups_;
  };
    
}

#endif
