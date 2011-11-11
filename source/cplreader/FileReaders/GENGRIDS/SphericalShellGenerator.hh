// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GENSPHERICALSHELL_2008
#define FILE_GENSPHERICALSHELL_2008

#include <vector>
#include <map>

#include "General/Environment.hh"

namespace CoupledField
{
  class SphericalShellGenerator 
  {
  public:
    SphericalShellGenerator();
    virtual ~SphericalShellGenerator() {};
    
    void GenSphericalShell(std::vector<double>& coords,
                           std::vector<UInt>& connect,
                           std::vector<UInt>& elemTypes,
                           UInt& maxNumElemNodes,
                           std::map<std::string, std::vector<UInt> >& regionElems,
                           std::map<std::string, std::vector<UInt> >& nodeGroups);
  private:
    int CalcHashVal(double x, double y, double z);
    
    void CreateBaseMesh();

    void RefineInnerQuads();

    void RefineQuad(int elemIdx,
                    std::vector<UInt>& connectOut,
                    std::vector<UInt>& regionsOut,
                    std::vector<UInt>& elemTypesOut,
                    bool generateElems = true);
    

    void TriangulateQuads(std::vector<UInt>& connectIn,
                          std::vector<UInt>& connectOut,
                          std::vector<UInt>& regionsIn,
                          std::vector<UInt>& regionsOut);
    

    void ProjectCoordsToSphere(std::vector<double>& coordsIn,
                               std::vector<double>& coordsOut,
                               double radius);
    

    void CreateOuterQuads();
    void CreateHexas();
    void CreatePyras();
    
    void MakeQuadraticQuads();
    
  private:
    UInt numLevels_;
    UInt numInnerPoints_;
    UInt numInnerElems_;
    Double innerRadius_;
    Double outerRadius_;
    UInt numLayers_;
    UInt maxNumElemNodes_;
    std::vector<std::string> regionNames_;
    std::vector<std::string> regionBaseNames_;
    bool octantRegions_;
    std::vector<Double> coords_;
    std::vector<UInt> connect_;
    std::vector<UInt> quadConnectFull_;
    std::vector<UInt> elemTypes_;
    std::vector<UInt> regions_;
    std::map<int, int> hashMap_;

    std::map<std::string, std::vector<UInt> > nodeGroups_;

    bool makeQuadratic_;
    bool makeSuperQuadratic_;
    
    
    
  };
    
}

#endif
