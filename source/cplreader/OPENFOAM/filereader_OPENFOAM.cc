
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream> 

#include "../params.hh"
#include "../settings.hh"
#include "filereader_OPENFOAM.hh"

namespace CoupledField
{

  FileReader_OPENFOAM::FileReader_OPENFOAM(const std::string& name,
                                           const UInt dim,
                                           const UInt numFiles) :
    FileReader(name, dim, numFiles)
  {
  }

  FileReader_OPENFOAM::~FileReader_OPENFOAM()
  {
  }

  void FileReader_OPENFOAM::Init()
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::Init" << std::endl;
  }


  void FileReader_OPENFOAM::ReadNodalCoords(std::vector<Double> & NODECOORD,
                                            const UInt partitionIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalCoords" << std::endl;

    // Hier werden die Koordinaten für eine Partition/Subdomain eingelesen.
    // partitionIdx kann hierbei von 0 ... n-1 gehen.
    // NODECOORD[0] = x1
    // NODECOORD[1] = y1
    // NODECOORD[2] = z1

    // NODECOORD[3] = x2
    // NODECOORD[4] = y2
    // NODECOORD[5] = z2

    // In ReadNodalValues müssen die Werte in der gleichen Reihenfolge
    // eingelesen werden.
  }

  void FileReader_OPENFOAM::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                         std::vector<UInt> & numNodesPerElem,
                                         std::vector<UInt> & elemTypes,
                                         const UInt partitionIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadTopology" << std::endl;

    // Hier werden die Elemente für eine Partition/Subdomain eingelesen.
    // partitionIdx kann hierbei von 0 ... n-1 gehen.
    // Für eine Übersicht über die Elementtypen siehe elementtypes.pdf

    // elemTypes[0] = ET_QUAD4
    // numNodesPerElem[0] = 4
    // TOPOLOGYDATA[0] = 1
    // TOPOLOGYDATA[1] = 2
    // TOPOLOGYDATA[2] = 3
    // TOPOLOGYDATA[3] = 4

    // elemTypes[1] = ET_TRIA3
    // numNodesPerElem[0] = 3
    // TOPOLOGYDATA[4] = 3
    // TOPOLOGYDATA[5] = 4
    // TOPOLOGYDATA[6] = 5
    //...
  }

  void FileReader_OPENFOAM::ReadNodalValues(std::vector<double> & flowdata,
                                            const UInt partitionIdx,
                                            const UInt timeStepIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalValues..." << std::endl;
    

    // flowdata layout:
    // flowdata[0] = Lighthill Quellterm 1
    // flowdata[1] = u1
    // flowdata[2] = v1 
    // flowdata[3] = w1
    // flowdata[4] = Druck1
    // flowdata[5] = N/A
    // flowdata[6] = N/A

    // flowdata[7] = Lighthill Quellterm 2
    // flowdata[8] = u2
    // flowdata[9] = v2
    // flowdata[10] = w2
    // flowdata[11] = Druck2
    // ...

  }

} // end of namespace
