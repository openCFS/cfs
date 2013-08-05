// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include "def_cplreader.hh"
#include "cplreader/FileReader.hh"
#include "cplreader/FlowDataTypes.hh"

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "integlib/elemIntegr.hh"

#include "FlowDataTypes.hh"
#include "OutputWriter.hh"
#include "CouplingHandler.hh"

#include "Meshes.hh"
#include "FlowDataLoader.hh"



//! Base class for reading topology and nodal data from results of fluid computations
/*!
  This class contains methods for the reading information
  from a fluid topology file and from the fluid results files.
  It reads and couples flow results from a set of files (or set of
  region files coming from a blocked structured fluid computation).
*/

using namespace CCM;
namespace CoupledField
{

  class FileReader_CCM : public FileReader
  {
  public:

    //! Constructor
	  FileReader_CCM(const std::string& name, const UInt dim, \
        const UInt numFiles, const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_CCM();

    virtual void Init();


    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD);


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                 std::vector<UInt> & elemTypes);

    //! get nodal values from the corresponding fluid datafile the new way
    virtual void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                 const std::vector<bool>& activeParts,
                                 const UInt timeStepIdx);

    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx);
    
  protected:
    CPLMesh cplmesh_;
    
    uint firstStep_;
    uint stepOffset_;
    uint digits_;
    
    FlowDataLoader* loader;
    
    std::string velocityName;
    std::string laplacePName;
    std::string pressureName;
    std::string divLHTName;
  };


} // end of namespace
