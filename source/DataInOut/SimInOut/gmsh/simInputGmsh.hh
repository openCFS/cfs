// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GMSHREADER_2009
#define FILE_GMSHREADER_2009

#include <map>
#include <string>
#include <vector>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/simInput.hh"
#include "Domain/elem.hh"
#include "General/defs.hh"
#include "PDE/basePDE.hh"
#include "Utils/StdVector.hh"
// #include <boost/bimap/set_of.hpp>
#include "boost/bimap.hpp"

namespace CoupledField {

  class CoordSystem;
class Grid;
  
  /**
   **/
  class SimInputGmsh : public SimInput
  {
  private:

  public:
    SimInputGmsh(std::string fileName, PtrParamNode inputNode);
    virtual ~SimInputGmsh();

    virtual void InitModule();

    virtual void ReadMesh(Grid *mi);
    
    // =======================================================================
    // GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    virtual UInt GetDim();
    
    //! Get total number of nodes in mesh
    virtual UInt GetNumNodes();
    
    //! Get total number of elements in mesh
    virtual UInt GetNumElems( const Integer );
    
    //! Get total number of regions
    virtual UInt GetNumRegions();

    //! Get total number of named nodes
    virtual UInt GetNumNamedNodes();

    //! Get total number of named elements
    virtual UInt GetNumNamedElems();

    //@}

    virtual void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                           std::map<UInt, UInt>& numSteps,
                                           bool isHistory );
  
    // =======================================================================
    // ENTITY NAME ACCESS
    // =======================================================================
    //@{ \name Entity Name Access
  
    //! Get vector with all region names in mesh
    
    //! Returns a vector with the names of regions in the mesh of all
    //! dimensions.
    //! \param regionNames (output) vector containing names of regions
    //! \note Since the regionIdType is guaranteed to be defined by
    //! a number type (UInt, uint32), the regionId of an element can
    //! be directly used as index to the regions-vector
    virtual void GetAllRegionNames( StdVector<std::string> & regionNames );

    //! Get vector with region names of given dimension

    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements.
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    virtual void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                      const UInt dim );

    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    virtual void GetNodeNames( StdVector<std::string> & nodeNames );
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of named elements.
    //! \param elemNames (output) vector with names of named elements
    virtual void GetElemNames( StdVector<std::string> & elemNames );

  private:

    void ElemType2FEType(UInt eType, Elem::FEType& feType, UInt& numNodes);
    void InitElemNodeMap();
    void GetPhysEnt2RegionMapFromXML();
    void LinearizeElem(Elem::FEType* elemType);
    void TransformNodes(CoordSystem& coordSys, double scaleFac);

    struct NodeStruct 
    {
      UInt nodeId;
      Double x,y,z;
    };

    StdVector<Double> nodeCoords_;
    std::map<UInt, UInt> nodeNumMap_;
    std::vector< Elem::FEType > elementTypes_;
    std::vector< UInt > elementPhysicsTypes_;
    std::vector< UInt > connectivity_;
    UInt numNodesPerElem_;
    UInt dim_;
    std::map<UInt, std::string> physEntities2RegionNames_;
    std::map<std::string, bool> linearizeRegions_;
    bool readOnlySomeRegions_;
    std::string coordSysId_;
    Double scaleFac_;
  }; 

} 

#endif
