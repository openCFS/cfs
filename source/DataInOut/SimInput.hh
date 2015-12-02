// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMINPUT_2006
#define FILE_SIMINPUT_2006

#include <string>
#include <vector>

#include "General/Environment.hh"
#include "Utils/tools.hh"
#include "General/Exception.hh"

#include "Domain/Results/BaseResults.hh"
#include "Domain/Mesh/Grid.hh"

#include "PDE/BasePDE.hh"

namespace CoupledField
{

  //! Forward class declaration
  class ParamNode;
  class BaseResult;

  //! Abstract base class for handling exceptions and errors
  class ErrorHandler {

  public:
    //! Constructor
    ErrorHandler() {};

    //! Destructor
    virtual ~ErrorHandler() {};

    //! Error (non-recoverable)
    virtual void Error( const Exception &exc ) = 0;

    //! WARN (recoverable)
    virtual void Warning( const Exception &exc ) = 0;
  };

  //! Abstract class for reading in mesh data

  //! This class defines an abstract interface for accessing 
  //! files containing geometric mesh information. 
  //!
  //! \note All mesh and geometric entities are counted one-based,
  //! whereas the acces to the C++ built in datatypes is zero-based!
  
  class SimInput
  {

  public:

    //! Define capabilities of specific input format
    typedef enum {NONE, MESH, MESH_RESULTS, HISTORY} Capability;

    // ========================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // ========================================================================
    //@{ \name Constructor / Initialization
  

    //! Constructor with name of mesh-file
      SimInput(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode ) :
          fileName_(fileName),
          mi_(NULL),
          myParam_(inputNode),
          myInfo_(infoNode),
          dim_(0)
      {};


    //! Destructor
    virtual ~SimInput() {};

    //! Get capabilites of interface
    const std::set<Capability>& GetCapabilities() const 
    { return capabilities_; }
    //@}

    virtual void InitModule() = 0;
      
    virtual void ReadMesh( Grid *mi ) = 0;

    //! Tries to copy all relevant mesh related information rather than creating a
    //! mesh of its own. USefull if different results are defined in different
    //! result files and one wants to combine those. Works only!!! for identical meshes.
    virtual void CopyMeshInfo( shared_ptr<SimInput> otherInput ){
      EXCEPTION("Copy of mesh data from other simInputs is not supported");
    }

    // ========================================================================
    // GENERAL MESH INFORMATION
    // ========================================================================
    //@{ \name General Mesh Information

    //! Get dimension of the mesh
    virtual UInt GetDim() = 0;

    /** @return includes the complete path */
    const std::string& GetFileName() const { return fileName_; }

    //! Get total number of nodes in mesh
    virtual UInt GetNumNodes() = 0;
 
    //! Get total number of elements in mesh
    virtual UInt GetNumElems( const Integer dim = -1 ) = 0;

    //! Get total number of regions
    virtual UInt GetNumRegions() = 0;

    //! Get total number of named nodes
    virtual UInt GetNumNamedNodes() = 0;

    //! Get total number of named elements
    virtual UInt GetNumNamedElems() = 0;

    //@}
  
    // =========================================================================
    // ENTITY NAME ACCESS
    // =========================================================================
    //@{ \name Entity Name Access
  
    //! Get vector with all region names in mesh
 
    //! Returns a vector with the names of regions in the mesh of all
    //! dimensions.
    //! \param regionNames (output) vector containing names of regions
    //! \note Since the RegionIdType is guaranteed to be defined by
    //! a number type (UInt, uint32), the regionId of an element can
    //! be directly used as index to the regionNames-vector
    virtual void GetAllRegionNames( StdVector<std::string> & regionNames ) = 0;

    //! Get vector with region names of given dimension

    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements.
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    virtual void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                      const UInt dim ) = 0;

    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    virtual void GetNodeNames( StdVector<std::string> & nodeNames ) = 0;
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of named elements.
    //! \param elemNames (output) vector with names of named elements
    virtual void GetElemNames( StdVector<std::string> & elemNames ) = 0;
    //@}

    // =========================================================================
    // GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information

    //! Return multisequence steps and their analysistypes
    virtual void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                           std::map<UInt, UInt>& numSteps,
                                           bool isHistory = false ) {
      analysis.clear();
      numSteps.clear();
    }

    //! Obtain list with result types in each sequence step
    virtual void GetResultTypes( UInt sequenceStep, 
                                 StdVector<shared_ptr<ResultInfo> >& infos,
                                 bool isHistory = false ) {
      infos.Clear();
    }
    
    //! Return list with time / frequency values and step for a given result
    virtual void GetStepValues( UInt sequenceStep,
                                shared_ptr<ResultInfo> info,
                                std::map<UInt, Double>& steps,
                                bool isHistory = false ) {
      steps.clear();
    }
    
    //! Return entitylist the result is defined on
    virtual void GetResultEntities( UInt sequenceStep,
                                    shared_ptr<ResultInfo> info,
                                    StdVector<shared_ptr<EntityList> >& list,
                                    bool isHistory = false) {
      list.Clear();
    }
    
    //! Fill pre-initialized results object with values of specified step
    virtual void GetResult( UInt sequenceStep,
                            UInt stepValue,
                            shared_ptr<BaseResult> result,
                            bool isHistory = false ) {
      EXCEPTION( "Not implemented in base class" );
    } 
    //@}

    shared_ptr<BaseResult> 
    GetResult( UInt sequenceStep,
               UInt stepValue,
               SolutionType solType,
               const std::string& regionName );
    
    template<typename TYPE>
    shared_ptr<FeFunction<TYPE> >
    GetFeFunction( UInt sequenceStep,
                   UInt stepValue,
                   SolutionType solType,
                   std::set<std::string> & regionNames );
    
    
  protected:

    //! Name of input file
    std::string fileName_;


    Grid *mi_;

    //! Capabilities of output class
    std::set<Capability> capabilities_;
    
    //! Parameter node for current output class
    PtrParamNode myParam_;
    
    //! Info node for current output class
    PtrParamNode myInfo_;

    UInt dim_;
    std::vector<UInt> numElemsOfDim_;
    std::map<UInt, StdVector<std::string> > regionNamesOfDim_;
    StdVector<std::string> namedElems_;
    StdVector<std::string> namedNodes_;
    
  };

}

#endif // FILE_FILETYPE

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
