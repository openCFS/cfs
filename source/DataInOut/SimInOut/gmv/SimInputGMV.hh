// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GMVREADER_2006
#define FILE_GMVREADER_2006

// This code is also based on the gmvread library from the official GMV
// website at: http://www-xdiv.lanl.gov/XCM/gmv/GMVHome.html

#include <map>

#include <DataInOut/SimInput.hh>

namespace CoupledField {

  /**
   **/
  class SimInputGMV : public SimInput
  {
  private:

    struct GridAttributeInfo
    {
    public:

      GridAttributeInfo() : dim_(0), elemAttrib_(false) 
      {
      }

      UInt dim_;
      bool elemAttrib_;
      std::vector<Double> data_;
    };

    Grid *mi_;
    std::vector< UInt > mElementMaterials;
    std::vector< UInt > mElementIds;
    std::vector< UInt > mVertexIds;
    StdVector< std::string > mMatNames;
    StdVector< RegionIdType > mMatNums;
    std::map< std::string, GridAttributeInfo > mGridAttributes;
    std::vector< UInt > mCycleNos;
    std::vector< Double > mProbTimes;
    static std::vector< std::string > mPossibleAttribs;

  public:
      SimInputGMV(std::string fileName, PtrParamNode inputNode, 
                  PtrParamNode infoNode );
      
    virtual ~SimInputGMV();

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

    // =========================================================================
    //  GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information

    //! Return multisequence steps and their analysistypes
    void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                   std::map<UInt, UInt>& numSteps,
                                   bool isHistory = false );

    //! Obtain list with result types in each sequence step
    void GetResultTypes( UInt sequenceStep, 
                         StdVector<shared_ptr<ResultInfo> >& infos,
                         bool isHistory = false );
    
    //! Return list with time / frequency values and step for a given result
    virtual void GetStepValues( UInt sequenceStep,
                                    shared_ptr<ResultInfo> info,
                                    std::map<UInt, Double>& steps,
                                    bool isHistory = false );
    
    //! Return entitylist the result is defined on
    void GetResultEntities( UInt sequenceStep,
                               shared_ptr<ResultInfo> info,
                               StdVector<shared_ptr<EntityList> >& list,
                               bool isHistory = false );

    //! Fill pre-initialized results object with values of specified step
    void GetResult( UInt sequenceStep,
                      UInt stepNum,
                      shared_ptr<BaseResult> result,
                      bool isHistory = false );
    //@}
  
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

    bool ProcessMesh();
    bool ProcessMaterials();
    bool ProcessGroups();
    bool ProcessVariables();
    bool ProcessVelocities();
    bool ProcessVectors();
    bool ProcessProbtime();
    bool ProcessCycleNo();
    bool SetupGridAndRegions();

    void Cleanup();
  }; 

} 

#endif
