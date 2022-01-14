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
    

    //! Return dimension of the mesh
    virtual UInt GetDim();
    

    // =========================================================================
    //  GENERAL SOLUTION INFORMATION
    // =========================================================================

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
