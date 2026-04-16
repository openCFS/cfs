// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMINPUT_2006
#define FILE_SIMINPUT_2006

#include <string>
#include <vector>
#include <map>

#include "General/Exception.hh"

#include "Domain/Mesh/Grid.hh"

#include "PDE/BasePDE.hh"

namespace CoupledField
{

  //! Forward class declaration
  class BaseResult;
  struct ResultInfo;
  template<class TYPE> class FeFunction;

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

    /** initialize the SimInputModule(), late constructor */
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
                            bool isHistory = false) {
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
    
    
    void SetTempRegionName(const std::string & name){
      tempRegionName_ = name;
    }

    void ResetTempRegionName(){
      tempRegionName_ = "";
    }

  protected:

    //! Name of input file
    std::string fileName_;

    Grid *mi_;

    //! Capabilities of input class
    std::set<Capability> capabilities_;
    
    //! Parameter node for current input class
    PtrParamNode myParam_;
    
    //! Info node for current input class
    PtrParamNode myInfo_;

    UInt dim_;
    std::vector<UInt> numElemsOfDim_;
    std::map<UInt, StdVector<std::string> > regionNamesOfDim_;
    StdVector<std::string> namedElems_;
    StdVector<std::string> namedNodes_;

    std::string tempRegionName_;

  private:  
    
    //! Class for identifying cached information in GetResult() function.
    class GetResultArguments {
    
    public:
      //! Contructor taking the definite arguments of GetResult() function
      GetResultArguments(UInt sequenceStep, SolutionType solType, std::string regionName);
      
      //! Less comparator for usage in std::map
      bool operator<(const GetResultArguments &r) const;
      
    private:
      UInt sequenceStep_;
      SolutionType solType_;
      std::string regionName_;
  
    };
    
    //! Maps the arguments of the GetResult() function against raw BaseResult objects to be used
    std::map< GetResultArguments, shared_ptr<BaseResult> > rawBaseResults_;
    
  };

}

#endif // FILE_FILETYPE

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
