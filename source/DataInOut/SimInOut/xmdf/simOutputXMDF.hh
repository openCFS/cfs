// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMOUTPUTXMDF_2006
#define FILE_SIMOUTPUTXMDF_2006

#include <set>
#include <map>

#include <Domain/grid.hh>
#include <Domain/resultInfo.hh>
#include <DataInOut/simOutput.hh>

#include "H5Cpp.h"
#include "Xmdf.h"

namespace CoupledField {

  class SimOutputXMDF: virtual public SimOutput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimOutputXMDF(std::string fileName, ParamNode * inputNode);
    
    //! Destructor
    virtual ~SimOutputXMDF();

    //@}


    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step, AnalysisType type );
    
    //! Register result (within one multisequence step)
    virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                 UInt saveBegin, UInt saveInc,
                                 UInt saveEnd );

    //! Begin single analysis step
    virtual void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    virtual void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    virtual void FinishStep( );

    //! End multisequence step
    virtual void FinishMultiSequenceStep( );

    //! Finalize the output
    virtual void Finalize();


    virtual void InitModule();

    virtual void WriteGrid();
  
  private:
    FEType XMDFElemType2ElemType( const Integer type );
    Integer ElemType2XMDFElemType( const FEType type );

    void ReorderConnectivity( const Integer eType,
                              const bool toXMDF,
                              const UInt* in,
                              UInt* out);

    void WriteRegions(const H5::Group& meshGroup);
    void WriteNamedNodes(const H5::Group& meshGroup);
    void WriteNamedElems(const H5::Group& meshGroup);

    void WriteAttribDescriptions();
    void WriteAttributes(const std::vector<std::string>& resultNames,
                         std::vector< Vector<Double>* >& results,
                         const std::string& regionStr,
                         const UInt numDOFs);
    void CreateExternalFile();

  private:

    bool gridWritten_;

    bool externalFiles_;

    H5::Group mainRoot_;
    H5::Group dataGroup_;
    H5::Group volDataGroup_;
    H5::Group currMSGroup_;
    H5::Group currAttrDescGroup_;
    H5::Group currStepGroup_;
    H5::H5File currStepFile_;

    UInt currMS_;
    UInt currStep_;

    AnalysisType currAnalysisType_;

    typedef std::map< std::string, std::vector< 
                                     boost::shared_ptr<BaseResult>
                                     >
                    > ResDescType;
    ResDescType registeredResults_;

    xid fileId_;
    UInt multiStep_, step_;
    bool msChange_;

    typedef struct region_desc_type { 
      char name[32]; 
      UInt dim; 
    };

    typedef struct named_entity_desc_type { 
      char name[32]; 
    };
    
  };

}

#endif // FILE_MOD_XMDF_2006
