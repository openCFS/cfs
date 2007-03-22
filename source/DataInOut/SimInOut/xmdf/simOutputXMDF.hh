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

  //! Class for reading in a mesh created by the ANSYS mkmesh-extension.

  //! Class, that is derived from class FileType for reading mesh-input data,
  //! which is produced by Ansys mkmesh-interface. 
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
    virtual void WriteData() {}
  
  protected:
    
    //! Transform type of elem in pointer to base class BaseFE

    //! This method maps the type number of an element - as given in the 
    //! mesh file - to a pointer to a reference finite element.
    //! \param itype (input) element type number as read in from the mesh

    FEType XMDFElemType2ElemType( const Integer type );
    Integer ElemType2XMDFElemType( const FEType type );

    void ReorderConnectivity( const Integer eType,
                              const bool toXMDF,
                              const UInt* in,
                              UInt* out);

    typedef std::vector< std::vector<UInt> > RegionElemType;
    typedef std::vector< std::set<UInt, std::less<UInt>, std::allocator<UInt> > > RegionNodeType;

    void WriteRegions(const H5::Group& meshGroup);
    void WriteNamedNodes(const H5::Group& meshGroup);
    void WriteNamedElems(const H5::Group& meshGroup);

    //@}

    // =======================================================================
    // CLASS ATTRIBUTES
    // =======================================================================
    //@{
    //! \name Attributes

    //@}

  private:
    void WriteAttribDescriptions();
    void WriteAttributes(const std::vector<std::string>& resultNames,
                         std::vector< Vector<Double>* >& results,
                         const std::string& regionStr,
                         const UInt numDOFs);
    void CreateExternalFile();

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
