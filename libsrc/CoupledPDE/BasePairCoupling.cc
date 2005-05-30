#include "BasePairCoupling.hh"


#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/SinglePDE.hh"
#include "DataInOut/MaterialData.hh"
#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/LoadMaterialDataFile.hh"


namespace CoupledField
{
  
  BasePairCoupling::BasePairCoupling(SinglePDE *pde1, SinglePDE *pde2 ) {
    ENTER_FCN( "BasePairCoupling::BasePairCoupling" );
    
    pde1_ = pde1;
    pde2_ = pde2;

    ptGrid_ = pde1_->ptgrid_;
  }

  
  BasePairCoupling::~BasePairCoupling() {
    ENTER_FCN( "BasePairCoupling::~BasePairCoupling" );

    pde1_ = NULL;
    pde2_ = NULL;
               
  }

  void BasePairCoupling::Init(UInt bcSequenceStep,
                              std::string  bcSequenceTag) {
    ENTER_FCN( "BasePairCoupling::Init" );
    
    bcSequenceTag_ = bcSequenceTag;
    bcSequenceIndex_ = bcSequenceStep;
    
    
    // get subdomains of coupling object
    StdVector<std::string> keyVec, attrVec, valVec, regionNames, regionTypes;
    StdVector<RegionIdType> regionIds;
    
    // we are looking in coupling section of current tag
    attrVec = "tag", "", "", "";
    valVec = bcSequenceTag, "", "", ""; 
    
    // get coupling region names
    keyVec = "couplingList", "direct", couplingName_, "coupling", "name";
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptGrid_->RegionNameToId( regionIds, regionNames );

    // get coupling region types
    keyVec = "couplingList", "direct", couplingName_, "coupling", "type";
    params->GetList( keyVec, attrVec, valVec, regionTypes );
 
    // check if there are as many region types as region names
    if ( regionIds.GetSize() != regionTypes.GetSize() ) {
      (*error) << "BasePairCoupling::Init: There have to be as many region "
               << "types as region Names for the" << couplingName_ 
               << "-Coupling!\nPlease Check your input file!";
      Error( __FILE__, __LINE__ );
    }

    for ( UInt i=0; i < regionIds.GetSize(); i++ ) {
     
      if ( regionTypes[i] == "region" )
        subdoms_.Push_back( regionIds[i] );
      else if ( regionTypes[i] == "interface" ) 
        surfRegions_.Push_back( regionIds[i] );
      else {
        (*error) << "BasePairCoupling::Init: The region type '" 
                 << regionTypes[i] << "' is not known for direct coupling '"
                 << couplingName_ << "'! Please check your input file!";
        Error( __FILE__, __LINE__ );
      }
    }
//     std::cerr << "BasePairCoupling: My subdoms are:\n" << subdoms_ << std::endl;
//     std::cerr << "BasePairCoupling: My surfRegions are:\n" << surfRegions_ << std::endl;

    //std::cerr << "Name of pde1 = " << pde1_->GetName() 
    //<< std::endl;
    

    // Get type of analysis and create according 
    // assemble object
    // -> copy simply from first pde
    std::string help;
    Enum2String((*pde1_).analysistype_ , help);
    //std::cerr << "Analysis of PDE is " 
    //<< help << std::endl;
    switch ( (*pde1_).analysistype_ ) {

    case STATIC:
      assemble_ = new StaticAssemble(algsys_, ptGrid_);
      break;
    case TRANSIENT:
      assemble_ = new TransientAssemble(algsys_, ptGrid_);
      break;
    case HARMONIC:
      assemble_ = new HarmonicAssemble(algsys_, ptGrid_);
      break;
    default:
      Error (" analysistype was not found" , __FILE__, __LINE__ );
      break;
    }

    // Set general parameter of assemble class
    assemble_->SetGeneralParams(couplingName_, 1, subdoms_,
                                surfRegions_, bcSequenceTag );
                                

    // set PDE Ids to assemble object
    PdeIdType id1 = (*pde1_).pdeId_;
    PdeIdType id2 = (*pde2_).pdeId_;
    assemble_->SetPDEId( id1, id2 );
    
   

    // initialize nonlinearities
    
    // get the equation objects the assigned pdes
    NodeEQN * eqn1 = (*pde1_).eqnData_;
    NodeEQN * eqn2 = (*pde2_).eqnData_;
    assemble_->SetPtr2EQNData(eqn1, eqn2);

    // Read in material data
    ReadMaterialData();
    
    // Define the integrators
    DefineIntegrators();

    // Read the results which have to be stored
    
  }

  void BasePairCoupling::ReadMaterialData() {
    ENTER_FCN( "BasePairCoupling::ReadMaterialData" );

    // -------------------
    // XMLPARAMS
    // -------------------

    std::string matFileName;

    // Allocate space to hold material data for each subdomain of this PDE
    materialData_ = new MaterialData[subdoms_.GetSize()];
  
    // Get list of subdomains and materials
    StdVector< std::string > subdomName;
    StdVector< RegionIdType> subdomId;
    StdVector< std::string > subdomMaterial;
    params->GetList( "name", subdomName, "domain", "region" );
    params->GetList( "material", subdomMaterial, "domain", "region" );
    ptGrid_->RegionNameToId( subdomId, subdomName );
        
    // Query name of file with material data
    params->Get( "file", matFileName, "materialData" );
    
    // Generate new material reader
    LoadMaterialDataFile loadMaterialFile( matFileName.c_str() );
    
    // Load material data for subdomains on which this PDE lives
    // from data file
    for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
      for( UInt k = 0; k <= subdomName.GetSize(); k++ ) {
        if( subdoms_[i] == subdomId[k] ){
          loadMaterialFile.GetMaterial( materialData_[i], subdomMaterial[k],
                                        materialClass_ );
          break;
        }
      }
    }
  }
  
  // ======================================================
  // METHODS FOR ASSEMBLING
  // ======================================================
  
  void BasePairCoupling::SetupMatrixGraph() {
    assemble_->SetupMatrixGraph();
  }

  void BasePairCoupling::AssembleMatrices() {
    assemble_->AssembleMatrices();
  }
    
  void BasePairCoupling::AssembleSrcRHS(const Double time) {
    assemble_->AssembleSrcRHS(time);
  }
  
  void BasePairCoupling::AssembleNLRHS(const Double time) {
    assemble_->AssembleNLRHS(time);
  }

  void BasePairCoupling::AssembleSprings(const Double time) {
    assemble_->AssembleSprings(time);
  }




} // end of namespace
