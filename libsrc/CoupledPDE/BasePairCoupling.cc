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
  }

  
  BasePairCoupling::~BasePairCoupling() {
    ENTER_FCN( "BasePairCoupling::~BasePairCoupling" );

    pde1_ = NULL;
    pde2_ = NULL;
	       
  }

  void BasePairCoupling::Init(Integer bcSequenceStep,
			      std::string  bcSequenceTag) {
    ENTER_FCN( "BasePairCoupling::Init" );
    
    bcSequenceTag_ = bcSequenceTag;
    bcSequenceIndex_ = bcSequenceStep;

    
    // get subdomains of coupling object
    StdVector<std::string> keyVec, attrVec, valVec;
    keyVec = "couplingList", "direct", couplingName_, "coupling", "name";
    attrVec = "tag", "", "", "";
    valVec = bcSequenceTag, "", "", "";
    params->GetList( keyVec, attrVec, valVec, subdoms_ );
 
    //std::cerr << "Name of pde1 = " << pde1_->GetName() 
    //<< std::endl;
    

    // Get type of analysis and create according 
    // assemble object
    // -> copy simply from first pde
    Grid *ptGrid = (*pde1_).ptgrid_;
    std::string help;
    Enum2String((*pde1_).analysistype_ , help);
    //std::cerr << "Analysis of PDE is " 
    //<< help << std::endl;
    switch ( (*pde1_).analysistype_ ) {

    case STATIC:
      assemble_ = new StaticAssemble(algsys_, ptGrid);
      break;
    case TRANSIENT:
      assemble_ = new TransientAssemble(algsys_, ptGrid);
      break;
    case HARMONIC:
      assemble_ = new HarmonicAssemble(algsys_, ptGrid);
      break;
    default:
      Error (" analysistype was not found" , __FILE__, __LINE__ );
      break;
    }

    // HARD CODED
    //
    //std::cerr << "couplingName = " << couplingName_ << std::endl;
    StdVector<std::string> surfdoms;
    assemble_->SetGeneralParams(couplingName_, 1, subdoms_,
				surfdoms, bcSequenceTag );
				

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
    StdVector< std::string > subdomMaterial;
    params->GetList( "name", subdomName, "domain", "region" );
    params->GetList( "material", subdomMaterial, "domain", "region" );
  
        
    // Query name of file with material data
    params->Get( "file", matFileName, "materialData" );
    
    // Generate new material reader
    LoadMaterialDataFile loadMaterialFile( matFileName.c_str() );
    
    // Load material data for subdomains on which this PDE lives
    // from data file
    for( Integer i = 0; i < subdoms_.GetSize(); i++ ) {
      for( Integer k = 0; k <= subdomName.GetSize(); k++ ) {
	if( subdoms_[i] == subdomName[k] ){
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

  void BasePairCoupling::AssembleMatrices(const Integer level) {
    assemble_->AssembleMatrices(level);
  }
    
  void BasePairCoupling::AssembleSrcRHS(const Integer level, 
					const Double time) {
    assemble_->AssembleSrcRHS(level, time);
  }
  
  void BasePairCoupling::AssembleNLRHS(const Integer level, 
				       const Double time) {
    assemble_->AssembleNLRHS(level, time);
 }

  void BasePairCoupling::AssembleSprings(const Integer level, 
					 const Double time) {
    assemble_->AssembleSprings(level, time);
 }




} // end of namespace
