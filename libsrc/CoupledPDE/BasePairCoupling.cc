#include "BasePairCoupling.hh"


#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/PlainMaterialHandler.hh"
#include "Materials/piezoMaterial.hh"
#include "Driver/assemble.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  BasePairCoupling::BasePairCoupling(SinglePDE *pde1, SinglePDE *pde2 )    
  {

    ENTER_FCN( "BasePairCoupling::BasePairCoupling" );
    
    // initialize pointers
    sol_          = NULL;
    solVec_       = NULL;
    pde1_         = NULL;
    pde2_         = NULL;
    ptGrid_       = NULL;
    algsys_       = NULL;

    pde1_   = pde1;
    pde2_   = pde2;
    ptGrid_ = pde1_->ptgrid_;

    results1_ = pde1_->GetResults();
    results2_ = pde2_->GetResults();

    isaxi_ = false;
    
    
  }


  // **************
  //   Destructor
  // **************
  BasePairCoupling::~BasePairCoupling() {

    ENTER_FCN( "BasePairCoupling::~BasePairCoupling" );

    pde1_ = NULL;
    pde2_ = NULL;

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      delete it->second;
    }
    materials_.clear();

  }


  // ********
  //   Init
  // ********
  void BasePairCoupling::Init( UInt bcSequenceStep,
                               std::string bcSequenceTag ) {

    ENTER_FCN( "BasePairCoupling::Init" );
    
    bcSequenceTag_ = bcSequenceTag;
    bcSequenceIndex_ = bcSequenceStep;

    // Determine analysistype
    analysisType_ = pde1_->GetAnalysisType();
 
    
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

    // Determine, if axisymmetric geometry is used
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );
    if( probGeo == "axi" ) {
      isaxi_ = true;
    } else {
      isaxi_ = false;
    }

    // Get type of analysis and create according 
    // assemble object
    // -> copy simply from first pde
    std::string help;
    Enum2String((*pde1_).analysistype_ , help);
    //std::cerr << "Analysis of PDE is " 
    //<< help << std::endl;

    if ( algsys_ == NULL ) {
      (*error) << "BasePairCoupling::Init: The pointer to the algebraic "
               << "system was not set yet! You must call 'SetAlgSys()' "
               << "before!";
      Error( __FILE__, __LINE__ );
    }

 //    switch ( (*pde1_).analysistype_ ) {

//     case STATIC:
//       assemble_ = new StaticAssemble( algsys_, ptGrid_, bcSequenceTag_ );
//       break;
//     case TRANSIENT:
//       assemble_ = new TransientAssemble( algsys_, ptGrid_, bcSequenceTag_ );
//       break;
//     case HARMONIC:
//       assemble_ = new HarmonicAssemble( algsys_, ptGrid_, bcSequenceTag_ );
//       break;
//     case EIGENFREQUENCY:
//       assemble_ = new TransientAssemble( algsys_, ptGrid_, bcSequenceTag_ );
//       break;
// //     case MULTIHARMONIC:
// //       assemble_ = new MHassemble( algsys_, ptGrid_, bcSequenceTag_ );
// //       break;
//     default:
//       std::string myType;
//       Enum2String( (*pde1_).analysistype_, myType );
//       (*error) << "Unsupported analysis type '" << myType << "' detected "
//                << "(in PDE1)";
//       Error ( __FILE__, __LINE__ );
//       break;
//     }

 
    //assemble_ = pde1_->getPDE_assemble();

    eqnMap1_ = pde1_->GetEqnMap();
    eqnMap2_ = pde2_->GetEqnMap();
    assert( eqnMap1_ != NULL);
    assert( eqnMap2_ != NULL);

    // Read in material data
    ReadMaterialData();
 
    // Define the integrators
    DefineIntegrators();


    // define which solution types have to be saved
    ReadStoreResults();

  }

  void BasePairCoupling::ReadMaterialData() {
    ENTER_FCN( "BasePairCoupling::ReadMaterialData" );

    // -------------------
    // XMLPARAMS
    // -------------------

    std::string matFileName, actRegionName;
    StdVector< std::string > keyVec, attrVec, valVec;


    // Obtain pointer to materialHandler
    MaterialHandler * matLoader = NULL;
    matLoader = domain->GetMaterialHandler();
    
    // Get list of subdomains and materials
    StdVector< std::string > subdomName;
    StdVector< RegionIdType> subdomId;
    StdVector< std::string > subdomMaterial, subdomCoordSys;
    params->GetList( "name", subdomName, "domain", "region" );
    params->GetList( "material", subdomMaterial, "domain", "region" );
    params->GetList( "refCoordSys", subdomCoordSys, "domain", "region" );
    ptGrid_->RegionNameToId( subdomId, subdomName );

    // Load material data for subdomains on which this PDE lives
    // from data file
    for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
      for( UInt k = 0; k < subdomName.GetSize(); k++ ) {
        if( subdoms_[i] == subdomId[k] ){
          // Check if region contains a material
          if ( subdomMaterial[k] != "" ) {

            actRegionName = ptGrid_->RegionIdToName( subdomId[k] );
            materials_[subdoms_[i]] = matLoader->
              LoadMaterial( subdomMaterial[k],materialClass_ );
            
            // Check if coordinate system is present
            if( subdomCoordSys[k] != "" ) {
              CoordSystem * actCoosy = 
                domain->GetCoordSystem(subdomCoordSys[k]);
              materials_[subdoms_[i]]->SetCoordSys( actCoosy );
            }
            
            // Fetch for each material rotation paramters
            StdVector<Double> rotVecX, rotVecY, rotVecZ;
            Vector<Double> rotVec (3);
            rotVec.Init();

            bool isRotated = false;

            attrVec = "", "name", "";
            valVec = "", actRegionName, "";

            // xAxis
            keyVec = "domain", "region", "rotation", "xAxis";
            params->GetList( keyVec, attrVec, valVec, rotVecX );
            if( rotVecX.GetSize() == 1) {
              rotVec[0] = rotVecX[0];
              isRotated = true;
            }

            // yAxis
            keyVec = "domain", "region", "rotation", "yAxis";
            params->GetList( keyVec, attrVec, valVec, rotVecY );
            if( rotVecY.GetSize() == 1) {
              rotVec[1] = rotVecY[0];
              isRotated = true;
            }
            
            // zAxis
            keyVec = "domain", "region", "rotation", "zAxis";
            params->GetList( keyVec, attrVec, valVec, rotVecZ );
            if( rotVecZ.GetSize() == 1) {
              rotVec[2] = rotVecZ[0];
              isRotated = true;
            }
            
            if( isRotated ) {
              materials_[subdoms_[i]]->
                RotateAllTensorsByRotationAngles( rotVec, true );
            }
            break;
          }
        }
      }
    }

     // -------------------
    // COMPOSITE MATERIALS
    // -------------------
    Double startHeight;
    StdVector<std::string> compMaterials, subdomComposite;
    StdVector<Double> compOrient, compThick;
    
    params->GetList( "composite", subdomComposite, "domain", "region" );

    std::ostringstream out;
    for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
      for( UInt k = 0; k < subdomId.GetSize(); k++ ) {
        if( subdoms_[i] == subdomId[k] ){
          
          // Check if region contains a material
          if ( subdomComposite[k] != "" ) {
            
            actRegionName = ptGrid_->RegionIdToName( subdomId[k] );
            out << "Composite material '" << subdomComposite[k] << "' for "
                << "region '" << actRegionName << "' (ID = " << subdomId[k]
                << ") follows:\n";
            Info->PrintF( couplingName_, out.str().c_str());
            out.str("");
            
            // Read data from parameter file
            
            // StartHeight
            keyVec = "domain", "composite", "startHeight";
            attrVec = "", "name";
            valVec = "", subdomComposite[k];
            params->Get( keyVec, attrVec, valVec, startHeight );
            
            // materials of laminae
            attrVec = "", "name", "";
            valVec = "", subdomComposite[k], "";
            keyVec = "domain", "composite", "lamina", "material";
            params->GetList( keyVec, attrVec, valVec, compMaterials );
            
            // thickness of laminae
            keyVec = "domain", "composite", "lamina", "thickness";
            params->GetList( keyVec, attrVec, valVec, compThick );
            
            // orientation of laminae
            keyVec = "domain", "composite", "lamina", "orientation";
            params->GetList( keyVec, attrVec, valVec, compOrient );
            
            // Check, if all vectors have same length
            if ( compMaterials.GetSize() != compThick.GetSize() ||
                 compMaterials.GetSize() != compOrient.GetSize() ) {
              Error( "Inconsistent definition of composite material",
                     __FILE__, __LINE__ );
            }
            
            // Create new lamina and fill ine materials and thicknesses
            Composite & myMat = compositeMaterials_[subdomId[k]];
            UInt numLaminae = compMaterials.GetSize();
            myMat.name = subdomComposite[k];
            myMat.zStart = startHeight;
            myMat.thickness = compThick;
            myMat.orientation = compOrient;
            
            // Read single materials and print data
            myMat.materials.Resize( numLaminae );
            for ( UInt iLam = 0; iLam < numLaminae; iLam++ ) {
              
              // Print information
              out << " --- Lamina " << iLam+1 << ": thickness = " 
                  << compThick[iLam]
                  << " m, orientation = " << compOrient[iLam] << "° ---\n";
              Info->PrintF( couplingName_, out.str().c_str());
              out.str("");
              
              myMat.materials[iLam] = matLoader->
                LoadMaterial( compMaterials[iLam], materialClass_ );
            } // over  laminae
          }
        }
      }
    }
  }


  PdeIdType BasePairCoupling::GetPdeId1() {
    ENTER_FCN( "BasePairCoupling::GetPdeId1" );
    return pde1_->GetPDEId();
  }

  PdeIdType BasePairCoupling::GetPdeId2() {
    ENTER_FCN( "BasePairCoupling::GetPdeId2" );
    return pde2_->GetPDEId();
  }

} // end of namespace
