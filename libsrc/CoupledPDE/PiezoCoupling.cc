#include "PiezoCoupling.hh"


#include "Driver/assemble.hh"
#include "DataInOut/MaterialData.hh"
#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

// integrator (bi-)linear forms
#include "Forms/linPiezoCoupling3D.hh"
#include "Forms/linPiezoCoupling2DAxi.hh"
#include "Forms/linPiezoCoupling2DPlaneStrain.hh"


namespace CoupledField
{

  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2 )
    : BasePairCoupling( pde1, pde2 ) {
    ENTER_FCN( "PiezoCoupling::PiezoCoupling" );

    couplingName_ = "piezoDirect";
    materialClass_ = "piezo";
  }


  PiezoCoupling::~PiezoCoupling() {
    ENTER_FCN( "PiezoCoupling::~PiezoCoupling" )
      }

  void PiezoCoupling::PostProcess() {
    ENTER_FCN( "PiezoCoupling::PostProcess" );
  }
  

  void PiezoCoupling::DefineIntegrators() {
    ENTER_FCN( "PiezoCoupling::DefineIntegrators" );
    

    piezoMaterialType realMatParameter = REALMATERIALPARAMETER;
    
    
    // iterate over all subdomains
    for ( UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {
       
      // ==============  add stiffness ========================================
      BaseForm * bilinearStiff = GetStiffIntegrator(&materialData_[actSD]);
      IntegratorDescriptor *actIntDescrStiff = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      bilinearStiff->SetPiezoMaterialType(realMatParameter);
      actIntDescrStiff->SetPiezoMaterialType(realMatParameter);
      actIntDescrStiff->SetPDEIds(pde1_, pde2_);      
      assemble_->AddIntegrator(actIntDescrStiff, subdoms_[actSD]);     
    }

  }
  BaseForm * PiezoCoupling::GetStiffIntegrator(MaterialData * actSDMat,
                                               Boolean reducedInt,
                                               Boolean isdamping) {

    ENTER_FCN( "PiezoCoupling::GetStiffIntegrator" );

    // Get problem geometry and mechanic subtype
    std::string probGeo, subType;
    params->Get( "subType", subType, "mechanic" );
    params->Get( "type", probGeo, "geometry" );

    BaseForm * bilinearStiff;

    if (subType == "planeStrain") {
      bilinearStiff = new linPiezoCoupling2DPlaneStrain();
    }
    else if (subType == "axi") { 
      bilinearStiff = new linPiezoCoupling2DAxi();
    }
    else if (subType == "3d") {
      bilinearStiff = new linPiezoCoupling3D(); 
    }
    else
      Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

    // Set pointer to material type
    bilinearStiff->SetMaterial( actSDMat );
    
    return bilinearStiff;
  }
  
  void PiezoCoupling::ReadStoreResults() {
    ENTER_FCN( "PiezoCoupling::ReadStoreResults" );
    
  }


} // end of namespace
