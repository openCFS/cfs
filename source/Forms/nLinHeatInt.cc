#include <stddef.h>
#include <string>

#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "Forms/linHeatCondInt.hh"
#include "Forms/linearForm.hh"
#include "General/environment.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/ApproxData.hh"
#include "Utils/nodestoresol.hh"
#include "laplaceInt.hh"
#include "nLinHeatInt.hh"

namespace CoupledField {



  // ==============================================
  //   Constructors: heat stiffness bilinear form
  // ==============================================
	
	
  nlinHeatStiffInt::nlinHeatStiffInt(BaseMaterial* matData, SubTensorType type,
                          bool coordUpdate ) 
    : linHeatCondInt(matData, type, coordUpdate ) {
    
    
    name_ = "nlinHeatStiffInt";
    isSolDependent_ = true;
  }




  void nlinHeatStiffInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 ) {
  
   
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrFncs); 
    elemMat.Init();

    // get solution of current element
    Matrix<Double> temp;
    sol_->GetElemSolutionAsMatrix( temp, ent1 );
    temp.ConvertToVec_AppendRows( temperature_ );

    //nonlineasr function for approximation
    nlinFnc_ = ptMaterial->GetNonlinFnc(HEAT_CONDUCTIVITY );
    
    Double heatConductivity = 1.0;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)  {
      //get shape functions at compute temperature at integration point
      ptelem->GetShFncAtIp(ShpFncAtIp, actIntPt, ent1.GetElem() );
      Double tempAtIP = temperature_ * ShpFncAtIp;

      //get nonlinear conductivity value
      heatConductivity = nlinFnc_->EvaluateFunc( tempAtIP);
      //      std::cout << "T: " << tempAtIP << ";  lambda=" << heatConductivity << std::endl;

      //get derivatives
      jacDet = 0;      
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                    jacDet, ent1.GetElem());

      xiDx.Transpose(xiDxTransp);
      partElemMat = xiDx * xiDxTransp;
        
      if (isaxi_)  {
        CoordAtIP = ptCoord_ * ShpFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] 
          * jacDet * heatConductivity * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * jacDet * heatConductivity;

      elemMat += partElemMat;
    }
    //    std::cout << "ke:\n" << elemMat << std::endl;

  }


  // ==============================================
  //   Constructors: heat mass bilinear form
  // ==============================================
	
	
  nlinHeatMassInt::nlinHeatMassInt(BaseMaterial* matData, SubTensorType type,
                                   bool coordUpdate ) 
    : linHeatCondInt(matData, type, coordUpdate ) {
    
    
    name_ = "nlinHeatMassInt";
    isSolDependent_ = true;
    matData->GetScalar(density_,DENSITY,Global::REAL);
  }




  void nlinHeatMassInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 ) {
  
   
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrFncs); 
    elemMat.Init();

    // get solution of current element
    Matrix<Double> temp;
    sol_->GetElemSolutionAsMatrix( temp, ent1 );
    temp.ConvertToVec_AppendRows( temperature_ );

    //nonlineasr function for approximation
    nlinFnc_ = ptMaterial->GetNonlinFnc(HEAT_CAPACITY );
    
    Double heatCapacity = 1.0;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
      
      ptelem->GetShFncAtIp(ShpFncAtIp, actIntPt, ent1.GetElem() );

      //tempertaure at integration point
      Double tempAtIP = temperature_ * ShpFncAtIp;
      
      //get nonlinear capacity value
      heatCapacity = nlinFnc_->EvaluateFunc( tempAtIP);

      partElemMat.DyadicMult(ShpFncAtIp);
      
      if (isaxi_) {
        CoordAtIP = ptCoord_ * ShpFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * heatCapacity * jacDet * CoordAtIP[0] * density_;
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * heatCapacity * jacDet * density_;
      
      elemMat += partElemMat;
    }

    //  std::cout << "ke:\n" << elemMat << std::endl;

  }


  //============== NL-RHS ===================================================================
  nLinHeat_linFormInt::nLinHeat_linFormInt( BaseMaterial* matData, 
                                            bool axi,
                                            bool coordUpdate)
    : LinearForm( matData ) {

    name_ = "nLinHeat_linFormInt";
    isSolDependent_ = true;

    isaxi_       = axi;
    coordUpdate_ = coordUpdate;

    // need nonlinear curve approximation for conductivity and capacity
  }
  
  nLinHeat_linFormInt::~nLinHeat_linFormInt()
  {
  }

  void nLinHeat_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );

    // get pointer to nonlinear curve approximation
    nlinFnc_ = ptMaterial->GetNonlinFnc(HEAT_CONDUCTIVITY);

    BaseForm * heatForm;
    if ( nlinFnc_ == NULL ) {
      //define the linear element matrix
      Double coeffstiff;
      ptMaterial->GetScalar(coeffstiff,HEAT_CONDUCTIVITY,Global::REAL);
      heatForm = new LaplaceInt(coeffstiff, isaxi_, coordUpdate_ ); 
    }
    else {
      //define the nonlinear element matrix
      SubTensorType tensorType = FULL;
      if ( ptelem->GetDim() == 2 ) {
        if ( isaxi_ )
          tensorType = AXI;
        else
          tensorType = PLANE_STRAIN;
      }
      heatForm = new nlinHeatStiffInt( ptMaterial, tensorType, coordUpdate_ );

      //important to set method to FixPoint, since we compute the RHS!!
      heatForm->SetNonLinMethod(FIXEDPOINT);
        
      //set the solution class to the operator
      heatForm->SetSolution( * sol_ );
    }

    heatForm->SetAnsatzFct( ansatzFct1_ );

    Matrix<Double> elemmat;
    heatForm->CalcElementMatrix(elemmat, ent, ent);

    // Get element solution
    Vector<Double> temp;
    sol_->GetElemSolution( temp, ent  );

    elemVec.Resize(numFncs);
    elemVec = -(elemmat * temp);

    //    std::cout << "RHS:\n" << "Temp.: " << temp << std::endl;
    //    std::cout << "elemVec:\n" << elemVec << std::endl;

    delete heatForm;
  }


}
