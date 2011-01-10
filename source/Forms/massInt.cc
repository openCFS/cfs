// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "massInt.hh"
#include "Domain/domain.hh"
#include "Optimization/Design/DesignSpace.hh"

#include "DataInOut/Logging/cfslog.hh"

DECLARE_LOG(forms)


namespace CoupledField {



  MassInt::MassInt(const Double aDensity, const UInt nrDofsPerNode, bool axi, bool coordUpdate)
    : BaseForm(NULL)
  {
    density_ = aDensity;
    Init(nrDofsPerNode, axi, coordUpdate);
  }

  MassInt::MassInt(BaseMaterial* mat, const MaterialDescriptor& md, const UInt nrDofsPerNode, bool axi, bool coordUpdate)
    : BaseForm(mat)
  {
    md_ = md;
    assert(md_.mat_1 == DENSITY);
    density_ = md_.GetScalar(mat);

    Init(nrDofsPerNode, axi, coordUpdate);
  }



  MassInt::~MassInt() {
  }


  void MassInt::Init(const UInt nrDofsPerNode, bool axi, bool coordUpdate)
  {
    nrDofsPerNode_ = nrDofsPerNode;
    diagMass_ = false;
    name_ = "MassInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = MASS;
    if ( coordUpdate ) 
      isSolDependent_ = true;
  }

  void MassInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1, 
                                   EntityIterator& ent2,
                                   DesignElement::Type direction)  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();

    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(numFncs);
    elemMat.Resize(numFncs);
    elemMat.Init();

    // applies topology optimization ersatz material to the density, including bimaterial!
    Double density = md_.GetErsatzMaterial(this, ent1.GetElem(), density_);
    // ErsatzMaterialMass is Bastian's stuff and not compatible with MaterialDescriptor::GetErsatzMaterial()
    if(domain->HasErsatzMaterialMass())
    {
      assert(domain->GetErsatzBiMaterial(ent1.GetElem(), md_.mat_class) == NULL);
      Double top_opt = domain->GetErsatzMaterial()->GetErsatzMaterialMass(ent1.GetElem(), direction);
      density *= top_opt;
    }

    Double factor = mParser_->Eval( mHandle_ );

    LOG_DBG3(forms) << GetName() << "::CEM(" << ent1.GetElem()->elemNum << ") density=" << density << " factor=" << factor;

    if (diagMass_ ) {
      Double mass = 0.0;
      for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );

        if (isaxi_) {
          ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
          CoordAtIP = ptCoord_ * shapeFncAtIp;
          mass += 2 * PI * intWeights[actIntPt-1] * density
                 * factor * jacDet * CoordAtIP[0];
        }
        else 
          mass += intWeights[actIntPt-1] * density * factor * jacDet;
      }

      for ( UInt i=0; i<numFncs; i++ ) 
        elemMat[i][i] = mass / (Double)numFncs;
    }
    else {
      for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );

        ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );

        partElemMat.DyadicMult(shapeFncAtIp);

        if (isaxi_) {
          CoordAtIP = ptCoord_ * shapeFncAtIp;
          partElemMat *= 2 * PI * intWeights[actIntPt-1] * density * factor* jacDet * CoordAtIP[0];
        }
        else 
          partElemMat *= intWeights[actIntPt-1] * density * factor * jacDet;

        elemMat += partElemMat;
      }
    }

    if (nrDofsPerNode_ > 1 ) {
      Matrix <Double> multDofMass;
      MassMultiDof(multDofMass, elemMat, nrDofsPerNode_);       
      elemMat = multDofMass;
    }

    LOG_DBG3(forms) << GetName() << "::CEM(" << ent1.GetElem()->elemNum << ") density=" << density << " factor=" << factor << " mat=" << elemMat.ToString();

  }

  void MassInt::MassMultiDof(Matrix<Double>& massMultDof, 
                             const Matrix<Double>& massMatSingleDof,  const UInt nrDofs)
  {
    
    const UInt singleDofSize = massMatSingleDof.GetNumRows();

    UInt i, j, actDof;
    
    massMultDof.Resize(nrDofsPerNode_*singleDofSize);
    massMultDof.Init();
    
    for (i=0; i < singleDofSize; i++) {
      for (j=0; j < singleDofSize; j++) {
        double v = massMatSingleDof[i][j];
        for (actDof=0; actDof < nrDofs; actDof++) {
          massMultDof[i*nrDofs + actDof][j*nrDofs + actDof] = v; 
        }
      }
    }
  }


  void MassInt::MassMultiDofZero(Matrix<Double>& massMultDofZero, const Matrix<Double>& massMatSingleDof)
  {
    
    const UInt singleDofSize = massMatSingleDof.GetNumRows();
    
    UInt i, j, actDof;
    
    massMultDofZero.Resize(nrDofsPerNode_*singleDofSize);
    massMultDofZero.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
        for (actDof=0; actDof < nrDofsPerNode_ ; actDof++)
          {
            massMultDofZero[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof] = 
              massMatSingleDof[i][j]; 
          }
  }
  
  double MassInt::GetErsatzMaterialMass(const Elem* elem, DesignElement::Type direction){
    // now check whether we have Param Mat, density is ignored then!!! 
    if(elem != NULL && domain->HasErsatzMaterialMass()){
      return domain->GetErsatzMaterial()->GetErsatzMaterialMass(elem, direction);
    }else{
      // for the harmonic topology optimzation case (or load ersatz material)
      double density = GetErsatzMaterialFactor(elem);
      LOG_DBG3(forms) << GetName() << "::CalcElementMatrix(" << elem->elemNum << ") -> density=" << density;    
      return density;
    }
  }
  
} // end namespace CoupledField



