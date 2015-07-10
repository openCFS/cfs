// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <stddef.h>
#include <fstream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "General/environment.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "massInt.hh"

namespace CoupledField {
class BaseMaterial;
}  // namespace CoupledField

DECLARE_LOG(forms)


namespace CoupledField {



  MassInt::MassInt(const Double aDensity, const UInt nrDofsPerNode, bool axi, bool coordUpdate)
    : BaseForm(NULL)
  {
    density_ = aDensity;
    Init(nrDofsPerNode, axi, coordUpdate);

    LOG_DBG(forms) << "MassInt::MassInt() dens=" << density_ << " dofs=" << nrDofsPerNode << " axi=" << axi;
  }

  MassInt::MassInt(BaseMaterial* mat, const MaterialDescriptor& md, const UInt nrDofsPerNode, bool axi, bool coordUpdate)
    : BaseForm(mat)
  {
    md_ = md;
    assert(md_.mat_1 == DENSITY);
    density_ = md_.GetScalar(mat);

    Init(nrDofsPerNode, axi, coordUpdate);
    LOG_DBG(forms) << "MassInt::MassInt() md! dens=" << density_ << " dofs=" << nrDofsPerNode << " axi=" << axi;
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

    Double density = GetErsatzMaterialMass(ent1.GetElem(), direction);
    
    Double factor = mParser_->Eval( mHandle_ );

    LOG_DBG3(forms) << GetName() << "::CEM(" << ent1.GetElem()->elemNum << ") density=" 
                    << density << " factor=" << factor << " diagMass=" << diagMass_;

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
      MassMultiDof(multDofMass, elemMat);       
      elemMat = multDofMass;
    }

    LOG_DBG3(forms) << GetName() << "::CEM(" << ent1.GetElem()->elemNum << ") density=" 
                    << density << " factor=" << factor << " mat=" << elemMat.ToString();

  }

  void MassInt::MassMultiDof(Matrix<Double>& massMultDof, 
                             const Matrix<Double>& massMatSingleDof)
  {
    
    const UInt singleDofSize = massMatSingleDof.GetNumRows();

    UInt i, j, actDof;
    
    massMultDof.Resize(nrDofsPerNode_*singleDofSize);
    massMultDof.Init();
    
    for (i=0; i < singleDofSize; i++) {
      for (j=0; j < singleDofSize; j++) {
        double v = massMatSingleDof[i][j];
        for (actDof=0; actDof < nrDofsPerNode_; actDof++) {
          massMultDof[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof] = v; 
        }
      }
    }
  }

  double MassInt::GetErsatzMaterialMass(const Elem* elem, DesignElement::Type direction){
    // now check whether we have Param Mat, density is ignored then!!! 
    if(elem != NULL && domain->HasErsatzMaterialTensor()){
      assert(domain->GetErsatzBiMaterial(elem, md_.mat_class) == NULL);
      return domain->GetErsatzMaterial()->GetErsatzMaterialMass(elem, direction);
    }else{
      // for the harmonic topology optimzation case (or load ersatz material)
      // applies topology optimization ersatz material to the density, including bimaterial!
      Double density = md_.GetErsatzMaterial(this, elem, density_);
      LOG_DBG3(forms) << GetName() << "::CalcElementMatrix(" << elem->elemNum << ") -> density=" << density;    
      return density;
    }
  }
  
} // end namespace CoupledField



