// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <fstream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "General/environment.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/mathParser/mathParser.hh"
#include "linGradBDBInt.hh"
#include "Optimization/Design/DesignSpace.hh"

DECLARE_LOG(forms)

namespace CoupledField {


  // ============
  //   calcBMat
  // ============
  void linGradBDBInt::CalcBMat( Matrix<Double> &bMat, UInt ip,
			     const Matrix<Double> &ptCoord ) {


    // Obtain info on number of elements' funtions
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( dim_, numFncs );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->SetAnsatzFct( ansatzFct1_ );
    if (isSetIntPoint_) {
      ptelem->GetGlobDerivShFnc( xiDx, intPoint_, ptCoord, it1_.GetElem() );
    } else {
      ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );
    }

    if ( subTensorType_ == FULL ) {
      // The matrix bMat can be seen as a 3 x numFncs block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the BDB product evaluated at the k-th node of the finite
      // element. 
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
        bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else  {
      // The matrix bMat can be seen as a 1 x numFncs block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the ADB product evaluated at the k-th node of the finite
      // element. We assume that the first coordinate equals y and the
      // second z.
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
      }
    }


// Matrix<Double> auxbMat;
// auxbMat.Init();
// bMat.Transpose (auxbMat);
//     std::cerr << "linGradBDBInt::bMat transpose = \n" << auxbMat << std::endl;
  }


  // ============
  //   calcDMat
  // ============
  void linGradBDBInt::calcDMat( Matrix<Double> &dMat, const Elem* elem, DesignElement::Type direction, double force_factor)
  {
    // consider standard FEM, SIMP and FMO
    if(elem == NULL ||
        !domain->HasErsatzMaterialDielecTensor() ||
        !domain->GetErsatzMaterial()->GetDielecTensor(dMat, elem, direction))
    {
      ptMaterial->GetTensor(dMat,matType_,matDataType_,subTensorType_);
      dMat *= mParser_->Eval( mHandle_ );

      Double density = elem != NULL ? GetErsatzMaterialFactor(elem) : 1.0;
      if(density != 1.0)
      {
        dMat *= density;

        // BiMaterial case only valid for "simpVar" scheme (rho^param MAT + (1-rho)^param) BIMAT)
        BaseMaterial* bm = elem != NULL ? domain->GetErsatzBiMaterial(elem,  ELECTROSTATIC) : NULL;

        if(bm != NULL)
        {
          Double bidensity = GetErsatzMaterialFactor(elem, true);
          Matrix<Double> tmp;
          bm->GetTensor(tmp, matType_, matDataType_, subTensorType_);
          // tmp *= (1.0 - density);
          tmp *= bidensity;
          dMat +=  tmp;
          LOG_DBG3(forms) << "linGradBDBInt::calcDMat: e=" << elem->elemNum << " bimat=" << tmp.ToString();
        }
      }
      LOG_DBG3(forms) << GetName() << "::calcDMat(" << (elem != NULL ? Integer(elem->elemNum) : -1) << ") -> density=" << density;
    }
    LOG_DBG3(forms) << GetName() << "::calcDMat(" << (elem != NULL ? Integer(elem->elemNum) : -1) << ") -> " << dMat.ToString();
  }


  void linGradBDBInt::SetFactor( const std::string& factor ) {
    
    mParser_->SetExpr( mHandle_, factor );
  }
  
  // returns B - matrix
  void linGradBDBInt::calcBMatOnly( Matrix<Double> &bMat, UInt ip,
      BaseFE* elem, Matrix<Double> &ptCoord ) {

    ptelem = elem;
    CalcBMat(bMat, ip, ptCoord);
  }



  // ================
  //   Constructors
  // ================


  linGradBDBInt::linGradBDBInt(BaseMaterial* matData,
                               MaterialType matType,
                               SubTensorType type,
                               bool coordUpdate ) 
  : BDBInt(matData, type, coordUpdate ) {

    matType_ = matType;
    name_ = "linGradBDBInt";
    if ( type == FULL ) {
      dim_ = 3;
    }
    else {
      dim_ = 2;
    }

    if ( type == AXI ) {
      isaxi_     = true;
    }
  }


}
