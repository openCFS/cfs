// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "nLinMassInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField {



  NlinMassInt::NlinMassInt(const Double aDensity,  const UInt nrDofsPerNode, 
                   bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode)
      
  {
    name_ = "MassInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = MASS;
    isSolDependent_ = true;
  }

 
  NlinMassInt::~NlinMassInt()
  {
  }

 

  void NlinMassInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1, 
                                   EntityIterator& ent2  )  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;


    // get the variable speed of sound, stored in the solution object
    Matrix<Double> temp;
    Vector<Double> sos; // speed of sound
    Double val = 0.0;

    sol_->GetElemSolutionAsMatrix( temp, ent1 );
    temp.ConvertToVec_AppendRows( sos );
    for ( UInt i=0; i<sos.GetSize(); i++ ) {
       val += sos[i];
    }
    val /= (Double)sos.GetSize();

    Double factor;
    factor = (1.0/val) * (1.0/val);
//    std::cout << "massFactor=" << factor << std::endl;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize(nrNodes);
    elemMat.Init();
    

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        CoordAtIP = ptCoord_ * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * factor
	  * jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * density_ * factor * jacDet;
        
      elemMat += partElemMat;
    }
  
  }

} // end namespace CoupledField



