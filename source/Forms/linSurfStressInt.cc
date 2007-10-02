// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "linSurfStressInt.hh"

#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField {
  
  
  SurfStress3DLinForm::SurfStress3DLinForm( Vector<Double> stressVec,
                                            shared_ptr<SurfElemList> surfElems )
    : LinearForm(), stress_(stressVec) {
    name_ = "SurfStress3DLinForm";

    isaxi_ = false;

    // Iterate over all elements in region
    EntityIterator it = surfElems->GetIterator();
    
    for( it.Begin(); !it.IsEnd(); it++ ) {
      SurfElem const * surfEl = it.GetSurfElem();
      //std::cerr << "Assigning volume element" << surfEl->ptVolElem1->elemNum
      //          << " surface elem" << surfEl->elemNum << std::endl;
      surfElems_[surfEl->ptVolElem1] = surfEl;
    }
      
  }
    
  
  
  SurfStress3DLinForm::~SurfStress3DLinForm() {

  }
  
  
  void SurfStress3DLinForm::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent ) {

    // Extract pointer to reference element and get coordinates
    ptelem = ent.GetElem()->ptElem;
   
    // Get related surface element
    if( surfElems_.find(ent.GetElem() ) == surfElems_.end() ) {
      *error << "Could not find matching surface element for "
             << "volume element " << ent.GetElem()->elemNum;
      Error( __FILE__, __LINE__ );
    }

    SurfElem const * surf = surfElems_[ent.GetElem()];
    BaseFE  * ptSurfFE = surf->ptElem;
    
    domain->GetGrid()->GetElemNodesCoord( ptCoord_, surf->connect, 
                                          coordUpdate_ );


    
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    ptSurfFE->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptSurfFE->GetNumIntPoints();
    const Vector<Double> & intWeights = ptSurfFE->GetIntWeights();
    Vector<Double> shapeFnc, normal, locVolCoord, elemStress(3);
    Vector<Double> * intPoints = ptSurfFE->GetIntPoints();

    elemVec.Resize(numFncs*3);
    elemVec.Init(0);
    //std::cerr << "Number of cuntions: " << numFncs << std::endl;
    // Calculate normal of surface element
    domain->GetGrid()->CalcSurfNormal(normal,*surf );
    normal *= (Double) surf->normalSign;

    // calculate element-local stress vector
    elemStress[0] = stress_[0] * normal[0] + stress_[5]*normal[1] 
                  + stress_[4] * normal[2];
    elemStress[1] = stress_[5] * normal[0] + stress_[1]*normal[1] 
                  + stress_[3] * normal[2];
    elemStress[2] = stress_[4] * normal[0] + stress_[3]*normal[1] 
                  + stress_[2] * normal[2];

    
    Double area = 0;
    
     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

       // Call following method with modified coordinates
       ptelem->GetLocalIntPoints4Surface(surf->connect, 
                                         ent.GetElem()->connect,
                                         intPoints[actIntPt-1],
                                         locVolCoord );
       //std::cerr << "2D coord: " << intPoints[actIntPt-1].Serialize() << std::endl;
       //std::cerr << "3D coord: " << locVolCoord.Serialize() << "\n\n";

       
       ptelem->GetShFnc( shapeFnc, locVolCoord, ent.GetElem() );
       //std::cerr << "shapeFnc =\n" << shapeFnc << std::endl;
       Double jacDet = ptSurfFE->CalcJacobianDetAtIp(actIntPt, ptCoord_,
                                                     surf );
       Double factor = intWeights[actIntPt-1] * jacDet;
       shapeFnc *= factor;    
       area += (jacDet * intWeights[actIntPt-1]);
    
       // now fill in entries for RHS vector
       for (UInt iFcn = 0; iFcn < numFncs; iFcn++)  {
         for( UInt iDof = 0; iDof < 3; iDof++ ) { 
           //std::cerr << "Writing to entry " << iFcn*3+iDof << std::endl;
           elemVec[iFcn*3+iDof] += shapeFnc[iFcn] * elemStress[iDof];
        
         } // dofs
       } // functions
     } // integration points

//      std::cerr << "area is " << area << std::endl;
//      std::cerr << "ElemVec = \n" << elemVec << std::endl;
    
     // Calculate force value in the middle
     Vector<Double> locCoord(3);
     locCoord.Init();
     
     
     Double sum = 0.0;
     Vector<Double> helpVec;
     
     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
       ptelem->GetLocalIntPoints4Surface(surf->connect, 
                                         ent.GetElem()->connect,
                                         intPoints[actIntPt-1],
                                         locVolCoord );
       
       ptelem->GetShFnc( shapeFnc, locVolCoord, ent.GetElem() );
       helpVec.Resize(elemVec.GetSize()) ;
       helpVec.Init();
       for (UInt iFcn = 0; iFcn < numFncs; iFcn++)  {
         for( UInt iDof = 0; iDof < 3; iDof++ ) { 
           helpVec[iFcn*3+iDof] = shapeFnc[iFcn];
         } 
       }        
       Double part = helpVec * elemVec;
       //std::cerr << "value at " << locVolCoord.Serialize();
       //std::cerr << " is " << part << std::endl;
       
       Double jacDet = ptSurfFE->CalcJacobianDetAtIp(actIntPt, ptCoord_,
                                                     surf );
       //std::cerr << "jacDet = " << jacDet << std::endl;
       Double factor = intWeights[actIntPt-1] * jacDet;
       sum += part * factor;
     }
       
     //     std::cerr << "helpVec = \n" << helpVec << std::endl;
  
     
     //Double midVal = helpVec * elemVec;

//      std::cerr << "Total value of force is " << sum << "[N]\n";
     //std::cerr << "\n\nValue at (" << locVolCoord.Serialize();


 
    
    
  }
  
}// end of namespace
