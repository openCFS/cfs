#include <iostream>
#include <fstream>

#include "laplaceInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField
{


  LaplaceInt::LaplaceInt(Double aVal, bool axi, bool coordUpdate )
    : BaseForm(NULL),laplVal_ (aVal)
  {
    ENTER_FCN( "LaplaceInt::LaplaceInt" );

    name_ = "LaplaceInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
  }


 
  LaplaceInt::~LaplaceInt()
  {
    ENTER_FCN( "LaplaceInt::~LaplaceInt" );
  }



  void LaplaceInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 ) {
  
    ENTER_FCN( "LaplaceInt::CalcElementMatrix" );
   
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
   
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes); 
    elemMat.Init();

    //     //check for material value
    //     if (materialArray_ != NULL) {
    //       laplVal_ = (*materialArray_)[actSD_][actElemNr_];
    //     }


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);

        xiDx.Transpose(xiDxTransp);

        partElemMat = xiDx * xiDxTransp;
        
        if (isaxi_)
          {
            ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
            CoordAtIP = ptCoord_ * ShpFncAtIp;
            partElemMat *= 2 * PI * intWeights[actIntPt-1] 
              * jacDet * laplVal_ * CoordAtIP[0];

          }
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * laplVal_;

        elemMat += partElemMat;
      }

    //    std::cout << "ElemMatLaplace:\n" << elemMat << std::endl;
  }



} // end namespace CoupledField
