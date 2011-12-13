// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <string>

#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "pyraFE.hh"

namespace CoupledField
{


  PyraFE::PyraFE()
  {

    Dim_        = 3;
    NumEdges_   = 8;
    NumFaces_   = 5;
    NumCorners_ = 5;
    //    numChilds_ = 8;
    MidPoint_.Resize(3);
    MidPoint_[0] = 0.0;
    MidPoint_[1] = 0.0;
    MidPoint_[2] = 1./4;

  }
  
  PyraFE::~PyraFE()
  {
  }


 void PyraFE::FillIntegrationPoints()
 {

      //"Pyramidal Elements"
      // F. Zgainski, J.L. Coulomb, Y. Marechal. IEEE Transactions on Magnetics, Vol. 32, No. 3, May 1996
      // converted from the old setIntPoints()
      static Double c2[][4] = { 
        {-0.506616303350116,  -0.506616303350116,  0.122514822655441,  0.100785882079825},
        {0.506616303350116,  -0.506616303350116,  0.122514822655441,  0.100785882079825},
        {0.506616303350116,  0.506616303350116,  0.122514822655441,  0.100785882079825},
        {-0.506616303350116,  0.506616303350116,  0.122514822655441,  0.100785882079825},
        {-0.263184055569884,  -0.263184055569884,  0.544151844011225,  0.232547451253508},
        {0.263184055569884,  -0.263184055569884,  0.544151844011225,  0.232547451253508},
        {0.263184055569884,  0.263184055569884,  0.544151844011225,  0.232547451253508},
        {-0.263184055569884,  0.263184055569884,  0.544151844011225,  0.232547451253508},
      };
      AddIntegrationPoints(CLASSICAL, 2, 8, (Double*) c2);
      
  }

  void PyraFE::CoordsInsideElem(const Matrix<Double> & localCoords,
                                  const Double tolerance,
                                  StdVector<bool> & coordsInside) const
  {
    UInt numPoints = localCoords.GetNumCols();
    double xi, eta, zeta;
    double threshold;

    coordsInside.Resize(numPoints);
    
    for(UInt i=0; i<numPoints; i++)
    {
      xi = localCoords[0][i];
      eta = localCoords[1][i];
      zeta = localCoords[2][i];
      threshold = 1 - zeta;
      

      coordsInside[i] = (zeta >= (0 - tolerance)) &&
                        (zeta <= (1.0 + tolerance)) &&
                        (  xi >= (-threshold - tolerance)) &&
                        ( eta >= (-threshold - tolerance)) &&
                        (  xi <= (threshold + tolerance)) &&
                        ( eta <= (threshold + tolerance));  
    }
  }
  
  void PyraFE::GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                         const StdVector<UInt> & volConnect,
                                         const Vector<Double> & surfIntPoint,
                                         Vector<Double> & volIntPoint)
  {

    // Try to find out, which vertices are in common with
    // the surface element. Then calculate the product of all four
    // and compare them
    //                                   zeta
    //             5                     ^  eta
    //             +                     |/
    //           // \                    0--> xi
    //          // \ \                      
    //         / / \  \                    
    //        / /   \  \
    //     2 +-/----\ ---+ 1
    //      / /     \   /
    //     / /      \  /     REFERENCE VOLUME ELEMENT
    //    //        \ /
    //  3+-----------+ 4
    //        

    // Check if surface element is triangle 
    // or quadrilateral
    if (surfConnect.GetSize() == 3 ||
        surfConnect.GetSize() == 6) 
    {
      // ---- Triangle Surface ---
      StdVector<Integer> commonIndex(3);
      Integer found = 0;
      Integer indexProduct = 0;
      std::string errMsg;

      volIntPoint.Resize(3);

      // loop over surface connect
      for (Integer iSurf=0; iSurf<3; iSurf++)
        // loop over volume connect
        for (Integer iVol=0; iVol<5; iVol++)
          if (surfConnect[iSurf] == volConnect[iVol]) {
            commonIndex[found++] = iVol+1;
          }
      indexProduct =  commonIndex[0] * commonIndex[1] * commonIndex[2];

      // Now we have to consider the following:
      // - The extension of the triangular element is from [0..1] in both
      //   local directions xi and eta
      // - The side length of the base-rectangular side is in each direction 
      //   [-1..+1], so we need a mapping in 
      // - The 

      // General rule:


      switch( indexProduct ) {
        case 10:
          // Surface[1,2,5] is common
          volIntPoint[0] = - 2.0 * (surfIntPoint[0] - 0.5);
          volIntPoint[1] =   1.0 - surfIntPoint[1];
          volIntPoint[2] =   surfIntPoint[1];
          break;
        case 30:
          // Surface[2,3,5] is common
          volIntPoint[0] =   surfIntPoint[1] - 1.0;
          volIntPoint[1] = - 2.0 * (surfIntPoint[0] - 0.5);
          volIntPoint[2] =   surfIntPoint[1];
          break;
        case 60:
          // Surface[3,4,5] is common
          volIntPoint[0] =   2.0 * (surfIntPoint[0] - 0.5);
          volIntPoint[1] =   surfIntPoint[1] - 1.0;
          volIntPoint[2] =   surfIntPoint[1];
          break;
        case 20:
          // Surface[4,1,5] is common
          volIntPoint[0] =  1.0 - surfIntPoint[1];
          volIntPoint[1] =  2.0 * (surfIntPoint[0] - 0.5);
          volIntPoint[2] =  surfIntPoint[1];
          break;
        default:
          EXCEPTION("PyraFE::GetLocalIntPoints4Surface: surface and volume element "
              << "have not three nodes in common. Check your mesh.");
      }
    } else {
      // ---- Quadrilateral Surface ---
      StdVector<Integer> commonIndex(4);
      Integer found = 0;
      Integer indexSum = 0;
      std::string errMsg;

      volIntPoint.Resize(3);

      // loop over surface connect
      for (Integer iSurf=0; iSurf<4; iSurf++)
        // loop over volume connect
        for (Integer iVol=0; iVol<5; iVol++)
          if (surfConnect[iSurf] == volConnect[iVol])
          {
            commonIndex[found++] = iVol+1;
          }
      indexSum =  commonIndex[0] + commonIndex[1] 
                                               + commonIndex[2] + commonIndex[3];

      // Safety check: Check, that the surface element is
      // really located on the bottom of the pyramid.
      if( indexSum != 10 ) {
        EXCEPTION("PyraFE::GetLocalIntPoints4Surface: surface and volume element "
            << "have not four nodes in common. Check your mesh.");
      }

      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = surfIntPoint[2];
      volIntPoint[2] = 0.0; // always on bottom
    }
  }  
} // end of namespace
