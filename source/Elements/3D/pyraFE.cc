// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
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

  void PyraFE :: CoordsInsideElem(const Matrix<Double> & localCoords,
                                  const Double tolerance,
                                  StdVector<bool> & coordsInside) const
  {
    UInt numPoints = localCoords.GetSizeCol();
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

     
} // end of namespace
