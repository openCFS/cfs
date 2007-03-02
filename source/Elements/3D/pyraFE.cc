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
    ENTER_FCN( "PyraFE::PyraFE" );

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
    ENTER_FCN( "PyraFE::~PyraFE" );
  }


 void PyraFE::FillIntegrationPoints()
 {
      ENTER_IFCN("PyraFE::FillIntegrationPoints");      

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
     
} // end of namespace
