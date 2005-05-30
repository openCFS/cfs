#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
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
    MidPoint_ = 0.0, 0.0, 1./5;
  
    std::string integtype;
    params->Get( "type", integtype, "integRules", "pyra" );
  
    IntegType = String2EnumIntegrationType(integtype.c_str());
  
    //  isSetAtCenter_=FALSE;
  }
  
  PyraFE::~PyraFE()
  {
    ENTER_FCN( "PyraFE::~PyraFE" );
  }

  void PyraFE::SetIntPoints()
  {
    ENTER_IFCN( "PyraFE::SetIntPoints" );

    switch(IntegType)
      {
      case GaussOrder1:

        Error("Integration type GaussOrder1 is not implemented for pyramids",__FILE__,__LINE__);
        break;

      case GaussOrder2:
        //"Pyramidal Elements"
        // F. Zgainski, J.L. Coulomb, Y. Marechal. IEEE Transactions on Magnetics, Vol. 32, No. 3, May 1996
        NumIntPoints_ = 8;
        DegreeInteg_  = 3;

        IntWeights_.Resize(NumIntPoints_);
        // weights are different for the lower and upper integration points
        for(UInt i=0; i<IntWeights_.GetSize(); i++)
          {
            if(i<(IntWeights_.GetSize()/2))
              {
                IntWeights_[i] = 0.100785882079825;
              }
            else
              {
                IntWeights_[i] = 0.232547451253508;
              }
          }
        
        if (!IntPoints_)
          IntPoints_ = new Vector<Double>[NumIntPoints_];

        for(UInt i=0; i<NumIntPoints_; i++)
          IntPoints_[i].Resize(Dim_);

        IntPoints_[0][0] = -0.506616303350116;
        IntPoints_[0][1] = -0.506616303350116;
        IntPoints_[0][2] = 0.122514822655441;
        
        IntPoints_[1][0] = 0.506616303350116;
        IntPoints_[1][1] = -0.506616303350116;
        IntPoints_[1][2] = 0.122514822655441;

        IntPoints_[2][0] = 0.506616303350116;
        IntPoints_[2][1] = 0.506616303350116;
        IntPoints_[2][2] = 0.122514822655441;
        
        IntPoints_[3][0] = -0.506616303350116;
        IntPoints_[3][1] = 0.506616303350116;
        IntPoints_[3][2] = 0.122514822655441;
        
        IntPoints_[4][0] = -0.263184055569884;
        IntPoints_[4][1] = -0.263184055569884;
        IntPoints_[4][2] = 0.544151844011225;

        IntPoints_[5][0] = 0.263184055569884;
        IntPoints_[5][1] = -0.263184055569884;
        IntPoints_[5][2] = 0.544151844011225;

        IntPoints_[6][0] = 0.263184055569884;
        IntPoints_[6][1] = 0.263184055569884;
        IntPoints_[6][2] = 0.544151844011225;

        IntPoints_[7][0] = -0.263184055569884;
        IntPoints_[7][1] = 0.263184055569884;
        IntPoints_[7][2] = 0.544151844011225;
        


        break;
      
      default:
        Error("Integration type is not implemented",__FILE__,__LINE__);
      }  
  }


} // end of namespace
