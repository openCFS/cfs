#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "hexaFE.hh"

namespace CoupledField
{

  HexaFE::HexaFE()
  {
    ENTER_FCN( "HexaFE::HexaFE" );

    Dim_        = 3;
    NumEdges_   = 12;
    NumFaces_   = 6;
    NumCorners_ = 8;
    numChilds_  = 8;
    MidPoint_   = 0.0, 0.0, 0.0;
    
    std::string integtype;
    params->Get( "type", integtype, "integRules", "hexa" );

    IntegType=String2EnumIntegrationType(integtype.c_str());

  }

  HexaFE::~HexaFE()
  {
    ENTER_FCN( "HexaFE::~HexaFE" );
    //  if (IntWeights) delete IntWeights;
  }

  void HexaFE::SetIntPoints()
  {
    ENTER_IFCN( "HexaFE::SetIntPoints" );
  
    Double factor = sqrt(3.0/5.0);

    switch(IntegType) 
      {
      case GaussOrder1:
          
        NumIntPoints_ = 1;
        DegreeInteg_  = 2;

        IntWeights_.Resize(NumIntPoints_);
        // all weights are 1.0
        for(UInt i=0; i<IntWeights_.GetSize(); i++)
          IntWeights_[i] = 1.0;
      
        if (!IntPoints_)
          IntPoints_ = new Vector<Double>[NumIntPoints_];
      
        for(UInt i=0; i<NumIntPoints_; i++)
          IntPoints_[i].Resize(Dim_);

        IntPoints_[0][0] = 0;
        IntPoints_[0][1] = 0;
        IntPoints_[0][2] = 0;      
      
        break;

      case GaussOrder2:

        NumIntPoints_=27;
        DegreeInteg_=3;

        IntWeights_.Resize(NumIntPoints_);
        IntWeights_[0]  = 0.171467764060357;
        IntWeights_[1]  = 0.274348422496571;
        IntWeights_[2]  = 0.171467764060357;
        IntWeights_[3]  = 0.274348422496571;
        IntWeights_[4]  = 0.438957475994513;
        IntWeights_[5]  = 0.274348422496571;
        IntWeights_[6]  = 0.171467764060357;
        IntWeights_[7]  = 0.274348422496571;
        IntWeights_[8]  = 0.171467764060357;
        IntWeights_[9]  = 0.274348422496571;
        IntWeights_[10] = 0.438957475994513;
        IntWeights_[11] = 0.274348422496571;
        IntWeights_[12] = 0.438957475994513;
        IntWeights_[13] = 0.702331961591221;
        IntWeights_[14] = 0.438957475994513;
        IntWeights_[15] = 0.274348422496571;
        IntWeights_[16] = 0.438957475994513;
        IntWeights_[17] = 0.274348422496571;
        IntWeights_[18] = 0.171467764060357;
        IntWeights_[19] = 0.274348422496571;
        IntWeights_[20] = 0.171467764060357;
        IntWeights_[21] = 0.274348422496571;
        IntWeights_[22] = 0.438957475994513;
        IntWeights_[23] = 0.274348422496571;
        IntWeights_[24] = 0.171467764060357;
        IntWeights_[25] = 0.274348422496571;
        IntWeights_[26] = 0.171467764060357;
  
        if (!IntPoints_)
          IntPoints_ = new Vector<Double>[NumIntPoints_];

        for(UInt i=0; i<NumIntPoints_; i++)
          IntPoints_[i].Resize(Dim_);

        IntPoints_[0][0] =  -factor;
        IntPoints_[0][1] =  -factor;
        IntPoints_[0][2] =  -factor;

        IntPoints_[1][0] =  -factor;
        IntPoints_[1][1] =  -factor;
        IntPoints_[1][2] =   0;

        IntPoints_[2][0] =  -factor;
        IntPoints_[2][1] =  -factor;
        IntPoints_[2][2] =   factor;

        IntPoints_[3][0] =  -factor;
        IntPoints_[3][1] =   0;
        IntPoints_[3][2] =  -factor;

        IntPoints_[4][0] =  -factor;
        IntPoints_[4][1] =   0;
        IntPoints_[4][2] =   0;

        IntPoints_[5][0] =  -factor;
        IntPoints_[5][1] =   0;
        IntPoints_[5][2] =   factor;

        IntPoints_[6][0] =  -factor;
        IntPoints_[6][1] =   factor;
        IntPoints_[6][2] =  -factor;

        IntPoints_[7][0] =  -factor;
        IntPoints_[7][1] =   factor;
        IntPoints_[7][2] =   0;

        IntPoints_[8][0] =  -factor;
        IntPoints_[8][1] =   factor;
        IntPoints_[8][2] =   factor;

        IntPoints_[9][0] =   0;
        IntPoints_[9][1] =  -factor;
        IntPoints_[9][2] =  -factor;

        IntPoints_[10][0] =   0;
        IntPoints_[10][1] =  -factor;
        IntPoints_[10][2] =   0;

        IntPoints_[11][0] =   0;
        IntPoints_[11][1] =  -factor;
        IntPoints_[11][2] =   factor;

        IntPoints_[12][0] =   0;
        IntPoints_[12][1] =   0;
        IntPoints_[12][2] =  -factor;


        IntPoints_[13][0] =   0;
        IntPoints_[13][1] =   0;
        IntPoints_[13][2] =   0;

        IntPoints_[14][0] =   0;
        IntPoints_[14][1] =   0;
        IntPoints_[14][2] =   factor;

        IntPoints_[15][0] =   0;
        IntPoints_[15][1] =   factor;
        IntPoints_[15][2] =  -factor;

        IntPoints_[16][0] =   0;
        IntPoints_[16][1] =   factor;
        IntPoints_[16][2] =   0;

        IntPoints_[17][0] =   0;
        IntPoints_[17][1] =   factor;
        IntPoints_[17][2] =   factor;

        IntPoints_[18][0] =   factor;
        IntPoints_[18][1] =  -factor;
        IntPoints_[18][2] =  -factor;

        IntPoints_[19][0] =   factor;
        IntPoints_[19][1] =  -factor;
        IntPoints_[19][2] =   0;

        IntPoints_[20][0] =   factor;
        IntPoints_[20][1] =  -factor;
        IntPoints_[20][2] =   factor;

        IntPoints_[21][0] =   factor;
        IntPoints_[21][1] =   0;
        IntPoints_[21][2] =  -factor;

        IntPoints_[22][0] =   factor;
        IntPoints_[22][1] =   0;
        IntPoints_[22][2] =   0;

        IntPoints_[23][0] =   factor;
        IntPoints_[23][1] =   0;
        IntPoints_[23][2] =   factor;

        IntPoints_[24][0] =   factor;
        IntPoints_[24][1] =   factor;
        IntPoints_[24][2] =  -factor;

        IntPoints_[25][0] =   factor;
        IntPoints_[25][1] =   factor;
        IntPoints_[25][2] =   0;

        IntPoints_[26][0] =   factor;
        IntPoints_[26][1] =   factor;
        IntPoints_[26][2] =   factor;

        break;
          
        //      NumIntPoints_=6;
        //       DegreeInteg_=3;


        //       IntWeights_.Resize(NumIntPoints_);
        //      // all weights are 4/3
        //       for(UInt i=0; i<IntWeights_.GetSize(); i++)
        //      IntWeights_[i] = 4.0/3.0;
      
        //       if (!IntPoints_)
        //      IntPoints_ = new Vector<Double>[NumIntPoints_];

        //       for(UInt i=0; i<NumIntPoints_; i++)
        //      IntPoints_[i].Resize(Dim_);

        //       IntPoints_[0][0] =  1.0/sqrt(6.0);
        //       IntPoints_[0][1] =  1.0/sqrt(2.0);
        //       IntPoints_[0][2] = -1.0/sqrt(3.0);

        //       IntPoints_[1][0] =  1.0/sqrt(6.0);
        //       IntPoints_[1][1] = -1.0/sqrt(2.0);
        //       IntPoints_[1][2] = -1.0/sqrt(3.0);
      
        //       IntPoints_[2][0] = -1.0/sqrt(6.0);
        //       IntPoints_[2][1] =  1.0/sqrt(2.0);
        //       IntPoints_[2][2] =  1.0/sqrt(3.0);

        //       IntPoints_[3][0] = -1.0/sqrt(6.0);
        //       IntPoints_[3][1] = -1.0/sqrt(2.0);
        //       IntPoints_[3][2] =  1.0/sqrt(3.0);

        //       IntPoints_[4][0] = -sqrt(2.0/3.0);
        //       IntPoints_[4][1] = 0.0;
        //       IntPoints_[4][2] = -1/sqrt(3.0);

        //       IntPoints_[5][0] = sqrt(2.0/3.0);
        //       IntPoints_[5][1] = 0.0;
        //       IntPoints_[5][2] = 1/sqrt(3.0);


        //       NumIntPoints_=8;
        //       DegreeInteg_=3;


        //      IntWeights_.Resize(NumIntPoints_);
        //      // all weights are 1.0
        //      for(UInt i=0; i<IntWeights_.GetSize(); i++)
        //        IntWeights_[i] = 1.0;
      
        //      if (!IntPoints_)
        //        IntPoints_ = new Vector<Double>[NumIntPoints_];

        //      for(UInt i=0; i<NumIntPoints_; i++)
        //        IntPoints_[i].Resize(Dim_);

        //       IntPoints_[0][0] = -0.57735026919;
        //       IntPoints_[0][1] = -0.57735026919;
        //       IntPoints_[0][2] = -0.57735026919;
      
        //       IntPoints_[1][0] = 0.57735026919;
        //       IntPoints_[1][1] = -0.57735026919;
        //       IntPoints_[1][2] = -0.57735026919;

        //       IntPoints_[2][0] = 0.57735026919;
        //       IntPoints_[2][1] = 0.57735026919;
        //       IntPoints_[2][2] = -0.57735026919;

        //       IntPoints_[3][0] = -0.57735026919;
        //       IntPoints_[3][1] = 0.57735026919;
        //       IntPoints_[3][2] = -0.57735026919;

        //       IntPoints_[4][0] = -0.57735026919;
        //       IntPoints_[4][1] = -0.57735026919;
        //       IntPoints_[4][2] = 0.57735026919;
      
        //       IntPoints_[5][0] = 0.57735026919;
        //       IntPoints_[5][1] = -0.57735026919;
        //       IntPoints_[5][2] = 0.57735026919;

        //       IntPoints_[6][0] = 0.57735026919;
        //       IntPoints_[6][1] = 0.57735026919;
        //       IntPoints_[6][2] = 0.57735026919;

        //       IntPoints_[7][0] = -0.57735026919;
        //       IntPoints_[7][1] = 0.57735026919;
        //       IntPoints_[7][2] = 0.57735026919;

        //      break;

      default:
        Error("Integration type is not implemented",__FILE__,__LINE__);
      }
  }

  void HexaFE::GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                         const StdVector<UInt> & volConnect,
                                         const Vector<Double> & surfIntPoint,
                                         Vector<Double> & volIntPoint)
  {
    ENTER_IFCN( "HexaFE::GetLocalIntPoints4Surface" );
  
    // Try to find out, which vertices are in common with
    // the surface element. Then calculate the product of all four
    // and compare them
    //
    //                    zeta 
    //                     ^ eta 
    //    8 +-------+ 7    |/
    //     /|      /|      0--> xi 
    //    / |     / |
    // 5 +--+----+6 |   
    //   |  +-- -|- + 3    
    //   | / 4   | /    REFERENCE VOLUME ELEMENT
    //   |/      |/
    // 1 +-------+ 2



    StdVector<UInt> commonIndex(4);
    UInt found = 0;
    UInt indexProduct = 0;
    std::string errMsg;
  
    volIntPoint.Resize(3);
  
    // loop over surface connect
    for (UInt iSurf=0; iSurf<4; iSurf++)
      // loop over volume connect
      for (UInt iVol=0; iVol<8; iVol++)
        if (surfConnect[iSurf] == volConnect[iVol])
          {
            commonIndex[found++] = iVol+1;
          }

    // std::cerr << std::endl << std::endl;
    //std::cerr << "commonIndex = " << std::endl << commonIndex << std::endl << std::endl;
    indexProduct =  commonIndex[0] * commonIndex[1];
    indexProduct *= commonIndex[2] * commonIndex[3];

    //std::cerr << "indexProduct = " << indexProduct << std::endl;
    switch(indexProduct)
      {
      case 24:
        // Surface[1,2,3,4] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = surfIntPoint[1];
        volIntPoint[2] = -1.0;
        break;

      case 1680:
        // Surface[5,6,7,8] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = surfIntPoint[1];
        volIntPoint[2] = 1.0;
        break;

      case 252:
        // Surface[2,3,7,6] is common
        volIntPoint[0] = 1.0;
        volIntPoint[1] = surfIntPoint[0];
        volIntPoint[2] = surfIntPoint[1];
        break;
      
      case 160:
        // Surface[1,5,8,4] is common
        volIntPoint[0] = -1.0;
        volIntPoint[1] = surfIntPoint[0];
        volIntPoint[2] = surfIntPoint[1];
        break;
      
      case 672:
        // Surface[4,3,7,8] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = 1.0;
        volIntPoint[2] = surfIntPoint[1];
        break;
      
      case 60:
        // Surface[1,2,6,5] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = -1.0;
        volIntPoint[2] = surfIntPoint[1];
        break;
      
      default:
        errMsg = "HexaFE::GetLocalIntPoints4Surface: surface and volume element ";
        errMsg = "have not four nodes in common. Check your .mesh-file.";
        Error(errMsg.c_str(), __FILE__, __LINE__);
      }
  }


  void HexaFE::GetMaxMinEdgeLength( Matrix<Double> &ptCoord, Double &Lmax, Double &Lmin )
  {
    ENTER_IFCN( "HexaFE::GetMaxMinEdgeLength" );

    Vector<Double> L(12);
    Double dx,dy,dz;

    //first face
    for (UInt i=0; i<3; i++ ) {
      dx   = ptCoord[0][i] - ptCoord[0][i+1];
      dy   = ptCoord[1][i] - ptCoord[1][i+1];
      dz   = ptCoord[2][i] - ptCoord[2][i+1];
      L[i] = sqrt( dx*dx + dy*dy + dz*dz);
    }
    dx   = ptCoord[0][0] - ptCoord[0][3];
    dy   = ptCoord[1][0] - ptCoord[1][3];
    dz   = ptCoord[2][0] - ptCoord[2][3];
    L[3] = sqrt( dx*dx + dy*dy + dz*dz);

    //second face
    for (UInt i=4; i<7; i++ ) {
      dx   = ptCoord[0][i] - ptCoord[0][i+1];
      dy   = ptCoord[1][i] - ptCoord[1][i+1];
      dz   = ptCoord[2][i] - ptCoord[2][i+1];
      L[i] = sqrt( dx*dx + dy*dy + dz*dz);
    }
    dx   = ptCoord[0][4] - ptCoord[0][7];
    dy   = ptCoord[1][4] - ptCoord[1][7];
    dz   = ptCoord[2][4] - ptCoord[2][7];
    L[7] = sqrt( dx*dx + dy*dy + dz*dz);

    // side faces
    UInt offset = 8;
    for (UInt i=0; i<4; i++ ) {
      dx   = ptCoord[0][i] - ptCoord[0][i+4];
      dy   = ptCoord[1][i] - ptCoord[1][i+4];
      dz   = ptCoord[2][i] - ptCoord[2][i+4];
      L[i+offset] = sqrt( dx*dx + dy*dy + dz*dz);
    }

    Lmax = L[0];
    Lmin = Lmax; 
    for ( UInt i=1; i<12; i++ ) {
      if ( L[i] > Lmax ) {
	Lmax = L[i];
      }
      if ( L[i] < Lmin ) {
	Lmin = L[i];
      }
    }
  }

} // end of namespace

