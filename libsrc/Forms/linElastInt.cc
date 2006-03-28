#include <iostream>
#include <fstream>

#include "linElastInt.hh"

namespace CoupledField
{

  // =====================
  //   CalcElementMatrix
  // =====================
  void linElastInt::CalcElementMatrix( Matrix<Double> &ptCoord,
				       Matrix<Double> &elemMat ) {

    ENTER_FCN( "linElastInt::CalcElementMatrix" );

    Boolean Softening;

    if (softeningModel_=="BK1") {
      softeningPart_ = "bendingBK1";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      Matrix<Double> helpMat = elemMat;
      // std::cout << "Bending Mat:\n" << helpMat << std::endl;

      softeningPart_ = "shearBK1";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      elemMat += helpMat;
      //     std::cout << "Total Mat:\n" << elemMat << std::endl;
    }
    else if (softeningModel_=="SRI") {
      softeningPart_ = "bendingSRI";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      Matrix<Double> helpMat = elemMat;
      // std::cout << "Bending Mat:\n" << helpMat << std::endl;

      softeningPart_ = "shearSRI";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      elemMat += helpMat;
      //     std::cout << "Total Mat:\n" << elemMat << std::endl;
    }
    else {
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
    }
  }

  // returns B - matrix for BDB
  void linElastInt::calcBMat( Matrix<Double> &bMat, UInt ip,
                              Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linElastInt::calcBMat" );

    const UInt nrNodes  = ptelem->GetNumNodes();
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  

    UInt actDim, actNode, j, k;
    
    
    bMat.Resize(getDimD(), nrNodes * nrDofs);
    
    // local shape functions derived after global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;

    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord);
    else
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);


    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < nrNodes; actNode++)
        bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
        j = 1;
        k = 0;
        
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
            bMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
          }

        if (isaxi_)
          {
            UInt idxtheta = getDimD();
            Vector<Double> ShpFncAtIp;
            Vector<Double> CoordAtIP;

            if (isSetIntPoint_) 
              ptelem->GetShFnc(ShpFncAtIp,intPoint_);
            else
              ptelem->GetShFncAtIp(ShpFncAtIp,ip);

            CoordAtIP = ptCoord * ShpFncAtIp;

            for (actNode = 0; actNode < nrNodes; actNode++)          
              bMat[idxtheta-1][actNode * spaceDim] =
                ShpFncAtIp[actNode] / CoordAtIP[0];
          }

        break;


      case 3:
        UInt actDim=spaceDim;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
          }

        actDim++;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
          }

        actDim++;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
          }
        break;
      }

    isSetIntPoint_ = FALSE;
  }

  void linElastInt:: calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linElastInt::calcDMat" );
    
    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,REAL, subTensorType_);

     //check for softening model
    if ( subTensorType_ == AXI ) {
      if (softeningPart_ == "bendingBK1" ) {
	Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
      }
      
      else if (softeningPart_ == "bendingSRI" ) {
	UInt idx = dimD_-2;
	dMat[idx][idx] = 0.0;
      }
      
      else if (softeningPart_ == "shearBK1" ) {
	Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
      }
      
      else if (softeningPart_ == "shearSRI" ) {
	for (UInt i=0; i<dimD_-2; i++) {
	  for (UInt j=0; j<dimD_-2; j++) {
	    dMat[i][j] = 0.0;
	  }
	}
	
	UInt rowidx=dimD_-1;
	for (UInt i=0; i<dimD_; i++) {
	  dMat[rowidx][i] = 0.0;
	}
	
	UInt colidx=dimD_-1;
	for (UInt i=0; i<dimD_; i++) {
	  dMat[i][colidx] = 0.0;
	}

      }
    }

    else if ( subTensorType_ ==  PLANE_STRAIN || subTensorType_ ==  PLANE_STRESS ) {
      //check for softening model
      if (softeningPart_ == "bendingBK1" ) {
	Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
      }
      
      else if (softeningPart_ == "bendingSRI" ) {
	UInt idx = dimD_-1;
	dMat[idx][idx] = 0.0;
      }
      
      else if (softeningPart_ == "shearBK1" ) {
	Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
      }
      
      else if (softeningPart_ == "shearSRI" ) {
	for (UInt i=0; i<dimD_; i++) {
	  for (UInt j=0; j<dimD_; j++) {
	    dMat[i][j] = 0.0;
	  }
	}
      }
    }
  
    else if (subTensorType_ == FULL ) {
      if (softeningPart_ == "bendingBK1" ) {
      
	if (  maxEdgeLength_ < 1e-15 ) {
	  Error("linElastInt::Calc3DMaterialMat: maxEdgeLength_=0",__FILE__,__LINE__);
	}
	
	Double f1, factor;
	f1 = (minEdgeLength_* minEdgeLength_);
	factor =  2*f1 / ( 2*f1 + maxEdgeLength_ * maxEdgeLength_);
	
	dMat[3][3] *= factor;   
	dMat[4][4] *= factor;   
	dMat[5][5] *= factor;   
      }
      
      else if (softeningPart_ == "bendingSRI" ) {
	for (UInt i=3; i<dimD_; i++) {
	  for (UInt j=3; j<dimD_; j++) {
	    dMat[i][j] = 0.0;
	  }
	}
      }
      
      else if (softeningPart_ == "shearBK1" ) {
	for (UInt i=0; i<3; i++) {
	  for (UInt j=0; j<3; j++) {
	    dMat[i][j] = 0.0;
	  }
	}
	
	Double f1, factor;
	f1 = maxEdgeLength_ * maxEdgeLength_;
	factor = f1 / ( f1 + 2* minEdgeLength_* minEdgeLength_ ); 
	dMat[3][3] *= factor;   
	dMat[4][4] *= factor;   
	dMat[5][5] *= factor; 
      }
      
      else if (softeningPart_ == "shearSRI" ) {
	for (UInt i=0; i<3; i++) {
	  for (UInt j=0; j<3; j++) {
	    dMat[i][j] = 0.0;
	  }
	}
      }
    }
    
  }


  void linElastInt::SetDimensions(SubTensorType type) {

    ENTER_FCN( "linElastInt::SetDimensions" );

    if ( type == FULL ) {
      dimD_   = 6;
      nrDofs_ = 3;
    }
    else if ( type == PLANE_STRAIN || type == PLANE_STRESS ) {
      dimD_   = 3;
      nrDofs_ = 2;
    }
    else if ( type == AXI ) {
      dimD_   = 4;
      nrDofs_ = 2;
      isaxi_  = TRUE;
    }
  }


  // =========================================================================
  // =============== standard con- and destructors (just for tracing) ========
  // =========================================================================

  // NOTE: We hardcode the orientation for 2D simulations to use the yz plane!

  linElastInt::linElastInt( BaseMaterial* matData, SubTensorType type) :
    BDBInt(matData,type) {
    ENTER_FCN( "linElastInt::linElastInt" );

    subTensorType_ = type;
    SetDimensions(type);
  }

  linElastInt::linElastInt( BaseFE *aptelem, BaseMaterial* matData, 
			    SubTensorType type) :
    BDBInt(matData,type) {

    ENTER_FCN( "linElastInt::linElastInt" );

    subTensorType_ = type;
    ptelem=aptelem;
    SetDimensions(type);
  }

  linElastInt::~linElastInt() {
    ENTER_FCN( "linElastInt::~linElastInt" );
  }

} // end namespace CoupledField
