#include <iostream>
#include <fstream>

#include "nrbcInt.hh"

namespace CoupledField {

  NrbcInt::NrbcInt(const Double aDensity,  const UInt nrDofsPerNode, bool axi)
    : BaseForm( NULL ), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode)
  {
    ENTER_FCN( "NrbcInt::NrbcInt" );

    name_ = "NrbcInt";

    factor_ = 1.0;
    isaxi_ = axi;
    nrbcMatType_ = 0;
    baseType_ = MASS;
  }

  NrbcInt::NrbcInt(const Double aDensity,  const UInt nrDofsPerNode, UInt nrbcMatType,
                   bool axi)
    : BaseForm( NULL ), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode),
      nrbcMatType_(nrbcMatType)
  {
    ENTER_FCN( "NrbcInt::NrbcInt" );

    name_ = "NrbcInt";
    factor_ = 1.0;
    isaxi_ = axi;

    if ((nrbcMatType_ == 1)||(nrbcMatType_ == 3))
      baseType_ = STIFFNESS;
    else if (nrbcMatType_ == 2)
      baseType_ = MASS;      
      
  }
 
  NrbcInt::~NrbcInt()
  {
    ENTER_FCN( "NrbcInt::~NrbcInt" );
  }

  void NrbcInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1, 
                                   EntityIterator& ent2 )
  {
    ENTER_FCN( "NrbcInt::CalcElemMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(numFncs);
    elemMat.Resize(numFncs);
    elemMat.Init();

     Matrix<Double> TransMattotest;    
//     //to test coord transformation in 3d surf
//      Matrix<Double>  ptCoordtotest;
//     Matrix<Double> ShellCoord;

//     ptCoordtotest.Resize(3,3);
//     ptCoordtotest[0][0]=2;
//     ptCoordtotest[1][0]=0;
//     ptCoordtotest[2][0]=2;

//     ptCoordtotest[0][1]=3;
//     ptCoordtotest[1][1]=0;
//     ptCoordtotest[2][1]=1;

//     ptCoordtotest[0][2]=2;
//     ptCoordtotest[1][2]=1;
//     ptCoordtotest[2][2]=0;
//     ptelem->CoordTrans( ptCoordtotest, TransMattotest, ShellCoord );
//     std::cout<<"ShellCoord= "<<std::endl;
//     std::cout<<ShellCoord<<std::endl;
//     std::cout<<"TransMattotest= "<<std::endl;
//     std::cout<<TransMattotest<<std::endl;
//     //end of coordTransf test
    
    // we can have two type of integrators, (NaNb) or (Na_tangDerivNb_tangDeriv)

    if (nrbcMatType_ == 3)//for integrator with tangential derivatives
      {
        Matrix<Double> xiDx;
        Matrix<Double> xiDxTransp;
        Vector<Double> normal;
        for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
          
          jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
          
          //ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);
          if(ptCoord_.GetSizeRow()==3)
            {
              Matrix<Double> ptCoord_onXY;
              Matrix<Double> TransMat;
              ptelem->CoordTrans( ptCoord_, TransMattotest, ptCoord_onXY );
              ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_onXY, 
                                            jacDet, ent1.GetElem() );
            }
          else
            {
              Matrix<Double> ptCoord_onXline;              
               Double length=dist_Mat(ptCoord_);
                //sqrt(pow((ptCoord[0][1]-ptCoord[0][0]),2)+pow((ptCoord[1][1]-ptCoord[1][0]),2));
               ptCoord_onXline.Resize(1,2);
               ptCoord_onXline.Init();
              ptCoord_onXline[0][0]=0.;
              ptCoord_onXline[0][1]=length;
              //              ptCoord_onXline[1][0]=0.;
              //ptCoord_onXline[1][1]=0.;
              ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, 
                                            ptCoord_onXline, ent1.GetElem());
 
              //ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord);

//                std::cout<<"Rows= "<<xiDx.GetSizeRow()<<std::endl;
//                std::cout<<"Columns= "<<xiDx.GetSizeCol()<<std::endl;
//                std::cout<<"xiDx= "<<std::endl;
//                std::cout<<xiDx<<std::endl;

               jacDet=length/2;
//                std::cout<<"Length: "<<length<<std::endl;

//               Matrix<Double> TransMat;
//               ptelem->CoordTrans2D( ptCoord, TransMattotest, ptCoord_onXline );
//               std::cout<<"ptCoord:\n "<<ptCoord<<" ptCoord_onXline:\n "<<ptCoord_onXline<<std::endl;
//               ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_onXline);



            }
          
          xiDx.Transpose(xiDxTransp);
          partElemMat = xiDx * xiDxTransp;         
          
          if (isaxi_) {
            Warning( "NRBC Q integrator is not implemented for axisymmetric mode", __FILE__, __LINE__ );
            //             CoordAtIP = ptCoord * shapeFncAtIp;
            //             partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * factor_* jacDet * CoordAtIP[0];
          }
          else 
            partElemMat *= intWeights[actIntPt-1] * density_ * factor_ * jacDet;
          //           std::cout<<"partElemMat= "<<std::endl;
          //           std::cout<<partElemMat<<std::endl;
          elemMat += partElemMat;
        }

        if(ptCoord_.GetSizeRow()==3)
          {
            LocaltoGlob( elemMat, TransMattotest  );
          }
 //        std::cout<<"elemMat= "<<std::endl;
//         std::cout<<elemMat<<std::endl;
      }
    
    else
      {
        for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
          
          jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
          
          ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
          
          partElemMat.DyadicMult(shapeFncAtIp);
          
          if (isaxi_) {
            CoordAtIP = ptCoord_ * shapeFncAtIp;
            partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * factor_* jacDet * CoordAtIP[0];
          }
          else 
            partElemMat *= intWeights[actIntPt-1] * density_ * factor_ * jacDet;
          
          elemMat += partElemMat;        
        }
      }
    
    if (nrDofsPerNode_ > 1 && nrbcMatType_ == 0 ) {
      Matrix <Double> multDofNRBC;
      NRBCMultiDof(multDofNRBC, elemMat, nrDofsPerNode_);       
      elemMat = multDofNRBC;
    }
    else if (nrDofsPerNode_ > 1 && nrbcMatType_ == 1 )
      {// 1 corresponds to R element matrix
      Matrix <Double> multDofNRBC;
      NRBCMultiDof_Rmat(multDofNRBC, elemMat);       
      elemMat = multDofNRBC;
      }
    else if (nrDofsPerNode_ > 1 && ((nrbcMatType_ == 2)||(nrbcMatType_ == 3)) )
      {// 2 corresponds to both P or Q element matrix
      Matrix <Double> multDofNRBC;
      NRBCMultiDof_PorQmat(multDofNRBC, elemMat);       
      elemMat = multDofNRBC;
      }

  }

  void NrbcInt::NRBCMultiDof(Matrix<Double>& nrbcMultDof, 
                             const Matrix<Double>& nrbcMatSingleDof,  const UInt nrDofs)
  {
    ENTER_FCN( "NrbcInt::NRBCMultiDof" );
    
    const UInt singleDofSize = nrbcMatSingleDof.GetSizeRow();

    UInt i, j, actDof;
    
    nrbcMultDof.Resize(nrDofsPerNode_*singleDofSize);
    nrbcMultDof.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
        for (actDof=0; actDof < nrDofs; actDof++)
          nrbcMultDof[i*nrDofs + actDof][j*nrDofs + actDof] = nrbcMatSingleDof[i][j]; 
  }

  void NrbcInt::NRBCMultiDof_Rmat(Matrix<Double>&nrbcMultDof_Rmat, const Matrix<Double>& nrbcMatSingleDof)
  {
    ENTER_FCN( "NrbcInt::NRBCMultiDofZero" );
    
    const UInt singleDofSize = nrbcMatSingleDof.GetSizeRow();
    
    UInt i, j, actDof;
    
    nrbcMultDof_Rmat.Resize(nrDofsPerNode_*singleDofSize);
    nrbcMultDof_Rmat.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
        for (actDof=0; actDof < (nrDofsPerNode_-1) ; actDof++)
          nrbcMultDof_Rmat[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof + 1] = 
                nrbcMatSingleDof[i][j]; 
  }

  void NrbcInt::NRBCMultiDof_PorQmat(Matrix<Double>& nrbcMultDof_PorQmat, const Matrix<Double>& nrbcMatSingleDof)
  {
    ENTER_FCN( "NrbcInt::NRBCMultiDofZero" );
    
    const UInt singleDofSize = nrbcMatSingleDof.GetSizeRow();
    
    UInt i, j, actDof;
    
    nrbcMultDof_PorQmat.Resize(nrDofsPerNode_*singleDofSize);
    nrbcMultDof_PorQmat.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
        for (actDof=0; actDof < (nrDofsPerNode_-1) ; actDof++)
          nrbcMultDof_PorQmat[i*nrDofsPerNode_ + actDof + 1][j*nrDofsPerNode_ + actDof] = 
                nrbcMatSingleDof[i][j]; 
  }

  void NrbcInt::LocaltoGlob( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat )
  {

    ENTER_FCN( "NrbcInt::LocaltoGlob" );

    int n = 0, i, j;
    const Integer row = ElemMat.GetSizeRow();
    const Integer spaceDim = 3;

    Matrix<Double> StiffTrans;
    StiffTrans.Resize( row );
    StiffTrans.Init();

    Matrix<Double> TFTrans;
    Matrix<Double> KTF;
    while( n < row / spaceDim )
      {
	for( i = 0; i < spaceDim; i++ )
	  for( j = 0; j < spaceDim; j++ )
	    {
	      StiffTrans[ n * spaceDim + i][ n * spaceDim + j] = TransMat[i][j];
	    }
	n++;
      }

#ifdef DEBUG
     (*debug) << std::endl << "The base transform is\n" << TransMat << std::endl;
     (*debug) << std::endl << "the transformation matrix is " << StiffTrans.GetSizeRow() <<
       " x " << StiffTrans.GetSizeCol() << std::endl;
     (*debug) << StiffTrans << std::endl;

     (*debug) << std::endl<< "before transform the element matrix is "<< ElemMat.GetSizeRow()<<
       " x " << ElemMat.GetSizeCol() << std::endl;
     (*debug) << ElemMat << std::endl;
#endif

    KTF = ElemMat * StiffTrans;

    StiffTrans.Transpose(TFTrans);

    ElemMat = TFTrans * KTF;

    // to check if the transform matrix is orthogonal
    KTF = TFTrans * StiffTrans;

#ifdef DEBUG
     (*debug) << std::endl<< "after transform the element matrix is "<< ElemMat.GetSizeRow()<<
       " x " << ElemMat.GetSizeCol() << std::endl;
     (*debug) << ElemMat << std::endl;
     (*debug) << std::endl << "KTF is\n" << KTF << std::endl;
#endif
  }















} // end namespace CoupledField



