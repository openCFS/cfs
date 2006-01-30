#include <iostream>
#include <fstream>
#include <string>

#include "basefe.hh"
#include "olas.hh"
#include <Utils/tools.hh> 
#include <Matrix/matrix.hh>
#include <DataInOut/WriteInfo.hh> 
 
namespace CoupledField
{

  BaseFE :: BaseFE()
  {
    ENTER_FCN( "entering BaseFE::BaseFE" );
  
    // initializing dynamic objects
    ShFncAtIp_      = NULL;
    ShFncDerivAtIp_ = NULL; 
    IntPoints_      = NULL; 
  
  }
 
  BaseFE :: ~BaseFE()
  {
    ENTER_FCN( "BaseFE::~BaseFE" );
 
    if( ShFncAtIp_ ) delete[] ShFncAtIp_;
    if( ShFncDerivAtIp_ ) delete[] ShFncDerivAtIp_;
    if( IntPoints_ ) delete[] IntPoints_;
  }

  void BaseFE :: GetShFnc(Vector<double> & S, 
                          const Vector<Double> & LCoord)
  {
    ENTER_FCN( "BaseFE::GetShFnc" );

    CalcShapeFnc(S, LCoord);

  }
  
  void BaseFE :: Local2GlobalCoord(Vector<Double> & globCoord,
                                   const Vector<Double> & locCoord, 
                                   const Matrix<Double> & coordMat ) {
    
    ENTER_FCN( "BaseFE::Local2GlobalCoord");
    
    // step 1: evaluate shape fncs. at local coordinate
    Vector<Double> shFnc;
    CalcShapeFnc(shFnc, locCoord);

    // step2: multiply shape fncs for each dimension with according matrix entries
    globCoord.Resize(Dim_);
    coordMat.Mult(shFnc,globCoord);
    
  }



  void BaseFE :: GetShFncAtIp(Vector<Double> & S, 
                              const UInt ip)
  {
    ENTER_FCN( "GetShFncAtIp" );

    S.Resize(NumNodes_);

    S = ShFncAtIp_[ip-1];
  
  }

  void BaseFE :: GetGlobDerivShFnc(Matrix<Double> & Deriv, 
                                   const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::GetGlobDerivShFnc" );

    Deriv.Resize(NumNodes_,Dim_);
    Matrix<Double> LDeriv, JInv;

    CalcLocalDerivShapeFnc(LDeriv, LCoord);

    CalcInvJacobian(JInv, LCoord, CornerCoords);
 
    Deriv = LDeriv * JInv;
  }

  void BaseFE :: GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
                                       const UInt ip,
                                       const Matrix<Double> & CornerCoords,
                                       Double & jacDet)
  {
    ENTER_FCN( "BaseFE::GetGlobDerivShFncAtIp" );

    std::string errMsg;

    //  Deriv.Resize(NumNodes_,Dim_);
    Matrix<Double> JInv;
    Double JInvDet;
  
    CalcInvJacobianAtIp(JInv, ip, CornerCoords);

    Deriv = ShFncDerivAtIp_[ip-1] * JInv;

    // det(A) = 1 / det(A^(-1))
    JInv.Determinant(JInvDet);
    jacDet = 1.0 / JInvDet;


    if ( jacDet < 0.0 ){
      errMsg = "BaseFE:GetGlobDerivShFncAtIp: ";
      errMsg += "Negative Jacobian Determinant!\n";
      errMsg += "The coordinates of the element are:\n";
      errMsg += CoordMatrix2String(CornerCoords);
      Error(errMsg.c_str(), __FILE__, __LINE__ );     
    }

 
  }


  void BaseFE :: GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
                                       const UInt ip,
                                       const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::GetGlobDerivShFncAtIp" );
  
    //  Deriv.Resize(NumNodes_,Dim_);
    Matrix<Double> JInv;
  
    CalcInvJacobianAtIp(JInv, ip, CornerCoords);
  
    Deriv = ShFncDerivAtIp_[ip-1] * JInv;
  }




  void BaseFE :: CalcJacobian(Matrix<Double> & J, 
                              const Vector<Double> & LCoord, 
                              const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::CalcJacobian" );

    //  J.Resize(Dim_,Dim_);

    Matrix<Double> LDeriv;

    CalcLocalDerivShapeFnc(LDeriv, LCoord);
    J = CornerCoords * LDeriv;
  }


  void BaseFE :: CalcJacobianAtIp(Matrix<Double> & J, 
                                  const UInt ip, 
                                  const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::CalcJacobianAtIp" );

    if (CornerCoords.GetSizeRow()==3 && Dim_==2) // Surface element in 3D
      J.Resize(CornerCoords.GetSizeRow(),Dim_);

    J = CornerCoords * ShFncDerivAtIp_[ip-1];
  }

  Double BaseFE :: CalcJacobianDet(const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::CalcJacobianDet" );

    std::string errMsg;
    Matrix<Double> J;
    Double jacDet;

    CalcJacobian( J, LCoord, CornerCoords );
    J.Determinant(jacDet);

    if ( jacDet < 0.0 ){
      errMsg = "BaseFE:CalcJacobianDet: Negative Jacobian Determinant!\n";
      errMsg += "The coordinates of the element are:\n";
      errMsg += CoordMatrix2String(CornerCoords);
      Error(errMsg.c_str(), __FILE__, __LINE__ );     
    }

    return jacDet;
  }


  Double BaseFE :: CalcJacobianDetAtIp(const UInt ip, 
                                       const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::CalcJacobianDetAtIp" );

    Matrix<Double> J;
    std::string errMsg;

    CalcJacobianAtIp( J, ip, CornerCoords);

    if (CornerCoords.GetSizeRow()==3 && Dim_==2)
      {
        Vector<Double> normal;
        normal.Resize(CornerCoords.GetSizeRow());  
        normal[0]= J[1][0]* J[2][1]- J[2][0]* J[1][1];
        normal[1]=J[2][0]*J[0][1]- J[0][0]* J[2][1];  
        normal[2]= J[0][0]* J[1][1]- J[1][0]*J[0][1];

        Double detJ = sqrt(sqr(normal[0])+sqr(normal[1])+sqr(normal[2]));

        if ( detJ < 0.0 ){
          errMsg = "BaseFE:CalcJacobianDet: Negative Jacobian Determinant!\n";
          errMsg += "The coordinates of the element are:\n";
          errMsg += CoordMatrix2String(CornerCoords);
          Error(errMsg.c_str(), __FILE__, __LINE__ );     
        }
        return detJ;
      }

    else  {
      Double jacDet ;
      J.Determinant(jacDet);
      if ( jacDet < 0.0 ){
        errMsg = "BaseFE:CalcJacobianDetAtIp: Negative Jacobian Determinant!\n";
        errMsg += "The coordinates of the element are:\n";
        errMsg += CoordMatrix2String(CornerCoords);
        Error(errMsg.c_str(), __FILE__, __LINE__ );     
      }

      return jacDet;
    }  
  }




  void BaseFE :: CalcInvJacobian(Matrix<Double> & JInv,
                                 const Vector<Double> & LCoord,
                                 const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::CalcInvJacobian" );
  
    Matrix<Double> J, LDeriv;

    //  J.Resize(Dim_,Dim_);

    CalcLocalDerivShapeFnc(LDeriv, LCoord);

    J = CornerCoords * LDeriv;

    J.Invert(JInv);
  }



 
  void BaseFE :: CalcInvJacobianAtIp(Matrix<Double> & JInv,
                                     const UInt ip,
                                     const Matrix<Double> & CornerCoords)
  {
    ENTER_FCN( "BaseFE::CalcInvJacobianAtIp" );

    JInv.Resize(Dim_,Dim_);

    Matrix<Double> J;

    //  J.Resize(Dim_,Dim_);

    J = CornerCoords * ShFncDerivAtIp_[ip-1];

    J.Invert(JInv);
  }


  void BaseFE :: SetShapeFncAtIp()
  {
    ENTER_FCN( "BaseFE::SetShapeFncAtIp" );
  
    if (!ShFncAtIp_)
      ShFncAtIp_ = new Vector<Double>[NumIntPoints_];


    for( UInt i=0; i<NumIntPoints_; i++ )      
      CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);
  }
  
  void BaseFE :: SetShapeFncDerivAtIp()
  {
    ENTER_FCN( "BaseFE::SetShapeFncDerivAtIp" );

    if( !ShFncDerivAtIp_)
      ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];

    for( UInt i=0; i<NumIntPoints_; i++ )
      CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);

  }



  enum IntegrationType BaseFE::String2EnumIntegrationType(const Char * inttype)
  {
    ENTER_FCN( "BaseFE::String2EnumIntegrationType" );

    enum IntegrationType result=null;
    if (!strcmp(inttype,"GaussOrder2")) result=GaussOrder2;
    if (!strcmp(inttype,"GaussOrder3")) result=GaussOrder3;
    if (!strcmp(inttype,"GaussOrder4")) result=GaussOrder4;
    if (!strcmp(inttype,"GaussOrder5")) result=GaussOrder5;
    if (!strcmp(inttype,"GaussOrder7")) result=GaussOrder7;

    if (result==null) 
      Error("Check your config file. Your integration type is wrong",
            __FILE__,__LINE__);

    return result;
  }



  void BaseFE::GetGlobalEdgeIndicesAtIP( Vector<Double> & globCoord, UInt ip,
                                         const Matrix<Double> & cornerCoords)
  {
    ENTER_FCN( "BaseFE::GetGlobalEdgeIndicesAtIP" );
    // cornerCoords: nrCorners x dim

    globCoord.Resize(Dim_);
  
    globCoord =  cornerCoords * ShFncAtIp_[ip-1];
  }




  // calculate global derivates of edge shape functions in integration point ip
  void BaseFE::GetEdgeGlobDerivShFncAtIp(StdVector< Matrix<Double>* > & deriv, 
                                         const UInt ip,
                                         const Matrix<Double> & cornerCoords)
  {
    ENTER_FCN( "BaseFE::GetEdgeGlobDerivShFncAtIp" );

    // vector of coordinates of the desired integration point
    Vector<Double> lCoord;
    lCoord.Resize(Dim_);

    for (UInt i=0; i<Dim_; i++)
      lCoord[i] = IntPoints_[ip-1][i];
  
    GetEdgeGlobalDerivShapeFnc(deriv, lCoord, cornerCoords);
  }

  void BaseFE::CalcEdgeShapeFncAtIp(Matrix<Double> & shape, 
                                    const UInt ip,
                                    const Matrix<Double> & cornerCoords)
  {
    ENTER_FCN( "BaseFE::CalcEdgeShapeFncAtIp" );

    Vector<Double> lCoord;
    lCoord.Resize(Dim_);

    for (UInt i=0; i<Dim_; i++)
      lCoord[i] = IntPoints_[ip-1][i];
  
    CalcEdgeShapeFnc(shape, lCoord, cornerCoords);
  }


  void BaseFE::GetGlobalEdgeIndices(StdVector<UInt> & globEdgeIndex,
                                    UInt * pDENodes, 
                                    BaseSystem * algsys)
  {
    ENTER_FCN( "BaseFE::GetGlobalEdgeIndices" );

    Error( "Edge functions are currently not supported!", __FILE__, __LINE__ );

    // define the global edge number
    //for(UInt actEdge=0; actEdge < globEdgeIndex.GetSize(); actEdge++)
    //  globEdgeIndex[actEdge] = 
    //    algsys->GetNode2Edge(pDENodes[ edgeVertices_[actEdge][0]],
    // pDENodes[ edgeVertices_[actEdge][1]]);
  }

  std::string BaseFE::CoordMatrix2String(const Matrix<Double> & coordMat)
  {
    std::string ret;
    for (UInt j=0; j<coordMat.GetSizeCol(); j++) {
      ret += "(";
      for (UInt i=0; i<coordMat.GetSizeRow()-1; i++) {
        ret += Info->GenStr(coordMat[i][j]);
        ret += ", ";
      }
      ret += Info->GenStr(coordMat[coordMat.GetSizeRow()-1][j]);
      ret +=")\n";
    }
    return ret;
  }

  void BaseFE::SetReducedIntegration()
  {
    ENTER_FCN( "BaseFE:SetReducedIntegration:" );

    IntegType =  GaussOrder1;

    SetIntPoints();
    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();  

  }

  void BaseFE::SetStandardIntegration()
  {
    ENTER_FCN( "BaseFE:SetStandardIntegration:" );

    IntegType =  GaussOrder2;

    SetIntPoints();
    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();  

  }

  Double BaseFE :: CalcVolume(const Matrix<Double> & CornerCoords, 
                              const Boolean isaxi)
  {
    ENTER_FCN( "BaseFE::CalcVolume" );

    Double elemVol = 0;
    Double  jacDet, partVol;
    for (UInt actIntPt=1; actIntPt <= NumIntPoints_; actIntPt++) {

      jacDet = CalcJacobianDetAtIp(actIntPt, CornerCoords);
        
      if (isaxi) {
        Vector<Double> shapeFncAtIp;
        Vector<Double> CoordAtIP;
        GetShFncAtIp(shapeFncAtIp, actIntPt);
        CoordAtIP = CornerCoords * shapeFncAtIp;
        partVol = 2 * PI * IntWeights_[actIntPt-1] * jacDet * CoordAtIP[0];
      }
      else 
        partVol = IntWeights_[actIntPt-1] * jacDet;
    
      elemVol += partVol;
    }
  
    return elemVol;
  }

  /////////////////////////////////////////////////////////////////////
  /////// For geometrical transformation using direction cosines //////

  void BaseFE::CoordTrans( const Matrix<Double> &ptCoord, 
			     Matrix<Double> &TransMat, 
			     Matrix<Double> &ShellCoord )
  {
    ENTER_FCN( "BaseFE::CoordTrans" );

    // Method based on book XXX and Master thesis from Xiong //
    const Integer spaceDim = 3;
    std::vector<Double> Vx, Vy, Vz;
    const Integer row = ptCoord.GetSizeRow();
    const Integer col = ptCoord.GetSizeCol();

    TransMat.Resize(spaceDim);
    TransMat.Init();

    Matrix<Double> NewCoord, temp;
    NewCoord.Resize( row, col );
    NewCoord.Init();
    temp.Resize( row, col );
    temp.Init();

    ShellCoord.Resize( row - 1, col );
    ShellCoord.Init();

    Double length;

    Vx.push_back( ptCoord[0][1] - ptCoord[0][0] );
    Vx.push_back( ptCoord[1][1] - ptCoord[1][0] );
    Vx.push_back( ptCoord[2][1] - ptCoord[2][0] );

    length = 1.0 / sqrt( Vx[0] * Vx[0] + Vx[1] * Vx[1] + Vx[2] * Vx[2] );
 
    //   cos(x_,x), cos(x_,y), cos(x_,z)
    TransMat[0][0] = Vx[0] * length;
    TransMat[0][1] = Vx[1] * length;
    TransMat[0][2] = Vx[2] * length;

    Vy.push_back( ptCoord[0][2] - ptCoord[0][0] );
    Vy.push_back( ptCoord[1][2] - ptCoord[1][0] );
    Vy.push_back( ptCoord[2][2] - ptCoord[2][0] );

    Vz.push_back( Vx[1] * Vy[2] - Vx[2] * Vy[1] );
    Vz.push_back( Vx[2] * Vy[0] - Vx[0] * Vy[2] );
    Vz.push_back( Vx[0] * Vy[1] - Vx[1] * Vy[0] );

    length = 1.0 / sqrt( Vz[0] * Vz[0] + Vz[1] * Vz[1] + Vz[2] * Vz[2] );

    TransMat[2][0] = Vz[0] * length;
    TransMat[2][1] = Vz[1] * length;
    TransMat[2][2] = Vz[2] * length;

    TransMat[1][0] = TransMat[0][2] * TransMat[2][1] - TransMat[0][1] 
      * TransMat[2][2];
    TransMat[1][1] = TransMat[0][0] * TransMat[2][2] - TransMat[0][2] 
      * TransMat[2][0];
    TransMat[1][2] = TransMat[0][1] * TransMat[2][0] - TransMat[0][0] 
      * TransMat[2][1];

    // transform geometry from real(global) coordinate to standard(local) 
    // coordinate

    for( int i = row - 1; i >= 0; i-- )
      for( int j = col - 1; j >= 0; j-- )
	//transform
	temp[i][j] = ptCoord[i][j] - ptCoord[i][0];

//     std::cout << "The new coordinate is\n" << temp << std::endl;

//     std::cout << "The base TransMatrix is\n" << TransMat << std::endl;

    NewCoord = TransMat * temp;

//     std::cout << "The 3D LocalCoord matrix is\n" << NewCoord << std::endl;

    //rotate
    for( int i = 0; i < row - 1; i++)
      for( int j = 0; j < col; j++)
	ShellCoord[i][j] = NewCoord[i][j];

  }

  void BaseFE::CoordTrans2D( const Matrix<Double> &ptCoord, 
			     Matrix<Double> &TransMat, 
			     Matrix<Double> &ShellCoord )
  {
    ENTER_FCN( "BaseFE::CoordTrans2D" );

    // Method based on book XXX and Master thesis from Xiong //
    const Integer spaceDim = 2;
    std::vector<Double> Vx, Vy;
    const Integer row = ptCoord.GetSizeRow();
    const Integer col = ptCoord.GetSizeCol();

    TransMat.Resize(spaceDim);
    TransMat.Init();

    Matrix<Double> NewCoord, temp;
    NewCoord.Resize( row, col );
    NewCoord.Init();
    temp.Resize( row, col );
    temp.Init();

    ShellCoord.Resize( row - 1, col );
    ShellCoord.Init();

    Double length;

    Vx.push_back( ptCoord[0][1] - ptCoord[0][0] );
    Vx.push_back( ptCoord[1][1] - ptCoord[1][0] );
 

    length = 1.0 / sqrt( Vx[0] * Vx[0] + Vx[1] * Vx[1] );
 
    //   cos(x_,x), cos(x_,y), cos(y_,x), cos(y_,y)
    TransMat[0][0] = Vx[0] * length;
    TransMat[0][1] = Vx[1] * length;
    TransMat[1][0] = Vx[1] * length;
    TransMat[1][1] = Vx[0] * length;


    // transform geometry from real(global) coordinate to standard(local) 
    // coordinate

    for( int i = row - 1; i >= 0; i-- )
      for( int j = col - 1; j >= 0; j-- )
	//transform
	temp[i][j] = ptCoord[i][j] - ptCoord[i][0];

     std::cout << "The new coordinate is\n" << temp << std::endl;

     std::cout << "The base TransMatrix is\n" << TransMat << std::endl;

    NewCoord = TransMat * temp;

     std::cout << "The 2D LocalCoord matrix is\n" << NewCoord << std::endl;

    //rotate
    for( int i = 0; i < row - 1; i++)
      for( int j = 0; j < col; j++)
	ShellCoord[i][j] = NewCoord[i][j];

  }

} // end namespace CoupledField
