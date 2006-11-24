 #include <iostream>
#include <fstream>
#include <string>

#include <string.h>
#include "basefe.hh"
#include "1D/line1fe.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#ifndef INTEGLIB
#include "olas.hh"
#endif
#include <Utils/tools.hh> 
#include <Matrix/matrix.hh>
#include <DataInOut/WriteInfo.hh> 
 
namespace CoupledField
{

  BaseFE :: BaseFE()
  {
    ENTER_FCN( "entering BaseFE::BaseFE" );
  
    // initializing dynamic objects
    ShFncAtIp_             = NULL;
    ShFncDerivAtIp_        = NULL; 
    ShFncICModesDerivAtIp_ = NULL; 
    IntPoints_             = NULL; 

    // Create dummy ansatzFct object for lagrange functions
    shared_ptr<AnsatzFct> fct( new LagrangeFct() );
    actFct_ = fct;
    actNumFcns_ = 0;

    ICModes_     = false;
    CalcICModes_ = false;
  }
 
  BaseFE :: ~BaseFE()
  {
    ENTER_FCN( "BaseFE::~BaseFE" );
 
    if( ShFncAtIp_ ) delete[] ShFncAtIp_;
    if( ShFncDerivAtIp_ ) delete[] ShFncDerivAtIp_;
    if( IntPoints_ ) delete[] IntPoints_;
    
    // delete our map. The content of normal integration rules are static array
    // the cartesian rules are complete dynamic
    std::map<const std::string, StdVector<Double*>*>::iterator iter; 
    for(iter = IntegrationPointsMap_.begin(); 
        iter != IntegrationPointsMap_.end(); iter++)  {
      StdVector<Double*>* data = iter->second; 
      if(data != NULL) {
        if(iter->first.find("Cartesian", 0 ) != std::string::npos)  { 
          for(UInt i = 0; i < data->GetSize(); i++)  {
            delete (*data)[i];
          } 
        }
        
        delete data;
      }    
    }
    
  }
  
  void BaseFE :: GetShFnc(Vector<double> & S, 
                          const Vector<Double> & LCoord,
                          const Elem* elem,
                          UInt dof  )
  {
    ENTER_FCN( "BaseFE::GetShFnc" );

    CalcShapeFnc(S, LCoord, elem, dof );

  }
  
  void BaseFE :: Local2GlobalCoord(Vector<Double> & globCoord,
                                   const Vector<Double> & locCoord, 
                                   const Matrix<Double> & coordMat,
                                   const Elem* elem ) {
    
    ENTER_FCN( "BaseFE::Local2GlobalCoord");
    
    // step 1: evaluate shape fncs. at local coordinate
    Vector<Double> shFnc;
    CalcShapeFnc(shFnc, locCoord, elem, 1, AnsatzFct::NODE);

    // step2: multiply shape fncs for each dimension with according matrix entries
    globCoord.Resize(Dim_);
    coordMat.Mult(shFnc,globCoord);
    
  }



  void BaseFE :: GetShFncAtIp(Vector<Double> & S, const UInt ip,
                              const Elem * elem, UInt dof )
  {
    ENTER_FCN( "GetShFncAtIp" );
    
    if( actFct_->GetType() == AnsatzFct::LAGRANGE ) {
      S = ShFncAtIp_[ip-1];
    } else {
      CalcShapeFnc( S, IntPoints_[ip-1], elem, dof, 
                    AnsatzFct::ALL );
    }
  }

  void BaseFE :: GetGlobDerivShFnc(Matrix<Double> & Deriv, 
                                   const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords,
                                   const Elem * elem, 
                                   UInt dof )
  {
    ENTER_FCN( "BaseFE::GetGlobDerivShFnc" );

    Deriv.Resize(NumNodes_,Dim_);
    Matrix<Double> LDeriv, JInv;

    CalcLocalDerivShapeFnc(LDeriv, LCoord, elem, dof);

    CalcInvJacobian(JInv, LCoord, CornerCoords, elem );
 
    Deriv = LDeriv * JInv;
  }

  void BaseFE :: GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
                                       const UInt ip,
                                       const Matrix<Double> & CornerCoords,
                                       Double & jacDet,
                                       const Elem * elem, 
                                       UInt dof )
  {
    ENTER_FCN( "BaseFE::GetGlobDerivShFncAtIp" );

    std::string errMsg;

    //  Deriv.Resize(NumNodes_,Dim_);
    Matrix<Double> JInv;
    Double JInvDet;
  
    CalcInvJacobianAtIp(JInv, ip, CornerCoords, elem);

    if( actFct_->GetType() == AnsatzFct::LAGRANGE ) {
      Deriv = ShFncDerivAtIp_[ip-1] * JInv;
    } else {
      Matrix<Double> lDeriv;
      CalcLocalDerivShapeFnc( lDeriv, IntPoints_[ip-1],
                              elem, dof, AnsatzFct::ALL );
      Deriv = lDeriv * JInv;
    }
        
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
                                       const Matrix<Double> & CornerCoords,
                                       const Elem* elem,
                                       UInt dof )
  {
    ENTER_FCN( "BaseFE::GetGlobDerivShFncAtIp" );
  
    //  Deriv.Resize(NumNodes_,Dim_);
    Matrix<Double> JInv;
  
    if ( CalcICModes_ && ICModes_ ) {
      //incompatible modes
      Vector<Double> LCoord( CornerCoords.GetSizeRow());
      LCoord.Init();
      CalcInvJacobian(JInv, LCoord, CornerCoords, elem);
      Deriv = ShFncICModesDerivAtIp_[ip-1] * JInv;
    }
    else {
      CalcInvJacobianAtIp(JInv, ip, CornerCoords, elem);

      if( actFct_->GetType() == AnsatzFct::LAGRANGE ) {
	Deriv = ShFncDerivAtIp_[ip-1] * JInv;
      } 
      else {
	Matrix<Double> lDeriv;
	CalcLocalDerivShapeFnc( lDeriv, IntPoints_[ip-1],
				elem, dof, AnsatzFct::ALL );
	Deriv = lDeriv * JInv;
      }
      //std::cerr << "Deriv = \n" << Deriv << std::endl;
    }
  }




  void BaseFE :: CalcJacobian(Matrix<Double> & J, 
                              const Vector<Double> & LCoord, 
                              const Matrix<Double> & CornerCoords,
                              const Elem* elem )
  {
    ENTER_FCN( "BaseFE::CalcJacobian" );

    //  J.Resize(Dim_,Dim_);

    Matrix<Double> LDeriv;

    CalcLocalDerivShapeFnc(LDeriv, LCoord, elem, AnsatzFct::NODE );
    J = CornerCoords * LDeriv;
  }


  void BaseFE :: CalcJacobianAtIp(Matrix<Double> & J, 
                                  const UInt ip, 
                                  const Matrix<Double> & CornerCoords,
                                  const Elem* elem)
  {
    ENTER_FCN( "BaseFE::CalcJacobianAtIp" );

    if (CornerCoords.GetSizeRow()==3 && Dim_==2) // Surface element in 3D
      J.Resize(CornerCoords.GetSizeRow(),Dim_);

    J = CornerCoords * ShFncDerivAtIp_[ip-1];
  }

  Double BaseFE :: CalcJacobianDet(const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords,
                                   const Elem* elem )
  {
    ENTER_FCN( "BaseFE::CalcJacobianDet" );

    std::string errMsg;
    Matrix<Double> J;
    Double jacDet;

    CalcJacobian( J, LCoord, CornerCoords, elem  );
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
                                       const Matrix<Double> & CornerCoords,
                                       const Elem* elem)
  {
    ENTER_FCN( "BaseFE::CalcJacobianDetAtIp" );

    Matrix<Double> J;
    std::string errMsg;

    CalcJacobianAtIp( J, ip, CornerCoords, elem);

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
                                 const Matrix<Double> & CornerCoords,
                                 const Elem* elem )
  {
    ENTER_FCN( "BaseFE::CalcInvJacobian" );
  
    Matrix<Double> J, LDeriv;

    //  J.Resize(Dim_,Dim_);

    CalcLocalDerivShapeFnc(LDeriv, LCoord, elem, AnsatzFct::NODE);

    J = CornerCoords * LDeriv;

    J.Invert(JInv);
  }



 
  void BaseFE :: CalcInvJacobianAtIp(Matrix<Double> & JInv,
                                     const UInt ip,
                                     const Matrix<Double> & CornerCoords,
                                     const Elem* elem)
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
  
    if (!ShFncAtIp_) {
      ShFncAtIp_ = new Vector<Double>[NumIntPoints_];
    } else{ 
      delete[] ShFncAtIp_ ;
      ShFncAtIp_ = new Vector<Double>[NumIntPoints_];
    }


    for( UInt i=0; i<NumIntPoints_; i++ )
      {
        CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i], 
                      NULL, 1, AnsatzFct::NODE );
      }
  }
  
  void BaseFE :: SetShapeFncDerivAtIp()
  {
    ENTER_FCN( "BaseFE::SetShapeFncDerivAtIp" );

    if( !ShFncDerivAtIp_) {
      ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_]; 
    } 
    else{ 
      delete[] ShFncDerivAtIp_ ;
      ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
    }

    for( UInt i=0; i<NumIntPoints_; i++ )
      CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i], 
                              NULL, 1, AnsatzFct::NODE );

    // check, if incompatible modes are used
    if ( ICModes_ ) {
      if( !ShFncICModesDerivAtIp_ ) {
	 ShFncICModesDerivAtIp_ = new Matrix<Double>[NumIntPoints_]; 
      } 
      else{ 
	delete[]  ShFncICModesDerivAtIp_ ;
	 ShFncICModesDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
      }
      for( UInt i=0; i<NumIntPoints_; i++ ) {
	CalcLocalICModesDerivShapeFnc(  ShFncICModesDerivAtIp_[i], IntPoints_[i], 
					NULL, 1, AnsatzFct::NODE );
      }
    }
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
        ret += GenStr(coordMat[i][j]);
        ret += ", ";
      }
      ret += GenStr(coordMat[coordMat.GetSizeRow()-1][j]);
      ret +=")\n";
    }
    return ret;
  }


  Double BaseFE :: CalcVolume(const Matrix<Double> & CornerCoords, 
                              const bool isaxi)
  {
    ENTER_FCN( "BaseFE::CalcVolume" );

    Double elemVol = 0;
    Double  jacDet, partVol;
    for (UInt actIntPt=1; actIntPt <= NumIntPoints_; actIntPt++) {

      jacDet = CalcJacobianDetAtIp(actIntPt, CornerCoords, NULL);
        
      if (isaxi) {
        Vector<Double> shapeFncAtIp;
        Vector<Double> CoordAtIP;
        GetShFncAtIp(shapeFncAtIp, actIntPt, NULL);
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
    const UInt row = ptCoord.GetSizeRow();
    const UInt col = ptCoord.GetSizeCol();

    TransMat.Resize(spaceDim,true);
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
    for( UInt i = 0; i < row - 1; i++)
      for( UInt j = 0; j < col; j++)
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
    const UInt row = ptCoord.GetSizeRow();
    const UInt col = ptCoord.GetSizeCol();

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

    for( UInt i = row - 1; i >= 0; i-- )
      for( UInt j = col - 1; j >= 0; j-- )
	//transform
	temp[i][j] = ptCoord[i][j] - ptCoord[i][0];

    std::cout << "The new coordinate is\n" << temp << std::endl;

    std::cout << "The base TransMatrix is\n" << TransMat << std::endl;

    NewCoord = TransMat * temp;

    std::cout << "The 2D LocalCoord matrix is\n" << NewCoord << std::endl;

    //rotate
    for( UInt i = 0; i < row - 1; i++)
      for( UInt j = 0; j < col; j++)
	ShellCoord[i][j] = NewCoord[i][j];

  }
  
   
  void BaseFE::MakeKey(IntegrationMethod type, int order, std::string &out)
  {
    ENTER_FCN( "BaseFE::MakeKey" );         
    Enum2String(type, out);       
    std::stringstream ss;
    ss << GetShapeName() << " (" << out << ") order " << order;
    out.assign(ss.str());
  }   

  /** private order encoder. Makes a check and exits on error */
  int BaseFE::EncodeCartesianOrder(int order1, int order2, int order3)
  {
    if(order1 > 99 || order2 > 99 || order3 > 99) 
      Error("Cartesian product numerical integration only with orders < 99", __FILE__, __LINE__ );
        
    return order1 + 100 * order2 + (order3 > 0 && Dim_ == 3 ? 10000 : 0) * order3;
  }

  /** private order decoder */
  void BaseFE::DecodeCartesianOrder(int encoded_order, int* order1, int* order2, int* order3)
  {
    if(encoded_order <= 11 || encoded_order > 999999) 
      Error("Invalid encoded cartesian integration order", __FILE__, __LINE__ );
           
           
    *order3 = (encoded_order >= 10000) ? encoded_order/10000 : 0;   
    encoded_order -= 10000 * (*order3);
        
    *order2 = encoded_order/100;
    encoded_order -= 100 * (*order2);        
        
    *order1 = encoded_order;
        
    // we cannot set it to 0 before as order2 needs to be < 19
    if(Dim_ != 3) *order3 = 0;
  }


  void BaseFE::MakeKey(int order1, int order2, int order3, std::string &out)
  {
    ENTER_FCN( "BaseFE::MakeKey(order1,order2,order3)" );         
    // this is only for debugging, we can encode the orders as defined in environment.hh
    int order = EncodeCartesianOrder(order1, order2, order3);
       
    MakeKey(CARTESIAN, order, out);
  }   



  /** check CreateCartesian... for similar implementation! */   
  void BaseFE::AddIntegrationPoints(IntegrationMethod type, int order, 
                                    int numberOfRows, Double* data)
  {
    ENTER_FCN( "BaseFE::AddIntegrationPoints" );   	   	

    std::string key;
    MakeKey(type, order, key);
        
    // check the key for uniqueness
    if(IntegrationPointsMap_.find(key) != IntegrationPointsMap_.end())
      {
        key.append(" is already in BaseFE::IntegrationPointsMap_");
        Error(key.c_str(), __FILE__, __LINE__ ); 
        exit(-1);         
      }

    // create the payload
    StdVector<Double*>*  value = new StdVector<Double*>(numberOfRows);
    for(int i = 0; i < numberOfRows; i++) (*value)[i]=(data + i * (Dim_+1)); // the data is Dim_ * coordintates + weight

    IntegrationPointsMap_[key] = value;  

    // sum up and print the weights for debugging         
    //Double sum = 0;
    //for(int i = 0; i < numberOfRows; i++) sum += *((*value)[i]+Dim_);
    //std::cout << key << " has weight sum " << sum << std::endl;
       
  }

  StdVector<Double*>* BaseFE::GetIntegrationPoints(IntegrationMethod type, int order, bool search_upwards, 
                                                   bool search_downwards, bool fallback)     
  {
    ENTER_FCN( "BaseFE::GetIntegrationPoints" );   	   	
    std::string key;
    int org_order = order;

    do
      {
        MakeKey(type, order, key);
          
        // store current value so the state is set to what we have
        IntegMethod = type;
        IntegOrder  = order;
   	      
        if(IntegrationPointsMap_.find(key) != IntegrationPointsMap_.end())
          {
            if(order != org_order)
              {
                *warning << "Use Integration " << key << " instead of the wanted order " << org_order << std::endl;
                Warning(__FILE__, __LINE__);
              }    
            return IntegrationPointsMap_[key];
          }   
   	      
        // nothing found
        if(search_upwards) order++;
        if(search_downwards) order--;   
        if(search_upwards && search_downwards) Error("searching up- and downward concurrently!! Stupid!",  __FILE__, __LINE__ );
      }    	
    while(order > 0 && order < 42 && (search_upwards || search_downwards));
    
    // this is a call from the fallback-part
    if(!fallback) return NULL;  

    // print msg;
    MakeKey(type, org_order, key);
    *warning << "No integration points defined for " << key << "! Available are: ";
    std::map<const std::string, StdVector<Double*>*>::iterator iter;
   
    for(iter = IntegrationPointsMap_.begin(); iter != IntegrationPointsMap_.end(); iter++) {
      *warning << iter->first << "\n";
    }                    
    Warning(__FILE__, __LINE__);
           
    // use default
    // restrict to order 2 but stay on 1 if desired
    int alt_o = IntegOrder >= 2 ? 2 : 1;

    // when no economical already use it - if economical doesn't work use classic
    IntegrationMethod alt_m = IntegOrder == ECONOMICAL ? CLASSICAL : ECONOMICAL;
    StdVector<Double*>* data = GetIntegrationPoints(alt_m, alt_o, false, false, false);
      
    // more tries
    if(data == NULL)
      {
        // try again with flipped method
        alt_m = alt_m == ECONOMICAL ? CLASSICAL : ECONOMICAL;          
        data = GetIntegrationPoints(alt_m, alt_o, false, false, false);
      }

    // special pyramid handling :(
    if(data == NULL && alt_o == 1)
      {
        alt_o = 2;
        alt_m = CLASSICAL;
        data = GetIntegrationPoints(alt_m, alt_o, false, false, false);
      } 

    // no more tries   
    if(data != NULL)
      {
        MakeKey(alt_m, alt_o, key);
        *warning << "Use fallback " << key << " for integration\n";
        Warning(__FILE__, __LINE__);
        return data;
      }
          
    Error("No default integration found", __FILE__, __LINE__ );  
    exit(-1);
  } 

  void BaseFE::CommonInit(IntegrationMethod method, int order)
  {
    ENTER_FCN( "BaseFE::CommonInit()" );             
    
#ifndef INTEGLIB
    // if undefined we have an empty constructor and search :)
    if(method == UNDEFINED)
      {  
        // when we habe a XML setting use it
        StdVector<std::string> list;
        // check defaults in case of legacy XML files
        params->GetList("type", list, "integRules");
        if(list.GetSize() > 0)    
          { 
            // we have values in the XML file
            String2Enum(list[0], IntegMethod );    
    
            std::string str;
           
            params->Get( "order", str, "integRules");
            IntegOrder = String2Int(str);    
          }    
        else
          {
            // The default Integration rules
            SetDefaultIntegration();
          }
      }
    else
      {
        IntegMethod = method;
        IntegOrder  = order;
      }
#else
    // The default Integration rules
    SetDefaultIntegration();
#endif        

    // first set integration points and corner coords ...
    // load all into the map
    FillIntegrationPoints();    
    // get the values by IntegMethod and IntegOrder
    SetIntPoints();
    SetCornerCoords();

    // ... then calc shape function values at integration points
    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();      
  }

  void BaseFE::DumpIntegrationPoints()
  {
    std::string str;
    Enum2String(IntegMethod, str);
    
    std::cout << "Current integration points for " << GetShapeName() << " method=" << str << " order=" << IntegOrder << std::endl;

    for(UInt ip=0; ip<NumIntPoints_; ip++)
      {         
        for(UInt dim = 0; dim < Dim_; dim++)
          {
            std::cout << ip << ":" << dim << "=" << IntPoints_[ip][dim] << "  ";
          }
          
        std::cout << "weight=" << IntWeights_[ip] << std::endl;
      }
  }

  void BaseFE::SetCartesianInteg(int order1, int order2, int order3, bool create_only)
  {
    // check if we already have
    std::string key;
    MakeKey(order1, order2, order3, key);

    if(IntegrationPointsMap_.find(key) == IntegrationPointsMap_.end())
      {
        // doesn't exist -> create
        Line1FE* line1 = new Line1FE(ECONOMICAL, order1);
        Line1FE* line2 = new Line1FE(ECONOMICAL, order2);
        Line1FE* line3 = Dim_ == 3 ? new Line1FE(ECONOMICAL, order3) : NULL;
          
        // store in map 
        CreateCartesianIntegration(line1, line2, line3);
          
        // destroy temp objects
        if(line3) delete line3;
        delete line2;
        delete line1;
      }   
      
    // set the state variables so SetIntPoints can load the data from the map
    IntegMethod = CARTESIAN;
    IntegOrder  = EncodeCartesianOrder(order1, order2, order3);
      
    // false only for the call from SetIntPoints()
    if(!create_only) SetIntPoints(); // this calls this method again but the map entry will already the set
  }
  

  /** Creates the integration points via cartesian product for 2D and 3D elements. The result is stored in the map.
   *  @param element1 IntegMethod and IntegOrder shall be set
   *  @param element2 IntegMethod and IntegOrder shall be set
   *  @param element3 the only optional element */
  void BaseFE::CreateCartesianIntegration(BaseFE* line1, BaseFE* line2, BaseFE* line3)
  {
    // we store in map format, as such we have dim corrdinates + weight
    int rows = line1->NumIntPoints_ * line2->NumIntPoints_ * (Dim_ == 3 ? line3->NumIntPoints_ : 1);

    std::string key;
    MakeKey(line1->IntegOrder, line2->IntegOrder, line3 != NULL ? line3->IntegOrder : 0, key);
        
    // check the key for uniqueness
    if(IntegrationPointsMap_.find(key) != IntegrationPointsMap_.end())
      {
        key.append(" is already in BaseFE::IntegrationPointsMap_");
        Error(key.c_str(), __FILE__, __LINE__ ); 
      }
     
    StdVector<Double*>* data = new StdVector<Double*>(rows);
    Double* row = NULL;
     
    int counter = 0;
    for(UInt i = 0; i < line1->NumIntPoints_; i++)
      {
        for(UInt j = 0; j < line2->NumIntPoints_; j++)
          {
            // in 2D run the k loop one time for every i,j iteration
            for(UInt k = 0; k < (Dim_ == 3 ? line3->NumIntPoints_ : 1); k++)
              {
                // create and store our data
                row = new Double[Dim_ + 1];
                (*data)[counter] = row;
                counter++;
                
                // store the coordinates and the weight
                *(row + 0)= line1->IntPoints_[i][0];
                *(row + 1)= line2->IntPoints_[j][0];                
                if(Dim_ == 3) *(row + 2)= line3->IntPoints_[k][0];                
                *(row + Dim_)= line1->IntWeights_[i] * line2->IntWeights_[j] * (Dim_  == 3 ? line3->IntWeights_[k] : 1.);
              }
          }
      }

    IntegrationPointsMap_[key] = data;  
  }

  /** Expicitly set and load the integration type  */
  void BaseFE::SetIntPoints(IntegrationMethod method, int order)
  {
    ENTER_FCN( "BaseFE:SetIntPoints(method,order):" );
    IntegMethod = method;
    IntegOrder  = order;

    SetIntPoints();

    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();  
  }


  
  // reads from the map but generates cartesian integration points on the fly when set in XML for the proper elements!
  void BaseFE::SetIntPoints()
  {
    ENTER_FCN( "BaseFE::SetIntPoints" );   	
      
    // if we are not a valid element for CARTESIAN (assuming reading from XML) the fallback will work
    if(IntegMethod == CARTESIAN 
       && (strcmp(GetShapeName(), "rectangle") == 0  
           || strcmp(GetShapeName(), "hexa") == 0))
      {        
        // create the map entry on the fly
        int order1, order2, order3;
        DecodeCartesianOrder(IntegOrder, &order1, &order2, &order3);
         
        // this is the cached version, true to omit calling this method recursively!
        SetCartesianInteg(order1, order2, order3, true);
      }
      
    // searched upwards, e.g. a quadrilateral for order 2,3 should be stored with higher order 3, so we find for 2
    StdVector<Double*>* data = GetIntegrationPoints(IntegMethod, IntegOrder, true, false, true);
      
    NumIntPoints_= data->GetSize();

    if(IntPoints_ == NULL) {
      IntPoints_ = new Vector<Double>[NumIntPoints_];
    } else {
      delete[] IntPoints_;
      IntPoints_ = new Vector<Double>[NumIntPoints_];
    }
      
    IntWeights_.Resize(NumIntPoints_);
      
    for(UInt ip=0; ip<NumIntPoints_; ip++)
      {         
        IntPoints_[ip].Resize(Dim_);      	
          
        for(UInt dim = 0; dim < Dim_; dim++)
          {
            IntPoints_[ip][dim] = *((*data)[ip]+dim);
          }
          
        IntWeights_[ip]= *((*data)[ip]+Dim_);
      }
      
    // DumpIntegrationPoints();
  }
  
  void BaseFE::DumpIntegrationPointsMap()
  {
    std::string key;
    MakeKey(IntegMethod, IntegOrder, key);
        
    std::cout << "Current Integration method with " << key << " and " << NumIntPoints_ << " integration points\n";
      
    std::map<const std::string, StdVector<Double*>*>::iterator iter;

    for(iter = IntegrationPointsMap_.begin(); iter != IntegrationPointsMap_.end(); iter++) {
      std::cout << iter->first << std::endl;
    }                    
    std::cout << std::flush;
  }  
 
  const char* BaseFE::GetShapeName() 
  {
    // the string representation of the FEType (environment.cc) does not
    // always match the names used on the XML file for integration type :(
    switch(feType()) 
      {
      case  LINE:
        return "line";
      case  TRIA:
        return "triangle";
      case  QUAD:
        return "rectangle";
      case  TET:
        return "tetra";
      case  HEX:
        return "hexa";
      case  PYR:
        return "pyra";
      case  WED:
        return "wedge";
      default:
        Error("type not implemented", __FILE__, __LINE__ );
        exit(-1);
      }
  }
  
  // =======================================================================
  // L E G E N D R E    P A R T
  // =======================================================================
  void  BaseFE::GetNumFncs(Vector<UInt>& numFcns,  
                           const shared_ptr<AnsatzFct>& fcnType, 
                           AnsatzFct::FctEntityType fctEntityType, 
                           UInt dof) {

    ENTER_FCN( "BaseFE::GetNumFcns" );
    
    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      numFcns.Resize( NumNodes_ );
      numFcns.Init(1);
    } else {
      *error << "In base class only implemented for Lagrange functions!";
      Error( __FILE__, __LINE__ );
    }
  }
  
  
  UInt BaseFE::GetNumFncs( const shared_ptr<AnsatzFct>& fcnType ) {

    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      return GetNumNodes();
    } else {
      *error << "In base class only implemented for Lagrange functions!";
      Error( __FILE__, __LINE__ );
    }

    return 0;
    
  }
  
  
  void BaseFE::EvalPolynom( Double& value, Double& deriv,
                            const UInt order, const Double* coeff, 
                            const Double xVal ) {
    ENTER_FCN( "BaseFE::EvalPolynom" );

    // Consider the following expression
    // f(xVal) = a0 * (a1*x^order + a2*x^(order-1) + .. + a(order+1))
    // The coefficients a0..a(order+1) are stored in the coeff-array

    value = coeff[1];
    deriv = 0.0;
    for( UInt i = 2; i < order+2; i++ ) {
      deriv = deriv * xVal + value;
      value = value * xVal + coeff[i];
    }
    // Multiply by pre-factor
    deriv *= coeff[0];
    value *= coeff[0];
  }
  
  
  // Define coefficients for legendre ansatz functions up to order 8
  Double  BaseFE::lCoeff_[9][10] = {
    {0.5                  ,   -1, 1,    0, 0,   0, 0,    0, 0, 0 },
    {0.5                  ,    1, 1,    0, 0,   0, 0,    0, 0, 0 },
    {0.25*sqrt(6.0)       ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
    {0.25*sqrt(10.0)      ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
    {1.0/16.0*sqrt(14.0)  ,    5, 0,   -6, 0,   1, 0,    0, 0, 0 },
    {3.0/16.0*sqrt(2.0)   ,    7, 0,  -10, 0,   3, 0,    0, 0, 0 },
    {1.0/32.0*sqrt(22.0)  ,   21, 0,  -35, 0,  15, 0,   -1, 0, 0 },
    {1.0/32.0*sqrt(26.0)  ,   33, 0,  -63, 0,  35, 0,   -5, 0, 0 },
    {1.0/256.0*sqrt(30.0) ,  429, 0, -924, 0, 630, 0, -140, 0, 5 }
  };
  
} // end namespace CoupledField
