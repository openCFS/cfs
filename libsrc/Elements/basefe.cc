#include <iostream>
#include <fstream>
#include <string>

#include "basefe.hh"
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

void BaseFE :: GetShFncAtIp(Vector<Double> & S, 
			    const Integer ip)
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
				     const Integer ip,
				     const Matrix<Double> & CornerCoords,
				     Double & jacDet)
{
  ENTER_FCN( "BaseFE::GetGlobDerivShFncAtIp" );

  //  Deriv.Resize(NumNodes_,Dim_);
  Matrix<Double> JInv;
  Double JInvDet;
  
  CalcInvJacobianAtIp(JInv, ip, CornerCoords);

  Deriv = ShFncDerivAtIp_[ip-1] * JInv;

  // det(A) = 1 / det(A^(-1))
  JInv.Determinant(JInvDet);
  jacDet = 1.0 / JInvDet;

  if ( jacDet < 0.0 ) {
    std::string msg = "Coordinates: ";
    Info->PrintMatrix(msg, CornerCoords);
    std::cout << "Jdet = " << jacDet << std::endl;
    Error( "Negative Jacobian determinante ", __FILE__, __LINE__ );
  }

}


void BaseFE :: GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				     const Integer ip,
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
				const Integer ip, 
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

  Matrix<Double> J;
  Double jacDet;

  CalcJacobian( J, LCoord, CornerCoords );
  J.Determinant(jacDet);
  if ( jacDet < 0.0 )
    Error( "Negative Jacobian determinante ", __FILE__, __LINE__ );
  return jacDet;
}

Double BaseFE :: CalcJacobianDetAtIp(const Integer ip, 
				     const Matrix<Double> & CornerCoords)
{
  ENTER_FCN( "BaseFE::CalcJacobianDetAtIp" );

  Matrix<Double> J;

  CalcJacobianAtIp( J, ip, CornerCoords);

  if (CornerCoords.GetSizeRow()==3 && Dim_==2)
    {
      Vector<Double> normal;
      normal.Resize(CornerCoords.GetSizeRow());  
      normal[0]= J[1][0]* J[2][1]- J[2][0]* J[1][1];
      normal[1]=J[2][0]*J[0][1]- J[0][0]* J[2][1];  
      normal[2]= J[0][0]* J[1][1]- J[1][0]*J[0][1];

      Double detJ = sqrt(sqr(normal[0])+sqr(normal[1])+sqr(normal[2]));
      if ( detJ < 0.0 )
	Error( "Negative Jacobian determinante ", __FILE__, __LINE__ );

      return detJ;
    }

  else
    {
      Double jacDet ;
      J.Determinant(jacDet);
      if ( jacDet < 0.0 )
	Error( "Negative Jacobian determinante ", __FILE__, __LINE__ );     
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
				   const Integer ip,
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


  for( Integer i=0; i<NumIntPoints_; i++ )      
      CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);
}
  
void BaseFE :: SetShapeFncDerivAtIp()
{
  ENTER_FCN( "BaseFE::SetShapeFncDerivAtIp" );

  if( !ShFncDerivAtIp_)
    ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];

  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);

}



enum IntegrationType BaseFE::String2EnumIntegrationType(const Char * inttype)
{
  ENTER_FCN( "BaseFE::String2EnumIntegrationType" );

 enum IntegrationType result=null;
 if (!std::strcmp(inttype,"GaussOrder2")) result=GaussOrder2;
 if (!std::strcmp(inttype,"GaussOrder3")) result=GaussOrder3;
 if (!std::strcmp(inttype,"GaussOrder4")) result=GaussOrder4;
 if (!std::strcmp(inttype,"GaussOrder5")) result=GaussOrder5;
 if (!std::strcmp(inttype,"GaussOrder7")) result=GaussOrder7;

 if (result==null) Error("Check your config file. Your integration type is wrong",__FILE__,__LINE__);

 return result;
}



void BaseFE::GetGlobalEdgeIndicesAtIP( Vector<Double> & globCoord, Integer ip,
				       const Matrix<Double> & cornerCoords)
{
  ENTER_FCN( "BaseFE::GetGlobalEdgeIndicesAtIP" );
  // cornerCoords: nrCorners x dim

  globCoord.Resize(Dim_);
  
  globCoord =  cornerCoords * ShFncAtIp_[ip-1];
}




// calculate global derivates of edge shape functions in integration point ip
void BaseFE::GetEdgeGlobDerivShFncAtIp(StdVector< Matrix<Double>* > & deriv, 
			       const Integer ip,
			       const Matrix<Double> & cornerCoords)
{
  ENTER_FCN( "BaseFE::GetEdgeGlobDerivShFncAtIp" );

  // vector of coordinates of the desired integration point
  Vector<Double> lCoord;
  lCoord.Resize(Dim_);

  for (Integer i=0; i<Dim_; i++)
    lCoord[i] = IntPoints_[ip-1][i];
  
  GetEdgeGlobalDerivShapeFnc(deriv, lCoord, cornerCoords);
}

void BaseFE::CalcEdgeShapeFncAtIp(Matrix<Double> & shape, 
				  const Integer ip,
				  const Matrix<Double> & cornerCoords)
{
  ENTER_FCN( "BaseFE::CalcEdgeShapeFncAtIp" );

  Vector<Double> lCoord;
  lCoord.Resize(Dim_);

  for (Integer i=0; i<Dim_; i++)
    lCoord[i] = IntPoints_[ip-1][i];
  
  CalcEdgeShapeFnc(shape, lCoord, cornerCoords);
}


void BaseFE::GetGlobalEdgeIndices(StdVector<Integer> & globEdgeIndex,
					  Integer * pDENodes, BaseSystem * algsys)
{
  ENTER_FCN( "BaseFE::GetGlobalEdgeIndices" );
  // define the global edge number
  for(Integer actEdge=0; actEdge < globEdgeIndex.GetSize(); actEdge++)
    globEdgeIndex[actEdge] = 
      algsys->GetNode2Edge(pDENodes[ edgeVertices_[actEdge][0]],
			   pDENodes[ edgeVertices_[actEdge][1]]);
}


} // end namespace CoupledField
