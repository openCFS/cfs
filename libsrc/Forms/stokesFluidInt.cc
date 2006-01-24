#include <iostream>
#include <fstream>

#include "stokesFluidInt.hh"

namespace CoupledField
{

  StokesFluidInt::StokesFluidInt(Double aVal, Double bVal)
    : BaseForm(),density_ (aVal),dynamicViscosity_ (bVal)
  {
    ENTER_FCN( "StokesFluidInt::StokesFluidInt" );
  }

  StokesFluidInt::StokesFluidInt()
    : BaseForm()
  {
    ENTER_FCN( "StokesFluidInt::StokesFluidInt" );
  }


 
  StokesFluidInt::~StokesFluidInt()
  {
    ENTER_FCN( "StokesFluidInt::~StokesFluidInt" );
  }


// The element matrix must be resorted because during matrix setup
// the local unknown vector is ordered as: 
//[U1 U2 U3 U4 V1 V2 V3 V4 P1 P2 P3 P4 Om1 Om2 Om3 Om4]^T
// for Node 1...4
// the needed order is
//[U1 V1 P1 Om1 U2 V2 P2 Om2 U3 V3 P3 Om3 U4 V4 P4 Om4 ]^T
  void StokesFluidInt::ResortElementMatrix(Matrix<Double> & sortedElemMat, 
                                      const Matrix<Double> & elemMat, 
                                      const UInt & nrOfNodes, 
                                      const UInt & dofPerNode) {
    ENTER_FCN( "StokesFluidInt::ResortElementMatrix" );

    UInt i, j, k, dofPerElem;
    Matrix<Double> auxMat;

    dofPerElem = nrOfNodes * dofPerNode;
    sortedElemMat.Resize(dofPerElem); sortedElemMat.Init();
    auxMat.Resize(dofPerElem); auxMat.Init();


    for (j=0; j<dofPerElem; j+=dofPerNode)
      {
	    for (i=0; i<dofPerElem; i++)
	      {
	        for (k=0; k<dofPerNode; k++)
	          {
                auxMat[i][j+k] = elemMat[i][j/dofPerNode + k*nrOfNodes];
	          }
	      }
      }

   for (i=0; i<dofPerElem; i+=dofPerNode)
     {
	   for (j=0; j<dofPerElem; j++)
	     {
	       for (k=0; k<dofPerNode; k++)
	         {
	           sortedElemMat[i+k][j] = auxMat[i/dofPerNode + k*nrOfNodes][j];
	         }
	     }
     }
  }


  void StokesFluidInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "StokesFluidInt::Print"); 
    (*out)<< "StokesFluid  stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
