#ifndef FILE_LINTHERMOMECHINT
#define FILE_LINTHERMOMECHINT


#include <Materials/baseMaterial.hh>
#include <General/environment.hh>

#include "Elements/basefe.hh"
#include "adbInt.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
	//! \f$ A \f$ is the matrix of the base function's derivatives of the mechanical field
	//! \f$ D \f$ is the matrix of thermal stress constants
	//! \f$ B \f$ is the matrix of the base functions of the heat conduction field

class LinThermoMechInt : public ADBInt
{
	public:

		/// Constructor
		LinThermoMechInt(BaseMaterial* matData, BaseMaterial *matDataMech, SubTensorType type);
		
		/// 
		~LinThermoMechInt();
		
	
		//!  Compute matrix \f$ A \f$ at given integration point.
		void calcAMat( Matrix<Double> &aMat, UInt ip,
							const Matrix<Double> &ptCoord );
	
		//!  Compute matrix \f$ D \f$ at given integration point.
		void calcDMat( Matrix<Double> &dMat );
	
		//!  Compute matrix \f$ B \f$ at given integration point.
		void calcBMat( Matrix<Double> &bMat, UInt ip,
							const Matrix<Double> &ptCoord );


		// 	QUERIES
	
	
		/*!   Query dimension of the matrix \f$ D \f$
	
		This method returns the dimensions of the data-matrix \f$ D \f$.
		Which corresponds to pyroelectric constants.
		In the case of 3D, a thermo-mechanic constants matrix has
		\f$ \mbox{dim}D=6\times 1 \f$
		In the case of 2D, a thermo-mechanic constants matrix has
		\f$ \mbox{dim}D=3\times 1 \f$ */
		void getDimD( UInt nRows, UInt nCols ) {
		//nRows = matDimRow_;
		//nCols = matDimCol_;
		};
	
		/*!   Query number of degrees of freedom for first physical quantity
	
		This method can be used to query the number of degrees of freedom at
		any node of a finite element for the physical quantity associated to
		the base functions whose derivatives form the matrix \f$ A \f$.
		In the case of 2D plane strain piezoelectric coupling the first
		physical quantity is the mechanical displacements.
		\return 2 */
		UInt getNumDofsA() {
		return numDofsA_;
		}
	
		/*!   Query number of degrees of freedom for second physical quantity
	
		This method can be used to query the number of degrees of freedom at
		any node of a finite element for the physical quantity associated to
		the base functions whose derivatives form the matrix \f$ B \f$.
		In the case of 2D plain strain piezoelectric coupling the second
		physical quantity is the electric potential.
		\return 1 */
		UInt getNumDofsB() {
		return numDofsB_;
		}
	
		//! Query material type for \f$ D \f$ tensor.
		MaterialType getDMaterialType();

		//! To say that the stiffness k of
		//! the coupling is thermo-mech
// 		bool getIsThermoMech(){
// 			return isThermoMech_;
// 		}

	private:

		UInt numDofsA_;
		UInt numDofsB_;
	
		// dimension of the D matrix
		UInt matDimRow_;
		UInt matDimCol_;
	
		// Stiffness Tensor from mech material pde
		Matrix<Double> cMatrix_;

		//bool isThermoMech_;


};


}

#endif // FILE_LINTHERMOMECHINT
