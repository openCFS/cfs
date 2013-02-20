#ifndef FILE_ICMODES_INT_HH
#define FILE_ICMODES_INT_HH

#include "BDBInt.hh"

namespace CoupledField {

  //! Bilinear form for linear elasticity with Incompatible Modes
  
  //! This Bilinear form calculates the linear elasticity bilinear form
  //! with additional incompatible modes of Taylor Wilson.
  //! Thus, it adds internally so-called incompatible bubble modes,
  //! which can be eliminated directly on the element level, i.e.
  //! the final element matrix has the same size as a standard approximation.
  //! \tparam COEF_DAATA_TYPE Data type of the material tensor  
  //! \tparam B_DATA_TYPE Data type of the differential operator
  template< class COEF_DATA_TYPE = Double, 
            class B_DATA_TYPE = Double>
  class ICModesInt : public BDBInt<COEF_DATA_TYPE, B_DATA_TYPE> {
  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor with differential operators
    
    //! This initializes the bilinearform
    //! \param bOp Standard differential operator (strainOp)
    //! \param gOp Strain-type differential operator for incompatible modes
    ICModesInt( BaseBOperator * bOp, 
                BaseBOperator * gOp,
                PtrCoefFct dData, MAT_DATA_TYPE factor,
                bool coordUpdate = false );

    //! Destructor
    virtual ~ICModesInt(){
      delete gOperator_;
    }
    
    //! Compute element matrix associated to AB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );

  protected:

    //! Differential operator for incompatible modes
    BaseBOperator* gOperator_;
    
    //! Store intermediate operator matrix for G
    Matrix<MAT_DATA_TYPE> gMat_;
    
    //! Store intermediate operator matrix for dG
    Matrix<MAT_DATA_TYPE> dgMat_;

  };
}


#endif
