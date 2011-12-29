// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINMAGHYSTINT2D
#define FILE_NLINMAGHYSTINT2D

#include "Forms/bdbInt.hh"
#include "General/environment.hh"

#include "General/defs.hh"
#include "General/exception.hh"
#include "MatVec/vector.hh"

namespace CoupledField {
  
  
  //! base class for calculation of linear elasticity
class BaseMaterial;
class CurlCurlNode2DInt;
class EntityIterator;
template <class TYPE> class Matrix;

  class nLinMagHystInt2D : public BDBInt {

  public:

    //! Constructor
    nLinMagHystInt2D( BaseMaterial* matData, bool axi, bool coordUpdate );

    //! Destructor
    virtual ~nLinMagHystInt2D();

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

    /** @see BaseForm::CalcBMat() */
    virtual void CalcBMat( Matrix<Double> &bMat, UInt ip,
                           const Matrix<Double> &ptCoord );


  protected:    

     //! calculates the material data 
    void calcDMat( Matrix<Double> &dMat, UInt elNr );

    //! returns B - matrix for BDB
    void calcMyBMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord );



    //! returns dimension of D matrix
    virtual UInt getDimD() { return 2; };

   //! returns nr. of degrees of freedom
    virtual UInt getNrDofs() { return 1; };

    //! Query material type for \f$D\f$ tensor
    virtual MaterialType getDMaterialType() {
      EXCEPTION("Not applicable for 2D magnetic with hysteresis");
      return NO_MATERIAL; }


  private:

    Double startmatVal_;       //!<  start value for reluctivity
    Vector<Double> magPot_;    //!< magnetic vector potential at nodes
    Double reluctivity0_;
    Vector<Double> Bfield_;

    CurlCurlNode2DInt*  ptBMat_; 
  };
  



} //end namespace

#endif // FILE_LINELASTINT
