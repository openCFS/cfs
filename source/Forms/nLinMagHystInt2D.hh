// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINMAGHYSTINT2D
#define FILE_NLINMAGHYSTINT2D

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Materials/baseMaterial.hh>
#include <General/environment.hh>
#include <Forms/curlCurlNodeInt.hh>

namespace CoupledField {
  
  
  //! base class for calculation of linear elasticity
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

    //! returns B, just for postprocessing
    virtual void calcBMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord );


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
