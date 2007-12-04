#ifndef SURFACENORMALINT_HH_
#define SURFACENORMALINT_HH_

// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <General/environment.hh>

namespace CoupledField 
{
  /** This class calculates the surface normal matrix as given by 
   * J. Du and N. Olhoff, "Minimization of sound radiation from vibrating bi-material structures 
   * using topology optimization", Struct Multidisc. Optim, 2007.
   * It Implements the element matrices $S_{ne}$ which are devined as
   * $S_{ne}=\int_{S_e} N^T n n^T N ds$ */
  class SurfaceNormalInt : public BDBInt 
  {
  public:
    SurfaceNormalInt( BaseMaterial* matData, SubTensorType type = FULL );

    virtual ~SurfaceNormalInt();

    /** Does a PostInit in the sense that we take a sample element to setup our 
     * corresponding volume base element reference */
    void PostInit(Elem* elem);
    
    //! Query material type for \f$D\f$ tensor
    MaterialType getDMaterialType() { return MECH_STIFFNESS_TENSOR; }

    void ExtractElemInfo(const EntityIterator& it)
    {
      BaseForm::ExtractElemInfo(const_cast<EntityIterator&>(it));
    }

     //! calculates the material data 
    void calcDMat( Matrix<Double> &dMat )
    {
      CalcNormalMatrix<Double>(dMat);
    }

     //! calculates the material data 
    void calcDMat( Matrix<Complex> &dMat )
    {
      CalcNormalMatrix<Complex>(dMat);
    }

    //! returns B - matrix for BDB
    virtual void calcBMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord );


    //! returns dimension of D matrix
    virtual UInt getDimD() {
      return dimD_; 
    };
    
    //! returns nr. of degrees of freedom
    virtual UInt getNrDofs() {
      return nrDofs_;
    };
    
  private:

   /** calculates normal matrix for the z-direction hardcoded! */
   template <class T>
   void CalcNormalMatrix( Matrix<T> &out );
    
    
    //dimension of Dmatrix
    UInt dimD_;
    
    //! number of degrees 
    UInt nrDofs_;

    //! subtype of tensor
    SubTensorType subTensorType_;

    //! number of incompatible modes
    UInt nrICModes_;
    
    UInt nrDofsPerNode_;   //!< degrees of freedom per node

    /** This is our corrsponding volume element set in PostInit */
    BaseFE* volElem;
    
  };
  



} //end namespace

#endif /*SURFACENORMALINT_HH_*/
