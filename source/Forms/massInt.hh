// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MASSINT_1
#define FILE_MASSINT_1

#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation  element mass matrix
  class MassInt : public BaseForm
  {
  public:

    // Constructor
    MassInt(const Double aDensity, const UInt nrDofsPerNode=1, 
            bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~MassInt();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2,
                            const DesignElement::Type direction);
    
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2) {
      CalcElementMatrix(elemMat, ent1, ent2, DesignElement::NO_DERIVATIVE);
    }
      
    //! Sets a multiplicative factor for element matrix
    void SetSecondFactor( const std::string& factor )
    {mParser_->SetExpr( mHandle_, factor );};

    //! define diagonal mass matrix
    void SetDiagMass() {
      diagMass_ = true;
    };

  protected:
    
    //! generates a multi-dof-matrix with similar entries for all dofs
    virtual void MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  
                              const UInt nrDofs);

    virtual void MassMultiDofZero(Matrix<Double>& massMultDofZero, 
                                  const Matrix<Double>& massMatSingleDof);
    
    double GetErsatzMaterialMass(const Elem* elem, DesignElement::Type direction);
  
  private:

    Double density_;          //!< multiplicative value for mass integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
    bool diagMass_;         //<! true, mass matrix is diagonal
  };

}

#endif // FILE_MASSINT
