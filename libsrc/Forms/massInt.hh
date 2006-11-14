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
                            EntityIterator& ent2 );
      

    //! define diagonal mass matrix
    void SetDiagMass() {
      diagMass_ = true;
    };

    //! for fractional damping model, is called in AcousticPDE::DefineIntegrators
    void SetFracDamping() 
    {isFracDamping_ = true;};

    //! for fractional damping model, is called in AcousticPDE::DefineIntegrators
    void UnsetFracDamping() 
    {isFracDamping_ = false;};
  
    //! set additional multiplicative factor of mass matrix
    void SetFactor(Double aFactor) 
    {factor_ = aFactor;};
   
  protected:
    
    //! generates a multi-dof-matrix with similar entries for all dofs
    virtual void MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  
                              const UInt nrDofs);

    virtual void MassMultiDofZero(Matrix<Double>& massMultDofZero, 
                                  const Matrix<Double>& massMatSingleDof);
  
  private:

    Double density_;          //!< multiplicative value for mass integrator
    Double factor_;           //!< yet another multiplicative value for mass integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
    bool diagMass_;         //<! true, mass matrix is diagonal
  };

}

#endif // FILE_MASSINT
