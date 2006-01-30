#ifndef FILE_NRBCINT_1
#define FILE_NRBCINT_1

#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation  element nrbc matrices 
  //! having a $(elemNod \times J-1) \times (elemNod \times J-1)$ dimension
  class NrbcInt : public BaseForm
  {
  public:
    // Constructor
    NrbcInt(BaseFE * aptelemt, Double aDensity, UInt nrbcMatType, Boolean axi=FALSE);

    // Constructor
    NrbcInt(const Double aDensity, const UInt nrDofsPerNode=1, Boolean axi=FALSE);

    // Constructor
    NrbcInt(const Double aDensity, const UInt nrDofsPerNode, UInt nrbcMatType, Boolean axi=FALSE);

    // Destructor
    virtual ~NrbcInt();

    //! Calculation of stiffmess matrix
    void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMa);
  
    //! set additional multiplicative factor of nrbc matrix
    void SetFactor(Double aFactor) 
    {factor_ = aFactor;};
   
    virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

  protected:
    //! generates a multi-dof-matrix with similar entries for all dofs
    virtual void NRBCMultiDof(Matrix<Double>& nrbcMultDof,
                              const Matrix<Double>& nrbcMatSingleDof,  
                              const UInt nrDofs);

    //! generates a multi-dof-matrix with the entries shifted (j + 1) and zero 
    //! entries for the last dof (j) rows
    virtual void NRBCMultiDof_Rmat(Matrix<Double>& nrbcMultDof_Rmat, 
                                  const Matrix<Double>& nrbcMatSingleDof);

    //! generates a multi-dof-matrix with the entries shifted (j - 1) and zero 
    //! entries for the first dof (j) rows
    virtual void NRBCMultiDof_PorQmat(Matrix<Double>& nrbcMultDof_PorQmat, 
                                  const Matrix<Double>& nrbcMatSingleDof);  

    void LocaltoGlob( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat );

  private:

    Double density_;          //!< multiplicative value for nrbc integrator
    Double factor_;           //!< yet another multiplicative value for nrbc integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
    UInt nrbcMatType_;     //!< value to know which type of matrix should assembled
    
  };





}

#endif // FILE_NRBCINT
