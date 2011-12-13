// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_POLMAT_RHS_INT
#define FILE_POLMAT_RHS_INT

#include "Forms/linSurfForm.hh"
#include "General/defs.hh"
#include "MatVec/vector.hh"

namespace CoupledField {
class EntityIterator;
struct Elem;
}  // namespace CoupledField

namespace CoupledField 
{
  
  //! Class implementing an inhomogeneous Neumann integrator for mech part of polarization matrix
  class PiezoPolarizationMatrixMechRHSInt : public LinearSurfForm
  {

  public:
    
    //! Ctor
    PiezoPolarizationMatrixMechRHSInt(const Vector<Double> &vals, const int num);
    
    //! Destructor
    virtual ~PiezoPolarizationMatrixMechRHSInt();

    /// Calculation of RHS vector for double entries, i.e. transient and static 
    void CalcElemVector( Vector<Double> & elemVec, EntityIterator& ent );

    /** number of the column of the material matrix */
    int number;
    
  protected:

    /// comprises part of CalcElemVector, which is equal for double and complex
    void PrepareElemVector( Vector<Double> & elemVec, EntityIterator& ent );
    
  private:
    // forbid standard ctor
    PiezoPolarizationMatrixMechRHSInt() {}
        
    Vector<Double> rhsvals;
    
    // pointer to the adjacent volume element
    Elem * ptVolElem_;

  };
  
  //! Class implementing an inhomogeneous Neumann integrator for elec part of polarization matrix
  class PiezoPolarizationMatrixElecRHSInt : public LinearSurfForm
  {

  public:
    
    //! Ctor
    PiezoPolarizationMatrixElecRHSInt(const Vector<Double> &vals, const int num);
    
    //! Destructor
    virtual ~PiezoPolarizationMatrixElecRHSInt();

    /// Calculation of RHS vector for double entries, i.e. transient and static 
    void CalcElemVector( Vector<Double> & elemVec, EntityIterator& ent );

    /** number of the column of the material matrix */
    int number;
        
  protected:

    /// comprises part of CalcElemVector, which is equal for double and complex
    void PrepareElemVector( Vector<Double> & elemVec, EntityIterator& ent );
    
  private:
    // forbid standard ctor
    PiezoPolarizationMatrixElecRHSInt() {}
        
    Vector<Double> rhsvals;
    
    // pointer to the adjacent volume element
    Elem * ptVolElem_;

  };

} // end of namespace
#endif
