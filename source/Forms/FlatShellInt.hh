// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FLATSHELLINT
#define FILE_FLATSHELLINT

#include "baseForm.hh"
#include "Domain/Composite.hh"


namespace CoupledField
{
  //! Class for flat shell Integrator
  class FlatShellInt : public BaseForm
  {

    public:

    //! Constructor with pointer to BaseMaterial
    FlatShellInt(BaseMaterial *matData, bool hasDrillDof );
   
   //!Constructor with pointer to Composite
    FlatShellInt(Composite * composite, bool hasDrillDof );
       
    //!  Destructor
    virtual ~FlatShellInt();

    //! Set Thickenss
    void SetThickness( Double thickness ) { thickness_ = thickness; }  
    
    //! Set Penalty value for drilling dof
    void SetPenaltyDof( Double penaltyDof ) { penaltyDof_ = penaltyDof; }
   
    //!returns nr. of degrees of freedom
    virtual Integer  getNrDofs(){
      if ( hasDrillDof_) {
        return 6;
      } else {
        return 5;
      }
    };
    

  protected:

    // Used for coordinate transformations
    void CoordTrans( const Matrix<Double> &ptCoord, Matrix<Double> &TransMat, Matrix<Double> &ShellCoord );
    
    //Used for going from local to Global 
    void LocaltoGlob( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat );
    
    //Used for going from local to Global for the Piezo Coupling terms 
    void LocaltoGlobPiezo( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat );
    
     //! Thickenss of the element
    Double thickness_;

    //! Penalty value for drilling dof
    Double penaltyDof_;
    
    //! Flags
    bool isComposite_;
    bool isOrthotropic_;
    
    //! Pointer to composite material
    Composite * composite_;

   // flag indicating if drilling dof is needed
   bool hasDrillDof_;
 };

}//end of namespace

#endif 
