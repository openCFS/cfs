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
   
    //! Constructor with pointer to BaseElem
    FlatShellInt(BaseMaterial *matData);
   
    //!  Destructor
    virtual ~FlatShellInt();

    //! Set Thickenss
    void SetThickness( Double thickness ) { thickness_ = thickness; }
   
    //! Set Penalty value for drilling dof
    void SetPenaltyDof( Double penaltyDof ) { penaltyDof_ = penaltyDof; }
   
    //!returns dimension of D matrix of bending stiffness
    UInt getDimD(){return 8;};
   
    //!returns nr. of degrees of freedom
    UInt getNrDofs(){return 6;};
    
    //!Function for calculation bdb matrix 
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
   
    //Enumeration for differentiation of the three mechanical parts Membrane, Bending, Shear
    typedef enum {MEMBRANE,BENDING,SHEAR,COUPLED1,COUPLED2} SubPartType;

   
  protected:

    // Used for coordinate transformations
    void CoordTrans( const Matrix<Double> &ptCoord, Matrix<Double> &TransMat, Matrix<Double> &ShellCoord );
    
    //Used for going from local to Global 
    void LocaltoGlob( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat );
    
    //!calculates the material Matrix
    void calcDMat( Matrix<Double> &dMat);
    
    //!returns B - matrix for BDB
    void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord, SubPartType parts );
           
    bool updateDMatInEveryIP_;
    
    StdVector<SubPartType> parts_;
    StdVector<bool> reducedPart_;
    
    //Parametres of the lamina structures
    Vector<Double> z_;
    Vector<Double> orAngle_;
    
    //! Thickenss of the element
    Double thickness_;
    
    //! Penalty value for drilling dof
    Double penaltyDof_;
  };

}//end of namespace

#endif 
