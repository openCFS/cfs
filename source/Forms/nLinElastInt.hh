// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINELASTINT_03
#define FILE_NLINELASTINT_03

#include "Forms/linElastInt.hh"

#include "General/defs.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"

namespace CoupledField {
class BaseFE;
class BaseMaterial;
class EntityIterator;
struct Elem;
template <class TYPE> class NodeStoreSol;
template <class TYPE> class Vector;
}  // namespace CoupledField

namespace CoupledField
{
  


  // =============================================================================
  // base class for nonlinear mechanics
  // =============================================================================

  /// base class for calculation of nonlinear elasticity
  class nLinElastInt : public linElastInt
  {
  public:

    /// Constructor
    nLinElastInt(BaseMaterial* matData);
  
    /// Destructor
    virtual ~nLinElastInt();  

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );


  //! Helper function for setting the current entity
  void SetActEntities( EntityIterator ent1,
                       EntityIterator ent2 );

protected:    
  
  virtual void calcDMat(Matrix<Double> & dMat)
  {
    EXCEPTION("child implementation missing");
  }
  
  /** We do not provide a SIMP variant!
   * @see linElastInt::calcDMat(Matrix<Double>, const Elem*) */
  void calcDMat(Matrix<Double> &dMat, const Elem* elem)
  {
    calcDMat(dMat); 
  };


    
    /** @see BaseForm::CalcBMat() */
    virtual void CalcBMat(Matrix<Double> & bMat, UInt ip,
                          const Matrix<Double> & ptCoord );

    /// displacement of all nodes of actual element
    Matrix<Double> elemDisp_;

  };
  
  
 

  // =============================================================================
  // 3D nonlinear mechanics
  // =============================================================================


  /// class for calculation of 3d nonlinear elasticity
                       // first part: nonlinear B-matrix
  class nLinMech3dInt_BNonLin : public nLinElastInt
  {
  public:

    /// Constructor
    nLinMech3dInt_BNonLin(BaseMaterial* matData);

  
    /// Destructor
    virtual ~nLinMech3dInt_BNonLin();  
  
  protected:  
    /// returns D - matrix for BDB
    virtual void calcDMat(Matrix<Double> & dMat);

    /// returns dimension of D matrix
    virtual UInt getDimD(){return 6;};
  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 3;};
  };


  /// class for calculation of 3d nonlinear elasticity
                                 // second part: regarding internal stresses
  class nLinMechInt_PiolaStress : public nLinElastInt
  {
  public:
    friend class nLinMech_linFormInt;
  
    /// Constructor
    nLinMechInt_PiolaStress(BaseMaterial* matData);
  
    /// Destructor
    virtual ~nLinMechInt_PiolaStress();  
  
  protected:  
    /// returns D - matrix for BDB (size: 9x9, contains the 2. Piola-Kirchhoff-Stress tensor!!)
    virtual void calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord);
  
    /** @see BaseForm::CalcBMat() */
    virtual void CalcBMat(Matrix<Double> & bMat, UInt ip, const Matrix<Double> & ptCoord);

    /// returns dimension of D matrix
    virtual UInt getDimD(){return piolaDimD_;};
  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs()=0;  


  protected:
    /// sets the size of d-matrix (needed for lin and nonlin B-matrix)
    virtual void setPiolaDimD(UInt actDim){piolaDimD_ = actDim;};

    /// calculates Piola-Kirchoff-stresses (vector notation)
    virtual void CalcStressVec(Vector<Double>& piolaStressVec, UInt ip, 
                               BaseFE *elem,
                               Matrix<Double> & ptCoord,
                               Matrix<Double> & elemDisp );  

    /// returns linear B - matrix
    virtual void calcLinBMat(Matrix<Double> & bMat, UInt ip, 
                             BaseFE *elem, Matrix<Double> & ptCoord );

    /// returns nonlinear B - matrix
    virtual void calcNonLinBMat(Matrix<Double> & bMat, UInt ip, 
                                BaseFE *elem, Matrix<Double> & ptCoord,
                                Matrix<Double> & elemDisp );

    /// returns material D-matrix for 3d mechanics
    virtual void calcMaterialDMat(Matrix<Double> & dMat);

    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize()=0;

    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize()=0;

    /// conversion of stress vector to stress tensor
    virtual void convertStressVecToTensor(Matrix<Double>& stressTensor,
                                               Vector<Double>& piolaStress);
  
  private:
  
    /// dimension of d-matrix (has to be changed for some dirty implementation features ... )
    UInt piolaDimD_;
  };
  


  /// class for calculation of 3d nonlinear elasticity
                                  // second part: regarding internal stresses
  class nLinMech3dInt_PiolaStress : public nLinMechInt_PiolaStress
  {
  public:

    /// Constructor
    nLinMech3dInt_PiolaStress(BaseMaterial* matData);
  
    /// Destructor
    virtual ~nLinMech3dInt_PiolaStress();  
  
  protected:  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 3;};


    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize(){return 6;};


    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize(){return 9;};
  };
  





  // =============================================================================
  // nonlinear plane strain mechanics (2D)
  // =============================================================================

  /// class for calculation of plane strain nonlinear elasticity
                                    // first part: nonlinear B-matrix
  class nLinMechPlaneStrainInt_BNonLin : public nLinElastInt
  {
  public:

    /// Constructor
    nLinMechPlaneStrainInt_BNonLin(BaseMaterial* matData);

  
    /// Destructor
    virtual ~nLinMechPlaneStrainInt_BNonLin();  
  
  protected:  
    /// returns D - matrix for BDB
    virtual void calcDMat(Matrix<Double> & dMat);

    /// returns dimension of D matrix
    virtual UInt getDimD(){return 3;};
  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};  
  };




  /// class for calculation of 2d nonlinear elasticity
                                         // second part: regarding internal stresses
  class nLinMechPlaneStrainInt_PiolaStress : public nLinMechInt_PiolaStress
  {
  public:

    /// Constructor
    nLinMechPlaneStrainInt_PiolaStress(BaseMaterial* matData);
  
    /// Destructor
    virtual ~nLinMechPlaneStrainInt_PiolaStress();  
  
  protected:  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};  

    /// conversion of stress vector to stress tensor
    virtual void convertStressVecToTensor(Matrix<Double>& stressTensor,
                                               Vector<Double>& piolaStress);

    /// returns material D-matrix for 3d mechanics
    virtual void calcMaterialDMat(Matrix<Double> & dMat);

    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize(){return 3;};

    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize(){return 4;};
  };
  







  // =============================================================================
  // nonlinear axisymmetrical mechanics
  // =============================================================================

  /// class for calculation of axisymmetric nonlinear elasticity
                                             // first part: nonlinear B-matrix
  class nLinMechAxiInt_BNonLin : public nLinElastInt
  {
  public:

    /// Constructor
    nLinMechAxiInt_BNonLin(BaseMaterial* matData);

  
    /// Destructor
    virtual ~nLinMechAxiInt_BNonLin();  
  
  protected:  
    
    /// returns D - matrix for BDB
    virtual void calcDMat(Matrix<Double> & dMat);
    
    /// returns dimension of D matrix
    virtual UInt getDimD(){return 4;};
  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};  
  };




  /// class for calculation of 2d axisymmetric nonlinear elasticity
                                 // second part: regarding internal stresses
  class nLinMechAxiInt_PiolaStress : public nLinMechInt_PiolaStress
  {
  public:

    /// Constructor
    nLinMechAxiInt_PiolaStress(BaseMaterial* matData);
  
    /// Destructor
    virtual ~nLinMechAxiInt_PiolaStress();  
  
  protected:  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};  

    /// conversion of stress vector to stress tensor
    virtual void convertStressVecToTensor(Matrix<Double>& stressTensor,
                                               Vector<Double>& piolaStress);

    /// returns material D-matrix for 3d mechanics
    virtual void calcMaterialDMat(Matrix<Double> & dMat);

  protected:
    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize(){return 4;};

    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize(){return 5;};  
  };
  






  // =============================================================================
  // prestress
  // =============================================================================



  /// class for regarding 3d prestress
  class PreStressInt : public nLinMechInt_PiolaStress
  {
  public:
    // preStressLinFormInt uses calcDMat from this class
    friend class PreStressLinFormInt;
  
    //! Constructor
    PreStressInt( BaseMaterial* matData );
  
    /// Destructor
    virtual ~PreStressInt();  

    //! Compute element matrix associated to BDB form
    //! needed in case of not constant prestress
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
    //! computation of D-matrix (consists of stress values)
    void calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord);
  
    /** @see BaseForm::CalcBMat() */
    void CalcBMat(Matrix<Double> & bMat, UInt ip, const Matrix<Double> & ptCoord);

    //! sets the prestresses (in case of given prestresses) 
    void SetPreStress(Vector<Double>& piolaStressVec );
    
    //!
    void convertStressVecToTensor( Matrix<Double> &stressTensor,
                                   Vector<Double> &piolaStress ) {
      EXCEPTION( "convertStressVecToTensor not implemented" );
    }
    
    //! sets the current displacement
    void SetMechDisp( NodeStoreSol<Double>& sol ) {
      sol_ = &sol;}

    //! sets the current displacement (we just need the real part)
    void SetMechDisp( NodeStoreSol<Complex>& sol ) {
      complexSol_ = &sol; 
      isSolComplex_ = true; }


  protected:

    /// calculates constant pre-stresses (vector notation)
    void CalcConstPreStressVec( Vector<Double>& piolaStressVec, 
                                Vector<Double>& piolaStressVal  );

    //! calculates the actual prestressing for the current element
    virtual void CalcActStressVec(Vector<Double>& StressVec, UInt ip, 
                                  BaseFE *elem, Matrix<Double> & ptCoord,
                                  Matrix<Double> & elemDisp );  

    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize(){
      EXCEPTION( "Not implemented here" );
      return 0;
    }

    /// returns dimension of D matrix
    virtual UInt getDimD() {
      EXCEPTION( "Not implemented here" );
      return 0;
    }
  
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){
      EXCEPTION( "Not implemented here" );
      return 0;
    }

    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize() {
      EXCEPTION( "Not implemented here" );
      return 0;
    }


  private: 
    // true, if constant prestress
    bool isPreStressConst_;

    //! stores the constant prestressing tensor
    Matrix<Double> preStressTensor_;

    //! yes, if displacement used for computing the prestress is complex
    bool isSolComplex_;

    //! stores complex mechanical displacement
    NodeStoreSol<Complex>* complexSol_;

    //! check, if CalcElementMatrix is used for the first time
    bool isVirgin_;
    
  };


  /// class for regarding 3d prestress
  class PreStressInt3D : public PreStressInt
  {
  public:
  
    //! Constructor
    PreStressInt3D(BaseMaterial* matData );
  
    /// Destructor
    virtual ~PreStressInt3D();  

    //!
    void convertStressVecToTensor(Matrix<Double>& stressTensor,
                                  Vector<Double>& piolaStress);

  protected:
    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize(){return 9;};

    /// returns dimension of D matrix
    virtual UInt getDimD(){return 6;};

    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 3;};

    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize(){return 6;};

  };


  /// class for regarding 2d prestress in plane strain case
  class PreStressIntPlaneStrain : public PreStressInt
  {
  public:
    //! Constructor
    PreStressIntPlaneStrain(BaseMaterial* matData );
  
    /// Destructor
    virtual ~PreStressIntPlaneStrain();  

    //!
    void convertStressVecToTensor(Matrix<Double>& stressTensor,
                                      Vector<Double>& piolaStress);

  protected:
    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize(){return 4;};

    /// returns dimension of D matrix
    virtual UInt getDimD(){return 3;};

    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};

    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize(){return 3;};

  };


  /// class for regarding 2d axi prestress
  class PreStressIntAxi : public PreStressInt
  {
  public:
  
    //! Constructor
    PreStressIntAxi(BaseMaterial* matData );
  
    /// Destructor
    virtual ~PreStressIntAxi();  

    //!
    void convertStressVecToTensor(Matrix<Double>& stressTensor,
                                      Vector<Double>& piolaStress);

  protected:
    /// returns the size of the full piola d-matrix
    virtual UInt getFullPiolaDMatSize(){return 5;};

    /// returns dimension of D matrix
    virtual UInt getDimD(){return 4;};

    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};

    /// returns the size of the material d-matrix
    virtual UInt getMaterialDMatSize(){return 4;};

  };

  
} //end namespace

#endif // FILE_XXX


