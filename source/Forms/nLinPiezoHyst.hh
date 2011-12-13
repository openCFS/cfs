// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NONLIN_PIEZO_HYST
#define FILE_NONLIN_PIEZO_HYST

#include "Forms/linGradBDBInt.hh"

#include "Domain/entityList.hh"
#include "Forms/linPiezoCoupling.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"


namespace CoupledField {

class BaseMaterial;
class EqnMap;
class Grid;
  class PiezoCoupling;
class StdPDE;
struct Elem;
struct ResultInfo;
template <class TYPE> class GradientFieldOp;
  template<typename> class MechStressStrain;
template <typename T> class NodeStoreSol;

  class nLinPiezoHystCouple : public linPiezoCoupling {

  public:

    // friend class declaration
   friend class ADBInt;

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Default constructor
    nLinPiezoHystCouple(BaseMaterial* matDataCouple, BaseMaterial* matDataMech, 
                        BaseMaterial* matDataElec, SubTensorType type,
                        bool isMech = true); 

    //! Destructor
   ~nLinPiezoHystCouple();

    //@}

   //! Compute the nonlinear data-matrix \f$D\f$
    void calcDMat(Matrix<Double> & dMat );

    /** We do not provide a SIMP variant!
     * @see linElastInt::calcDMat(Matrix<Double>, const Elem*) */
    void calcDMat(Matrix<Double> &dMat, const Elem* elem)
    {
      calcDMat(dMat); 
    }; 
    
   //! set objects for computation of E-field
   void Set4NonLinMaterial( Grid* ptGrid, 
                            StdPDE* ptPDE2Elec,
                            shared_ptr<EqnMap> eqnMapElec,
                            shared_ptr<ResultInfo> resultElec,
                            PiezoCoupling* ptPiezoCoupling);

    //!
    void CalcElementMatrix( Matrix<Double>& elemMat,
                           EntityIterator& ent1, 
                           EntityIterator& ent2 );

    //! set solution object of PDE elec
    void SetSolutionElec(NodeStoreSol<Double>& solhelp){
      solElec_= &solhelp;
    }

    //! set previous solution object of PDE elec
    void SetPrevSolutionElec(NodeStoreSol<Double>& solPrevHelp){
      solPrevElec_= &solPrevHelp;
    }
   

    //! set solution object of PDE mech
    void SetSolutionMech(NodeStoreSol<Double>& solhelp){
      solMech_= &solhelp;
    }

    //! set prveious solution object of PDE mech
    void SetPrevSolutionMech(NodeStoreSol<Double>& solPrevHelp){
      solPrevMech_= &solPrevHelp;
    }
   

    void getDimD( UInt nRows, UInt nCols ) {
      nRows = matDimRow_;
      nCols = matDimCol_;
    };

    UInt getNumDofsA() {
      return numDofsA_;
    }


    UInt getNumDofsB() {
      return numDofsB_;
    }


  private:

    //!
    bool isMech_;
    //!
    PiezoCoupling* piezoCoupling_;
   
    /// scalar electric potential of all nodes of actual element
    Vector<Double> elemPot_;
    Matrix<Double> elemDispl_;
    Vector<Double> elemPotPrev_;
    Matrix<Double> elemDisplPrev_;
    
    GradientFieldOp<Double>* EfieldOp_;
    GradientFieldOp<Double>* EfieldPrevOp_;
    MechStressStrain<Double>* mechStrainOp_;

    // local copy of entity iterator
    EntityIterator ent1_;
    
    UInt numDofsA_;
    UInt numDofsB_;
    
    UInt matDimRow_;
    UInt matDimCol_;
    
    NodeStoreSol<Double> * solElec_;
    NodeStoreSol<Double> * solMech_;
    NodeStoreSol<Double> * solPrevElec_;
    NodeStoreSol<Double> * solPrevMech_;
    
    SubTensorType subTensorType_;

    BaseMaterial* matDataCouple_;
    BaseMaterial* matDataMech_;
    BaseMaterial* matDataElec_;
  };



  //! electric stiffnes bilinear form in case of piezoelectric hysteresis
  class nLinPiezoHystElec : public linGradBDBInt {

    
  public:

    //@{ \name Construction and destruction

    //! Default constructor
    nLinPiezoHystElec(BaseMaterial* matDataCouple, BaseMaterial* matDataMech, 
                      BaseMaterial* matDataElec, SubTensorType type); 

    //! Destructor
    ~nLinPiezoHystElec();

    //@}

   //! Compute the nonlinear data-matrix \f$D\f$
    void calcDMat(Matrix<Double> & dMat );

    /** We do not provide a SIMP variant!
     * @see linElastInt::calcDMat(Matrix<Double>, const Elem*) */
    void calcDMat(Matrix<Double> &dMat, const Elem* elem)
    {
      calcDMat(dMat); 
    }; 
    
    
   //! set objects for computation of E-field
   void Set4NonLinMaterial(Grid* ptGrid, 
                           StdPDE* ptPDE2Elec,
                           shared_ptr<EqnMap> eqnMapElec,
                           shared_ptr<ResultInfo> resultElec,
                           PiezoCoupling* ptPiezoCoupling);

   //!
   void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

    //! set solution object of PDE elec
    void SetSolutionElec(NodeStoreSol<Double>& solhelp){
      solElec_= &solhelp;
    }

    //! set previous solution object of PDE elec
    void SetPrevSolutionElec(NodeStoreSol<Double>& solPrevHelp){
      solPrevElec_= &solPrevHelp;
    }
   

    //! set solution object of PDE mech
    void SetSolutionMech(NodeStoreSol<Double>& solhelp){
      solMech_= &solhelp;
    }

    //! set prveious solution object of PDE mech
    void SetPrevSolutionMech(NodeStoreSol<Double>& solPrevHelp){
      solPrevMech_= &solPrevHelp;
    }
   
  private:
    //!
    PiezoCoupling* piezoCoupling_;
   
    /// scalar electric potential of all nodes of actual element
    Vector<Double> elemPot_;
    Matrix<Double> elemDispl_;
    Vector<Double> elemPotPrev_;
    Matrix<Double> elemDisplPrev_;
    
    GradientFieldOp<Double>* EfieldOp_;
    GradientFieldOp<Double>* EfieldPrevOp_;
    MechStressStrain<Double>* mechStrainOp_;

    // local copy of entity iterator
    EntityIterator ent1_;
    
    NodeStoreSol<Double> * solElec_;
    NodeStoreSol<Double> * solMech_;
    NodeStoreSol<Double> * solPrevElec_;
    NodeStoreSol<Double> * solPrevMech_;
    
    SubTensorType subTensorType_;

    BaseMaterial* matDataCouple_;
    BaseMaterial* matDataMech_;
    BaseMaterial* matDataElec_;
  };

}

#endif
