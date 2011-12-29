// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NONLIN_PIEZO_COUPLING
#define FILE_NONLIN_PIEZO_COUPLING


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
class StdPDE;
struct ResultInfo;
template <class TYPE> class GradientFieldOp;
template <typename T> class NodeStoreSol;

  class nLinPiezoCoupling : public linPiezoCoupling {

    
  public:

    // friend class declaration
   friend class ADBInt;

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Default constructor
   nLinPiezoCoupling( BaseMaterial* matData, BaseMaterial* matDataMech, 
                      BaseMaterial* matDataElec, SubTensorType type); 

    //! Destructor
   ~nLinPiezoCoupling();

    //@}

   //! Compute the nonlinear data-matrix \f$D\f$
   void calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord);

   //! set objects for computation of E-field
   void Set4NonLinMaterial(Grid* ptGrid, 
                           StdPDE* ptPDE,
                           shared_ptr<EqnMap> eqnMap,
                           shared_ptr<ResultInfo> result);

   //!
   void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

   void SetSolutionElec(NodeStoreSol<Double> * solhelp){
     solElec_= solhelp;
     sol2_=solhelp;
   }
   
    void SetSolutionMech(NodeStoreSol<Double> * solhelp){
      solMech_= solhelp;
      sol1_=solhelp;
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

    /// scalar electric potential of all nodes of actual element
    Vector<Double> elemPot_;
    Matrix<Double> elemDispl_;
    
    GradientFieldOp<Double> * EfieldOp_;
    
    // local copy of entity iterator
    EntityIterator ent1_;
    
    UInt numDofsA_;
    UInt numDofsB_;
    
    UInt matDimRow_;
    UInt matDimCol_;
    
    NodeStoreSol<Double> * solElec_;
    NodeStoreSol<Double> * solMech_;
    
    BaseMaterial* matDataMech_;
    BaseMaterial* matDataElec_;
    
    bool isHysteresis_; 
    UInt dirP_;   //< direction of polarization
    Double Psat_; //< maximum value of saturation
  };

}

#endif
