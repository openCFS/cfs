// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MICROPIEZORHS
#define FILE_MICROPIEZORHS

#include "Forms/linearForm.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"

namespace CoupledField {
class BaseMaterial;
class EntityIterator;
class EqnMap;
class Grid;
class PiezoCoupling;
class StdPDE;
class linPiezoCoupling;
struct ResultInfo;
template <class TYPE> class GradientFieldOp;
template <class TYPE> class MechStressStrain;
template <class TYPE> class Vector;
template <typename T> class NodeStoreSol;
}  // namespace CoupledField
#ifndef INTEGLIB
#endif

namespace CoupledField
{
  // forward class declaration


  
  // =============================================================================
  // piezoelectric polarization
  // =============================================================================
  
  /// base class for calculation of right hand side due to micro-model for
  /// polarization P
  class MicroPiezoPolarizationRhsInt : public LinearForm
  {
  public:

    ///
    MicroPiezoPolarizationRhsInt( BaseMaterial* matDataCouple, 
                             BaseMaterial* matDataMech, 
                             BaseMaterial* matDataElec, 
                             SubTensorType type); 
    
    ///
    virtual ~MicroPiezoPolarizationRhsInt();
    
    /// Calculation of vector of right hand side 
    virtual void CalcElemVector( Vector<Double> & result,
				 EntityIterator& ent ) {
      EXCEPTION("CalcElemVector not implemented in base class");
    }
    
    //! set objects for computation of E-field
    void Set4NonLinMaterial( Grid* ptGrid, 
                             StdPDE* ptPDE2Elec,
                             shared_ptr<EqnMap> eqnMapElec,
                             shared_ptr<ResultInfo> resultElec,
                             PiezoCoupling* ptPiezoCoupling);

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

  protected:
    //!
    PiezoCoupling* piezoCoupling_;

    //! object for computing electric field
    GradientFieldOp<Double>* EfieldOp_;
    GradientFieldOp<Double>* EfieldPrevOp_;

    //! object for computing mechanical strain
    MechStressStrain<Double>* mechStrainOp_;

    //!object for B-operators
    linPiezoCoupling* piezoBilinearForm_;

    //! pointer to material data
    BaseMaterial* matDataCouple_;
    BaseMaterial* matDataMech_;
    BaseMaterial* matDataElec_;

    NodeStoreSol<Double> * solElec_;
    NodeStoreSol<Double> * solMech_;
    NodeStoreSol<Double> * solPrevElec_;
    NodeStoreSol<Double> * solPrevMech_;

  };

  
  /// class for calculation of right hand side due to a micro-model for
  /// polarization P: RHS of electric equation
  class MicroPiezoPolarizationElecRhsInt : public MicroPiezoPolarizationRhsInt  {

  public:

    MicroPiezoPolarizationElecRhsInt( BaseMaterial* matDataElec,
				 BaseMaterial* matDataPiezo,
				 BaseMaterial* matDataMech,
				 SubTensorType type); 
    
    ///
    virtual ~MicroPiezoPolarizationElecRhsInt();
    
    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
    
  private:

  };

  
  
  /// class for calculation of right hand side due to a micro-model for
  /// polarization P: RHS of mechanical equation
  class MicroPiezoPolarizationMechRhsInt : public MicroPiezoPolarizationRhsInt
  {
  public:

    ///
    MicroPiezoPolarizationMechRhsInt( BaseMaterial* matDataCouple, 
                                 BaseMaterial* matDataMech, 
                                 BaseMaterial* matDataElec, 
                                 SubTensorType type); 
        
    ///
    virtual ~MicroPiezoPolarizationMechRhsInt();
    
    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
    
  private:

  };




} // end of namespace

#endif // FILE_LINEARFORM
