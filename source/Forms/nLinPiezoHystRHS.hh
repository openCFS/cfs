// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PIEZOHYSTRHS
#define FILE_PIEZOHYSTRHS

#include "linearForm.hh"
#include "nLinElastInt.hh"
#include "Utils/ApproxData.hh"
#include "gradfieldop.hh"
#include "linPiezoCoupling.hh"
#include "curlCurlNodeInt.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField
{
  // forward class declaration
  class CoordSystem;
  class PiezoCoupling;
  template<typename> class MechStressStrain;

  
  // =============================================================================
  // piezoelectric polarization
  // =============================================================================
  
  /// base class for calculation of right hand side due to a hysteresis model with
  /// polarization P
  class PiezoPolarizationRhsInt : public LinearForm
  {
  public:

    ///
    PiezoPolarizationRhsInt( BaseMaterial* matDataCouple, 
                             BaseMaterial* matDataMech, 
                             BaseMaterial* matDataElec, 
                             SubTensorType type); 
    
    ///
    virtual ~PiezoPolarizationRhsInt();
    
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

  
  /// class for calculation of right hand side due to a hysteresis model with
  /// polarization P: RHS of electric equation
  class PiezoPolarizationElecRhsInt : public PiezoPolarizationRhsInt  {

  public:

    PiezoPolarizationElecRhsInt( BaseMaterial* matDataElec,
				 BaseMaterial* matDataPiezo,
				 BaseMaterial* matDataMech,
				 SubTensorType type); 
    
    ///
    virtual ~PiezoPolarizationElecRhsInt();
    
    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
    
  private:

  };

  
  
  /// class for calculation of right hand side due to a hysteresis model with
  /// polarization P: RHS of mechanical equation
  class PiezoPolarizationMechRhsInt : public PiezoPolarizationRhsInt
  {
  public:

    ///
    PiezoPolarizationMechRhsInt( BaseMaterial* matDataCouple, 
                                 BaseMaterial* matDataMech, 
                                 BaseMaterial* matDataElec, 
                                 SubTensorType type); 
        
    ///
    virtual ~PiezoPolarizationMechRhsInt();
    
    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
    
  private:

  };




} // end of namespace

#endif // FILE_LINEARFORM
