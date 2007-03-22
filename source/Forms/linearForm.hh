// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINEARFORM_2
#define FILE_LINEARFORM_2

#include "baseForm.hh"
#include "Forms/nLinElastInt.hh"
#include "Utils/ApproxData.hh"
#include <Forms/gradfieldop.hh>
#include "linPiezoCoupling.hh"

#ifndef INTEGLIB
#include "Utils/mathParser/mathParser.hh"
#endif

namespace CoupledField
{
  // forward class declaration
  class CoordSystem;


  /// base class class for calculation right hand side
  class LinearForm : public BaseForm
  {
  public:

    ///
    LinearForm();

    ///
    virtual ~LinearForm();

    /// Calculation of RHS vector for double entries, i.e. transient and static 
    virtual void CalcElemVector( Vector<Double> & result,
                                 EntityIterator& ent );
    
    /// Calculation of RHS vector for double entries, i.e. harmonic 
    virtual void CalcElemVector( Vector<Complex> & result,
                                 EntityIterator& ent );
  };



  // =============================================================================
  // edge integration
  // =============================================================================


  /// class for calculation of right hand side of edge elements
  class LinearEdgeInt : public LinearForm
  {
  public:
    ///
    LinearEdgeInt( const std::string& val, UInt direction, 
                   Vector<Double> * coilMidPoint = NULL);

    ///
    virtual ~LinearEdgeInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

  private:

  
    /// direction of source
    /*! 1: x-direction, 2: y-direction, 3: z-direction
     */
    UInt direction_;

    /// midpoint of coil (needed for circular coils to calculate the current dirction)
    Vector<Double> * coilMidPt_;  
  };




  // =============================================================================
  // volume source integration
  // =============================================================================


  /// class for calculation of right hand side of an volume source
  class VolumeSrcInt : public LinearForm
  {
  public:
    ///
    VolumeSrcInt( const std::string& val,
                  bool isaxi, bool coordUpdate = false );

    ///
    virtual ~VolumeSrcInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

    //! Sets a multiplicative factor for element matrix
    void SetFactor( const std::string& factor );

  private:
    
  };


  // =============================================================================
  // permanent magnets in 2D
  // =============================================================================


  /// class for calculation of right hand side of a permanent magnet
  class MagPerm2DInt : public LinearForm
  {
  public:
    ///
    MagPerm2DInt(Vector<Double> val, Double rel, bool isaxi,
                 bool coordUpdate = false );

    ///
    virtual ~MagPerm2DInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

  private:
    //! magnetization
    Vector<Double> perm_;

    //!reluctivity
    Double reluctivity_;
  };


  // =============================================================================
  // nonlinear magnetics
  // =============================================================================

  /// class for calculation of right-hand-side of nonlinear magnetics
  class nLinMagNode2D_linFormInt : public LinearForm
  {
  public:

    /// constructor
    nLinMagNode2D_linFormInt(ApproxData *nlinFnc, Double startVal, 
                             bool axi=false, bool coordUpdate = false );

    /// destructor
    virtual ~ nLinMagNode2D_linFormInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
  
  protected:
    /// material data
    BaseMaterial* matData_;

    Double startmatVal_;
    ApproxData *nlinFnc_;
  };


  // =============================================================================
  // nonlinear mechanics
  // =============================================================================


  /// class for calculation of right-hand-side of nonlinear mechanics
  class nLinMech_linFormInt : public LinearForm
  {
  public:

    /// constructor
    nLinMech_linFormInt(BaseMaterial* matData, bool axi=false);

    /// destructor
    virtual ~nLinMech_linFormInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

  
  protected:
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 3;};

    /// material data
    BaseMaterial* matData_;

  };

  /// class for calculation of right-hand-side of prestress
  class PreStressLinFormInt : public nLinMech_linFormInt
  {
  public:
    /// constructor
    PreStressLinFormInt(BaseMaterial* matData, 
                        const std::string &  aPreStressVal, 
                        Directions stressDir);
  

    /// destructor
    virtual ~PreStressLinFormInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

  private: 
    ///

    ///
    Directions preStressDir_;
  };



  // =============================================================================
  // add mechanical stress as forces on RHS
  // =============================================================================
  

  /// add mecganical stress as forces on RHS
  class  AddStressRHSInt : public LinearForm
  {
  public:

    /// constructor
    AddStressRHSInt(BaseMaterial* matData,
                    Vector<Double>& addStrainVec,
                    SubTensorType type);
    
    /// destructor
    virtual ~AddStressRHSInt();
    
    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
    
    
  protected:
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 3;};

    /// material data
    BaseMaterial* matData_;

    //! strain vector
    Vector<Double> addStress_;

    //!
    SubTensorType subTensorType_;
  };


  // =============================================================================
  // calculation for right hand side of flownoise problem
  // =============================================================================

  /// class for right hand side of flownoise problem
  class LinearFlowNoiseInt : public LinearForm
  {
  public:
    ///
    LinearFlowNoiseInt(BaseFE * aptelem);

    ///
    virtual ~LinearFlowNoiseInt();

    /// Calculation of vector of right hand side for the surface elements on the obstacle (dipole)
    void CalcElemVector4Dip(Matrix<Double>& ptCoord, 
                            const StdVector<UInt> & connecth, 
                            Vector<Double> & Result, 
                            const Vector<Double> gradN_x_P);
    
    /// Calculation of vector of right hand side using a dTijdi vector given by another program 
    void CalcElemVec_withdTijdi(const Matrix<Double>& ptCoord,
                                const Matrix<Double>& dTijdi,
                                Vector<Double> & Result);
    
    /// Calculation of Lighthills tensor diverg. (dTijdi) using velocity values from vortex source
    void CalcElemVec_withVortexVel(const Matrix<Double>& ptCoord,
                                   const Matrix<Double> & NodalVel,
                                   Vector<Double> & Result);
    
    /// Calculation of vector of right hand side given from quadrupole contribution
    void CalcElemVector4Quad(Matrix<Double>& ptCoord, 
                             const StdVector<UInt> & connecth,
                             const Matrix<Double> & FlowData, 
                             Vector<Double> & Result);
    
    /// Extraction of element velocity values from total flowdata matrix to a matrix (connecth, dim)
    void GetQttiesOfElement(Matrix<Double>& elVec, 
                            const Matrix<Double>& FlowData,
                            const StdVector<UInt>& connecth, 
                            UInt matrixRow);
    
  private:

  };


  // =============================================================================
  // compute nonlinear rhs for acoustics
  // =============================================================================

  /// calculation of RHS in nonlinear acoustics using Kuznetsov's equation
  class nLinKuznetsovRHSInt : public LinearForm
  {
  public:
    ///
    nLinKuznetsovRHSInt( const std::string &  val, bool isaxi);

    ///
    virtual ~nLinKuznetsovRHSInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

    //! the solution is needed for setting up the matrix
    virtual void SetActElemSol(Vector<Double>& asol) 
    { sol_ = asol;};

    //! the first time derivative is needed for setting up the matrix
    virtual void SetActElemSolDeriv1(Vector<Double>& solderiv1) 
    { solderiv1_ = solderiv1;};

    //! the second time derivative is needed for setting up the matrix
    virtual void SetActElemSolDeriv2(Vector<Double>& solderiv2) 
    { solderiv2_ = solderiv2;};

    //! set multiplicative factor for matrix
    virtual void SetFactor( const std::string& factor ) 
    { factorN1_ = factor;};

    //! set multiplicative factor for matrix
    virtual void SetSecondFactor( const std::string& factor ) 
    { factorN2_ = factor;};

  private:
    Vector<Double> sol_;        //!< solution
    Vector<Double> solderiv1_;  //!< first time derivative of solution
    Vector<Double> solderiv2_;  //!< second time derivative of solution

    std::string factorN1_;      //!< multiplicative value for integrator
    std::string  factorN2_;     //!< multiplicative value for integrator

    /// source factor
    std::string val_;
  };

  /// calculation of RHS in nonlinear acoustics using Westervelt's equation
  class nLinWesterveltRHSInt : public LinearForm
  {
  public:
    ///
    nLinWesterveltRHSInt( const std::string &  val, bool isaxi);

    ///
    virtual ~nLinWesterveltRHSInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

    //! the solution is needed for setting up the matrix
    virtual void SetActElemSol(Vector<Double>& asol) 
    { sol_ = asol;};

    //! the first time derivative is needed for setting up the matrix
    virtual void SetActElemSolDeriv1(Vector<Double>& solderiv1) 
    { solderiv1_ = solderiv1;};

    //! the second time derivative is needed for setting up the matrix
    virtual void SetActElemSolDeriv2(Vector<Double>& solderiv2) 
    { solderiv2_ = solderiv2;};
    
    virtual void SetFactor( const std::string& factor);

  private:

    Vector<Double> sol_;        //!< solution
    Vector<Double> solderiv1_;  //!< first time derivative of solution
    Vector<Double> solderiv2_;  //!< second time derivative of solution

    std::string  val_;          //!< source factor
  };



  /// class for calculation of right hand side for recovery procedure
  class RHSForRecoveryProcedure : public LinearForm
  {
  public:
    ///
    RHSForRecoveryProcedure();

    ///
    virtual ~RHSForRecoveryProcedure();

    // Calculation of vector of right hand side
    /*
      /\f[
      \int \phi \frac{\partial u^{FEM}}{\partial x_i}
      \f]
      \param ptCoord (input) Matrix with coordinates of the element
      \param fncNodesElem (input) value of solution at nodes of element
      \param aComponent (input) number of the gradient's component 
      \param elemVec (output) vector with result
    */
    void CalcElemVectorRHSForSPR(BaseFE * aptelem,
                                 Matrix<Double>& ptCoord,
                                 Vector<Double> & fncNodesElem,
                                 const UInt aComponent,
                                 Vector<Double> & elemVec);

  private:
    
  };

  // =============================================================================
  // electric polarization
  // =============================================================================

  
  /// class for calculation of right hand side due to a hysteresis model with
  /// polarization P
  class ElecPolarizationInt : public LinearForm
  {
  public:
    ///
    ElecPolarizationInt(bool isaxi, bool coordUpdate = false );
    
    ///
    virtual ~ElecPolarizationInt();
    
    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
    
    ///
    void SetSrcVec(Vector<Double> vec) {D_ = vec;};
  private:
    
    //!value of polarization
    Vector<Double> D_;
  };
  
  
  // =============================================================================
  // piezoelectric polarization
  // =============================================================================
  
  /// base class for calculation of right hand side due to a hysteresis model with
  /// polarization P
  class PiezoPolarizationRhsInt : public LinearForm
  {
  public:

    ///
    PiezoPolarizationRhsInt( BaseMaterial* matDataElec, SubTensorType type ); 
    
    ///
    virtual ~PiezoPolarizationRhsInt();
    
    /// Calculation of vector of right hand side 
    virtual void CalcElemVector( Vector<Double> & result,
				 EntityIterator& ent ) {
      Error("CalcElemVector not implemented in base class",__FILE__,__LINE__);
    }
    
    //! set objects for computation of E-field
    void Set4NonLinMaterial(Grid* ptGrid, 
			    StdPDE* ptPDE,
			    shared_ptr<EqnMap> eqnMap,
			    shared_ptr<ResultInfo> result);

    //! compute normalized irreversible strain
    void ComputeNormalizedSirr(SubTensorType type, UInt dirP,
			       Vector<Double>& Sirr);

  protected:

    //! object for computing electric field
    GradientFieldOp<Double>* EfieldOp_;

    //!object for B-operators
    linPiezoCoupling* piezoBilinearForm_;

    //! vector for normalized irreversible strains
    Vector<Double> baseVecSirr_;

    //! direction of polarization
    UInt dirP_;   //< direction of polarization

    //! pointer of electric material data
    BaseMaterial* matDataElec_;

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
    BaseMaterial* matDataPiezo_;
    BaseMaterial* matDataMech_;

  };

  
  
  /// class for calculation of right hand side due to a hysteresis model with
  /// polarization P: RHS of mechanical equation
  class PiezoPolarizationMechRhsInt : public PiezoPolarizationRhsInt
  {
  public:

    ///
    PiezoPolarizationMechRhsInt( BaseMaterial* matDataElec,
				 BaseMaterial* matDataMech,
				 SubTensorType type); 
    
    ///
    virtual ~PiezoPolarizationMechRhsInt();
    
    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
    
  private:
    BaseMaterial* matDataMech_;

  };


  // =========================================================================
  // mechanic rhs volume integrator
  // =========================================================================
  
  //! Linear integrator for mechanic volume load sources. This integrator
  //! can also work on different local coordinate systems
  class MechVolForceInt : public LinearForm {
    
  public:
    
    //! Constructor
    MechVolForceInt(UInt numDof, const std::string& phase, bool isaxi);
    
    //! Destructor
    virtual ~MechVolForceInt();
    
    //! Set the volume force vector
    //! \param volForce vector with volume force w.r.t. coordSys
    //! \param coordSys pointer to reference coordinate system
    void SetVolForceVector(StdVector<std::string> & volForce, 
                           const CoordSystem * coordSys,
                           bool isUnit, Double volume);
    
    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Complex> & result,
                         EntityIterator& ent );

  protected:

    //! Helper function for calculating part element vector
    template<class TYPE>
    void CalcPartVector( Vector<TYPE>& elemVec, 
                         Vector<TYPE>& loadVec,
                         EntityIterator& ent );
    
    //! Number of degrees of freedom
    UInt numDofs_;
    
    //! Vector with volume force (local coordinate system)
    StdVector<std::string> locForce_;

    //! Phase of force
    std::string phase_;

    //! Reference coordinate system
    const CoordSystem * coordSys_;

    //! Volume of region
    Double volume_;

    //! Flag if force is unit value
    bool isUnitValue_;
    
  };

} // end of namespace

#endif // FILE_LINEARFORM
