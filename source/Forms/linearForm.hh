// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINEARFORM_2
#define FILE_LINEARFORM_2

#include "baseForm.hh"
#include "nLinElastInt.hh"
#include "Utils/ApproxData.hh"
#include "gradfieldop.hh"
#include "linPiezoCoupling.hh"
#include "curlCurlNodeInt.hh"

#ifndef INTEGLIB
#include "Utils/mathParser/mathParser.hh"
#endif

namespace CoupledField
{
  // forward class declaration
  class CoordSystem;
  class linElastInt;


  /// base class class for calculation right hand side
  class LinearForm : public BaseForm
  {
  public:

    ///
    LinearForm( BaseMaterial* matData = NULL );

    ///
    virtual ~LinearForm();

    /// Calculation of RHS vector for double entries, i.e. transient and static 
    virtual void CalcElemVector( Vector<Double> & result,
                                 EntityIterator& ent );
    
    /// Calculation of RHS vector for double entries, i.e. harmonic 
    virtual void CalcElemVector( Vector<Complex> & result,
                                 EntityIterator& ent );

  protected:

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
  // permanent magnets in 3D
  // =============================================================================


  /// class for calculation of right hand side of a permanent magnet
  class MagPerm3DInt : public LinearForm
  {
  public:
    ///
    MagPerm3DInt(Vector<Double> val, Double rel, bool isaxi,
                 bool coordUpdate = false );

    ///
    virtual ~MagPerm3DInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

  private:
    //! number of dofs
    UInt nrDofs_;

    //! magnetization
    Vector<Double> perm_;

    //!reluctivity
    Double reluctivity_;

    CurlCurlNode3DInt *curlOp_;
  };


  // =============================================================================
  // nonlinear magnetics
  // =============================================================================

  /// class for calculation of right-hand-side of nonlinear magnetics
  class nLinMagNode2D_linFormInt : public LinearForm
  {
  public:

    /// constructor
    nLinMagNode2D_linFormInt( BaseMaterial* matData, bool axi=false, 
                              bool coordUpdate = false );

    /// destructor
    virtual ~ nLinMagNode2D_linFormInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
  
  protected:

    Double startmatVal_;
    ApproxData *nlinFnc_;
  };


  /// class for calculation of right-hand-side of nonlinear magnetics
  class nLinMagNode3D_linFormInt : public LinearForm
  {
  public:

    /// constructor
    nLinMagNode3D_linFormInt( BaseMaterial* matData, 
                             bool coordUpdate = false );

    /// destructor
    virtual ~ nLinMagNode3D_linFormInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
  
  protected:
    Double startmatVal_;    //!< star value for reluctivity
    ApproxData *nlinFnc_;   //!< pointer to approximation object for BH curve
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
    
    //! linElastInt pointer
    linElastInt *bilinearStiff_;
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
    
    /// Calculation of Lighthill tensor diverg. (dTijdi) using velocity values from the CFS-CFD solver
    void CalcElemVec_withCFSVel(const Matrix<Double>& ptCoord,
                                const Matrix<Double> & NodalVel,
                                Vector<Double> & Result,
                                Double density); 

    /// Calculation of vector of right hand side given from quadrupole contribution
    void CalcElemVector4Quad(Matrix<Double>& ptCoord, 
                             const StdVector<UInt> & connecth,
                             const Matrix<Double> & FlowData, 
                             Vector<Double> & Result);

    /// Calculation of vector of right hand side using nodal dTijdxj values
    void CalcElemVec4Quad(const Matrix<Double>& ptCoord,
                          const  Matrix<Double> & NodaldTijdxj,
                          Vector<Double> & Result);

    /// Calculation of vector of right hand side using nodal velocity values
    void CalcElemVec4QuadwithVel(const Matrix<Double>& ptCoord, const Matrix<Double> & NodalVel,
                                 Vector<Double> & Result, Vector<Double> & nodalLoadDensity,
                                 Vector<Double>& divLHTensor, const Elem* elem, Double density);

    
    /// Extraction of element velocity values from total flowdata matrix to a matrix (connecth, dim)
    void GetQttiesOfElement(Matrix<Double>& elVec, 
                            const Matrix<Double>& FlowData,
                            const StdVector<UInt>& connecth, 
                            UInt matrixRow);
    
  /// Calculation of vector of right hand side using nodal pressure values
  void CalcElemVec4QuadwithPress(const Matrix<Double>& ptCoord, const Matrix<Double> & NodalPress,
                                 Vector<Double> & Result, const Elem* elem);


  //================= Combustion Noise ========================================================

  void CalcElemVec4CombustionTij(const Matrix<Double>& ptCoord,
				 const Matrix<Double>& NodalVel,
				 const Vector<Double>& NodalRho,
				 Vector<Double>& Result, const Elem* elem);

  void CalcElemVec4CombustionVec(const Matrix<Double>& ptCoord,
				 const Matrix<Double>& NodalVec,
				 Vector<Double>& Result, const Elem* elem);

  void CalcElemVec4CombustionScalar(const Matrix<Double>& ptCoord,
				    const Vector<Double>& NodalScalar,
				    Vector<Double>& Result, const Elem* elem);

  void CalcElemVec4CombustionTijOnSurface(const Matrix<Double>& ptCoord,
					  const Matrix<Double>& NodalVel,
					  const Vector<Double>& NodalRho,
					  Vector<Double>& Result, const Elem* elem);

  void CalcElemVec4CombustionVectorOnSurface(const Matrix<Double>& ptCoord,
					     const Matrix<Double>& NodalVel,
					     Vector<Double>& Result, const Elem* elem);

  void ComputeDerivOfTij( const Double RhoAtIP, const Vector<Double>& RhoDerAtIP,
			  const Vector<Double>& VelAtIP, const Matrix<Double>& VelDerAtIP,
			  const Integer dim, Vector<Double>& helpVec);

  void ComputeNormalVec( const Matrix<Double>& ptCoord, Vector<Double>& nVec);

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
    nLinKuznetsovRHSInt( bool isaxi);

    ///
    virtual ~nLinKuznetsovRHSInt();

    /// Calculation of RHS vector 
    void CalcElemVector( Vector<Double>& result,
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
    virtual void SetFactor( const Double& factor ) 
    { factorN1_ = factor;};

    //! set multiplicative factor for matrix
    virtual void SetSecondFactor( const Double& factor ) 
    { factorN2_ = factor;};

  private:
    Vector<Double> sol_;        //!< solution
    Vector<Double> solderiv1_;  //!< first time derivative of solution
    Vector<Double> solderiv2_;  //!< second time derivative of solution

    Matrix<Double> xiDx_; //!< derivative of shape fncts after global coordinates
    Vector<Double> ShpFncAtIp_; //!< shape fnct at int point
    Vector<Double> CoordAtIp_;
    Vector<Double> solGradAtIp_, solDeriv1GradAtIp_; //!< gradient of sol at int point
    Double solDeriv1AtIp_, solDeriv2AtIp_; //!< 1st + 2nd deriv at int point

    Double factorN1_;      //!< multiplicative value for integrator
    Double factorN2_;     //!< multiplicative value for integrator
  };

  /// calculation of RHS in nonlinear acoustics using Westervelt's equation
  class nLinWesterveltRHSInt : public LinearForm
  {
  public:
    ///
    nLinWesterveltRHSInt( bool isaxi );

    ///
    virtual ~nLinWesterveltRHSInt();

    /// Calculation of RHS vector 
    void CalcElemVector( Vector<Double>& result,
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
    
    virtual void SetFactor( const Double& factor)
    { factor_ = factor;};

  private:

    Vector<Double> sol_;        //!< solution
    Vector<Double> solderiv1_;  //!< first time derivative of solution
    Vector<Double> solderiv2_;  //!< second time derivative of solution

    Vector<Double> ShpFncAtIp_; //!< shape fnct at int point
    Vector<Double> CoordAtIp_;
    Double solAtIp_, solDeriv1AtIp_, solDeriv2AtIp_; //<! sol at integration point

    Double factor_;      //!< multiplicative value for integrator
  };
  
  // =============================================================================
  // acoustic source term in heat conduction PDE
  // =============================================================================
  
  /// compute acoustic power source in heat conduction PDE
  class linAcouPowerSourceInt : public LinearForm
  {
  public:
    ///
    linAcouPowerSourceInt( bool isaxi, bool isharmonic,
                           const std::string& readerId, const std::string& regionName,
                           Double density);
  
    ///
    virtual ~linAcouPowerSourceInt();
    
    //! calculate RHS source vector
    void CalcElemVector( Vector<Double>& elemPower, EntityIterator& ent) {
      
      if( isharmonic_ == true )
        CalcElemVectorComplex( elemPower, ent );
      else
        CalcElemVectorDouble( elemPower, ent );
    }
    
  private:
    
    //! transient acoustic - transient heat conduction coupling
    void CalcElemVectorDouble(Vector<Double>& elemPower, EntityIterator& ent );
    
    //! harmonic acoustic - transient heat conduction coupling
    void CalcElemVectorComplex( Vector<Double>& elemPower, EntityIterator& ent);

    std::string readerId_;
    StdVector<std::string> regionNames_;
    bool isharmonic_;
    Double density_;
    UInt actStep_, lastStep_;
    UInt sequenceStep_;
    shared_ptr<NodeStoreSol<Complex> > res_C_, resD1_C_;
    shared_ptr<NodeStoreSol<Double> > res_R_, resD1_R_;
    
    Vector<Double> ShpFncAtIp_; //! shape function at integration point
    Matrix<Double> xiDx_; //! derivation of shape functions after global coordinates 
    Vector<Double> CoordAtIp_; //! coordinate at integration point
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
  
  
  // =========================================================================
  // rhs volume integrator
  // =========================================================================
  
  //! Linear integrator for volume load sources. This integrator
  //! can also work on different local coordinate systems
  class VolForceInt : public LinearForm {
    
  public:
    
    //! Constructor
    VolForceInt(UInt numDof, const std::string& phase, bool isaxi);
    
    //! Destructor
    virtual ~VolForceInt();
    
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
  
  // =============================================================================
   // edge integration
   // =============================================================================

   /// class for calculation of right hand side of edge elements
    class LinearEdgeSrcInt : public VolForceInt {

    public:
      ///
      LinearEdgeSrcInt( UInt numDof, const std::string& phase, bool isaxi);

      ///
      virtual ~LinearEdgeSrcInt();

      //! Calculation of vector of right hand side 
      void CalcElemVector( Vector<Double> & result,
                           EntityIterator& ent );

      //! Calculation of vector of right hand side 
      void CalcElemVector( Vector<Complex> & result,
                           EntityIterator& ent );

      /// Calculation of vector of right hand side 
      template<class TYPE>
      void CalcPartVector( Vector<TYPE>& elemVec, 
                           Vector<TYPE>& loadVec,
                           EntityIterator& ent );

    private:

    };
    
    /** add mechanical strain as forces on RHS 
     *  needed for homogenization
     *  
     *  implementation is based on Bendsoe/Sigmund: Topology Optimization, p. 122ff */
    class  AddStrainRHSInt : public LinearForm
    {
    public:
      /// constructor
      AddStrainRHSInt(BaseMaterial* matData,
                      const Vector<double>& addStrainVec,
                      SubTensorType type);
      
      /// destructor
      virtual ~AddStrainRHSInt();
      
      /// Calculation of vector of right hand side 
      void CalcElemVector(Vector<double> & result, EntityIterator& ent);
      
      
    protected:
      /// returns nr. of degrees of freedom
      virtual UInt getNrDofs() { return 3; }

      /// material data
      BaseMaterial* matData_;

      //! strain vector
      Vector<Double> addStrain_;

      //! tensor type (plane_strain, full)
      SubTensorType subTensorType_;
      
      //! linElastInt pointer
      linElastInt *bilinearStiff_;
    };

} // end of namespace

#endif // FILE_LINEARFORM
