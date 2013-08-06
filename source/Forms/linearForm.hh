// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINEARFORM_2
#define FILE_LINEARFORM_2

#include <stddef.h>
#include <string>

#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Utils/StdVector.hh"
#include "Utils/result.hh"
#include "baseForm.hh"

namespace CoupledField {
class ApproxData;
class BaseFE;
class BaseMaterial;
class CurlCurlNode3DInt;
class CurlCurlEdgeInt;
class EntityIterator;
class linGradBDBInt;
struct Elem;
template <class TYPE> class NodeStoreSol;
}  // namespace CoupledField

#ifndef INTEGLIB
#endif

namespace CoupledField
{
  class BiotSavart;
  // forward class declaration
  class CoordSystem;
  class LinMagStrictInt;
  class linElastInt;


  /// base class class for calculation right hand side
  class LinearForm : public BaseForm
  {
  public:

    typedef std::map<UInt, UInt > UintMap;

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


  //! Class for permanent magnets in Edge formulation class MagPermEdgeInt : public LinearForm {
  class MagPermEdgeInt : public LinearForm {

   public:
     //! Constructor
    MagPermEdgeInt( Vector<Double> vecVal, Double rel,
                    bool isaxi, bool coordUpdate = 0 );

     //! Destructor
     virtual ~MagPermEdgeInt();
     
    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

   private:

    //! magnetization
    Vector<Double> perm_;

    //!reluctivity
    Double reluctivity_;

     //! Bilinearform for curl calculation
     CurlCurlEdgeInt* op_;
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
  
  /// class for calculation of right-hand-side of nonlinear magnetics using edge elements
  class nLinMagEdge_linFormInt : public LinearForm
  {
  public:

    /// constructor
    nLinMagEdge_linFormInt( BaseMaterial* matData, 
                            bool coordUpdate = false );

    /// destructor
    virtual ~ nLinMagEdge_linFormInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

  protected:
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

    /// Calculation of vector of right hand side for aeroacouSource PDE in fespace branch
    /// i.e. for extraction of fluidmech pressure
    void CalcElemVec4AeroAcouSrc(const Matrix<Double>& ptCoord,
                                 const Matrix<Double> & NodalVel,
                                 const Matrix<Double> & NodalMeanVel,
                                 Vector<Double> & Result,
                                 const Elem* elem);

    void CalcElemVecSurfForce(const Matrix<Double> & coordMat,
                          const Matrix<Double>& NodalForce,
                          Vector<Double>& elemvec,
                          const Elem* elem);

    /// Calculation of vector of right hand side using nodal velocity values at
    /// the centre of an element
    void CalcElemVec4QuadwithVelCentre(const Matrix<Double>& ptCoord, 
                                       const Matrix<Double> & NodalVel,
                                       Vector<Double> & Result, 
                                       Vector<Double> & nodalLoadDensity,
                                       Vector<Double>& divLHTensor, 
                                       const Elem* elem, Double density);


    /// Calculation of vector of right hand side using nodal velocity values
    void CalcElemVec4QuadwithVel(const Matrix<Double>& ptCoord, 
                                 const Matrix<Double> & NodalVel,
                                 Vector<Double> & Result, 
                                 Vector<Double> & nodalLoadDensity,
                                 Vector<Double>& divLHTensor, 
                                 const Elem* elem, Double density);
    
    /// Calculation of vector of right hand side using nodal velocity values
    void CalcElemVec4QuadwithDivTij(const Matrix<Double>& ptCoord, 
                                 const Matrix<Double> & NodalDivTij,
                                 Vector<Double> & Result, 
                                 Vector<Double> & nodalLoadDensity,
                                 Vector<Double>& divLHTensor, 
                                 const Elem* elem, Double density);

    void CalcLighthillSurfaceTermVel(const Elem* volElem,
                                 const Elem* surfElem,
                                 const Matrix<Double>& ptVolCoord,
                                 const Matrix<Double>& ptSurfCoord,
                                 const Matrix<Double> & volumeVel,
                                 Vector<Double> & surfNormal,
                                 Double density,
                                 Vector<Double> & Result,
                                 Vector<Double> & ResultLHTens);


    void CalcLighthillSurfaceTermVelCenter(const Elem* VolElem,
                                     const Elem* surfElem,
                                     const Matrix<Double>& ptVolCoord,
                                     const Matrix<Double>& ptSurfCoord,
                                     const Matrix<Double> & volumeVel,
                                     Vector<Double> & surfNormal,
                                     Double density,
                                     Vector<Double> & Result,
                                     Vector<Double> & ResultLHTens){
      Exception("CalcLighthillSurfaceTermCenter: not implemented!");
    }

    
    /// Extraction of element velocity values from total flowdata matrix to a matrix (connecth, dim)
    void GetQttiesOfElement(Matrix<Double>& elVec, 
                            const Matrix<Double>& FlowData,
                            const StdVector<UInt>& connecth, 
                            UInt matrixRow);
    
  /// Calculation of vector of right hand side using nodal pressure values
  void CalcElemVec4QuadwithPress(const Matrix<Double>& ptCoord, const Matrix<Double> & NodalPress,
                                 Vector<Double> & Result, const Elem* elem);

  /// Calculation of vector of right hand side for wave equation using laplacian pressure values
  void CalcElemVecLHwithPress(const Matrix<Double>& ptCoord, const Vector<Double> & NodalPress,
                              Vector<Double> & Result, Vector<Double>& nodalLoadDensity,const Elem* elem);

  /// Calculation of vector of right hand side for wave equation for given Lapl. Press (basically a mass integrator)
  void CalcElemVecWaveWithPressD2(const Matrix<Double>& ptCoord, const Vector<Double> & NodalPress,
                                      Vector<Double> & Result, Vector<Double>& nodalLoadDensity,const Elem* elem);

  /// Calculation of vector of right hand side using nodal mean pressure values
  /// Computes the integration of the total differential
  void CalcElemVec4QuadwithPress(const Matrix<Double>& ptCoord,
                                 const Vector<Double> & NodalPress,
                                 const Vector<Double> & NodalPresTDeriv,
                                 const Matrix<Double> & NodalMeanVelocity,
                                 Vector<Double> & Result, const Elem* elem,
                                 Double density);

  ///Calcualte aeroacoustic source term based on lamb vector
  void CalcElemVecWithLamb(const Matrix<Double>& ptCoord,
                           const Matrix<Double> & NodalVelocity,
                           const Matrix<Double> & NodalMeanVelocity,
                           Vector<Double> & Result, Vector<Double> & elemLambvec,
                           const Elem* elem,
                           Double density);

  void CalcLambSurfaceTermVel(const Elem* volElem,
                              const Elem* surfElem,
                              const Matrix<Double>& ptVolCoord,
                              const Matrix<Double>& ptSurfCoord,
                              const Matrix<Double> & volumeVel,
                              const Matrix<Double> & volumeMeanVel,
                              Vector<Double> & surfNormal,
                              Double density,
                              Vector<Double> & Result);


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
  // acoustic source term in CAA: flow pressure
  // =============================================================================
  
  /// compute acoustic power source in heat conduction PDE
  class linFlowPressureRHSInt : public LinearForm
  {
  public:
    ///
    linFlowPressureRHSInt( bool isaxi, 
                        const std::string& readerId, 
                        const std::string& regionName );
  
    ///
    virtual ~linFlowPressureRHSInt();
    
    //! calculate RHS source vector
    void CalcElemVector( Vector<Double>& rhsvec, EntityIterator& ent);
      

  private:
    
    std::string readerId_;
    StdVector<std::string> regionNames_;
    UInt actStep_, lastStep_;
    UInt sequenceStep_;
    shared_ptr<NodeStoreSol<Double> > resPress_;
    
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
  // rhs maxwell homogenization integrator
  // =========================================================================

  //! Linear integrator for volume charge sources for maxwell homogenization.
  //! This integrator can also work on different local coordinate systems
  class VolChargeHomInt : public LinearForm {

  public:

    //! Constructor
    VolChargeHomInt(BaseMaterial* matData, Global::ComplexPart matDataType, const std::string& phase, bool isaxi);

    //! Destructor
    virtual ~VolChargeHomInt();

    //! Set the volume force vector
    //! \param volForce vector with volume force w.r.t. coordSys
    //! \param coordSys pointer to reference coordinate system
    void SetVolChargeVector(Vector<double> & volCharges, const CoordSystem * coordSys, bool isUnit, Double volume, const int dim);

    //! Calculation of vector of right hand side
    //   template <typename T>
    void CalcElemVector( Vector<double> & result,
        EntityIterator& ent );

    //! Calculation of vector of right hand side
    void CalcElemVector( Vector<Complex> & result,
        EntityIterator& ent );

    linGradBDBInt* GetBDBInt() {return bilinearStiff_;};

  protected:

    //! Helper function for calculating part element vector
    /*    template<class TYPE>
      void CalcPartVector( Vector<TYPE>& elemVec,
                           Vector<TYPE>& chargeVec,
                           EntityIterator& ent );*/

    //! Number of degrees of freedom
    UInt numDofs_;

    //! Vector with volume charge (local coordinate system)
    Vector<double> locCharges;

    //! Phase of force
    std::string phase_;

    Global::ComplexPart locMatDataType_;

    //! Reference coordinate system
    const CoordSystem * coordSys_;

    //! Volume of region
    Double volume_;

    //! Dimension of Domain
    UInt dim_;

    //! Flag if force is unit value
    bool isUnitValue_;

    //! Material Data pointer
    BaseMaterial* matData_;

    //! linGradBDBInt pointer
    linGradBDBInt *bilinearStiff_;

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
      
      /** overwrites the virtual base function */
      void CalcElemVector(Vector<double> & result, EntityIterator& ent) { CalcElemVector(result, ent, NULL); }

      /** Calculation of the RHS.
       * @param test_strain optional, overwrites the value addStrainVec set in the constructor */
      void CalcElemVector(Vector<double> & result, EntityIterator& ent, const Vector<Double>* test_strain = NULL);
      
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
    
    // =============================================================================
    // BiotSavartSourceInt (RHS integrator for Biot-Savart coil of magnetic Scalar PDE)
    // =============================================================================

    //! Calculates the RHS for Biot-Savart coils for scalar potential formulation
    
    //! This class implements the RHS excitation for a magnetic-scalar formulation
    //! of the magnetic field. The integrator calcualates
    //! \f [ \int \nabla N_a (\mu - \mu_0) \mathbf{H}_s d\Omega \f ]
    class BiotSavartSourceInt : public LinearForm {
    public:
      
      //! Constructor
      BiotSavartSourceInt( BaseMaterial* matData,
                           SubTensorType tensorType,
                           shared_ptr<BiotSavart>& bs,
                           bool isaxi,
                           bool coordUpdate = false );

      //! Destructor
      virtual ~BiotSavartSourceInt();

      //! Calculation of vector of right hand side 
      void CalcElemVector( Vector<Double> & result,EntityIterator& ent );

      //! Set additional multiplicative factor
      void SetFactor( const std::string& factor ); 
      
    private:
      
      //! Biot-Savart object for coil excitation
      shared_ptr<BiotSavart> biotSavart_;
          
    };
    
    // ============================================================================
    // BiotSavarMechCouplingRHSInt
    // =============================================================================

    //! Calculates the RHS coupling vector of magnetostriction for fundamentaal field

    //! This class implements the RHS coupling vector of the mechanic PDE
    //! to the fundamental field \f [H_s \f ], as generated by the law
    //! of Biot-Savart. The integrator has the form
    //! \f [$ \int ( {\cal B} \mathbf{u}') [\mathbf{e}]^T \mathbf{H}_s d\Omega$ \f ] 
    class BiotSavartMechCouplingInt : public LinearForm {
    public:
      //! Constructor
      BiotSavartMechCouplingInt( BaseMaterial* matDataCouple,
                                 SubTensorType type,
                                 shared_ptr<BiotSavart> & bs,
                                 bool isaxi,
                                 bool coordUpdate = false  );

      //! Destructor
      virtual ~BiotSavartMechCouplingInt();

      //! Calculation of vector of right hand side 
      void CalcElemVector( Vector<Double> & result,EntityIterator& ent );

    private:

      //! Biot-Savart object for coil excitation
      shared_ptr<BiotSavart> biotSavart_;

      //! Magnetostriction bilinearform. Needed for B-operator.
      LinMagStrictInt* magStrictForm_;

    };



    //! add mechanical strain as forces on RHS 
    //! due to magnetostrictive effect  
    class  AddMagStrictStrainRHSInt : public LinearForm
    {
    public:
      //! constructor
      AddMagStrictStrainRHSInt(BaseMaterial* matData,
                               const std::string& readerId,
                               const std::string& regionName,
                               UintMap elemGlobalLocal,
                               SubTensorType type);
      
      //! destructor
      virtual ~AddMagStrictStrainRHSInt();
      
      //! overwrites the virtual base function 
      void CalcElemVector(Vector<double> & result, EntityIterator& ent);

    protected:
      //! returns nr. of degrees of freedom
      virtual UInt getNrDofs() { return 3; }

      //! computes irreversible strain
      void computeStrainIrr(Vector<double>& vecB );

      //! material data
      BaseMaterial* matData_;

      //! map for global to local element number
      UintMap globalToLocalElemeNr_;

      //! strain vector
      Vector<Double> strainIrr_;

      //! tensor type (plane_strain, full)
      SubTensorType subTensorType_;
      
      //! linElastInt pointer
      linElastInt *bilinearStiff_;

      //result of magnetic computation
      shared_ptr<BaseResult> result_;
      std::string readerId_;
      StdVector<std::string> regionNames_;
      UInt actStep_, lastStep_;
      UInt sequenceStep_;
      shared_ptr<NodeStoreSol<Double> > res_R_;
      UInt offsetElemLevel_;
};
    

} // end of namespace

#endif // FILE_LINEARFORM
