#ifndef FILE_NEWMARKFRACDAMP_2001
#define FILE_NEWMARKFRACDAMP_2001

#include <General/environment.hh>
#include <Domain/grid.hh>
#include "StdPDE.hh"

#include "timestepping.hh"

namespace CoupledField {

  //! class for time stepping of hyperbolic PDE: method is NewmarkFracDamp

  class NewmarkFracDamp: public TimeStepping
  {
  public:
    //! constructor
    /*!
      \param apdename name of PDE
      \param algebraicsystem pointer to algebraic system used by PDE
      \param ptEQN
      \param aptgrid
      \param aptStdPDE pointer of class from which NewmarkFracDamp is initiated
      \param asubdomainList list of subdomains
      \param adampingList list damping description for subdomains
      \param afracMemory number of stored function values
      \param ainType descriptor for interpolation of past function values
      \param isaxi axisymmetric setup
    */
    NewmarkFracDamp( BaseSystem * algebraicsystem,
                     UInt rhsSize,
                     const PdeIdType apdeId,
                     NodeEQN * ptEQN, 
                     Grid * aptgrid, StdPDE * aptStdPDE, 
                     StdVector<RegionIdType> asubdomainList,
                     StdVector<DampingType> adampingList,
                     UInt afracMemory, InterpolType ainType, 
                     Boolean isaxi );
  
    //! deconstructor
    virtual ~NewmarkFracDamp();
  
    //! initilization
    void Init( std::map<FEMatrixType,Double> & matrix_factors, 
               Double dt );
    
    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS();

    //! compute parameters for multiplication
    void CalcParameters(Double dt);

    //! get beta coefficient from Newmark time stepping scheme
    Double GetNewmarkBeta()
    { return beta_;};

  private:

    //! get element solution, for assembling RHS in fractional damping model
    void GetElemSolution (const Vector<Double>& sol, 
                          Vector<Double>& elemsol, 
                          const StdVector<Integer> & connectPDE);

    //! compute Weights for Gruenwald-Letnikov formula
    void GLWeights(UInt memory, Double y);

    //! compute Weights for Luise Blanks frac diff spline collocation formula
    void BlankWeights(UInt memory, Double y, Boolean full);

    //! print solMemoryVal_ in .info file
    void PrintSolMemoryVal();

    //! name of the pde
    std::string pdename_;

    //! algsys identifier of the pde
    PdeIdType pdeId_;

    //@{
    //! integration parameters
    Double alpha_, gamma_, beta_;
    //@}

    //@{
    //! coefficients from NewmarkFracDamp method
    Double a0_,a1_,a2_,a3_,a4_;
    //@}
    
    //! predictor for nodal solution
    Vector<Double> solpred_;

    //!predictor for derivative of solution
    Vector<Double>  solderiv1pred_;

    //! pointer to grid
    Grid * ptgrid_;
    
    //! pointer to pde
    StdPDE * ptStdPDE_;

    //! pointer to equation object
    NodeEQN * ptEQN_;

    //! last calculated time step
    UInt laststepcalc_;
   
    //! number of timesteps with which frac deriv is calculated
    UInt calclimit_;

    //! describes used damping model for whole domain
    //DampingType dampType_; 

    //! damping type for all regions
    StdVector<DampingType> dampingList_;

    //! all names of subdomains
    StdVector<RegionIdType> subdoms_;

    //! flag indicating axisymmetric model
    Boolean isaxi_;

    //@{ \name Fractional Damping Model

    //! weights of BDF formula
    std::vector<Double> coeff_;

    //! number of stored solution values
    UInt fracMemory_;

    //! storing of solution values
    Vector<Double> *solMemory_; 

    //! describes storing in solmemory_
    std::vector<InterpolType> solMemoryVal_; 

    //! type of interpolation of solution values used
    InterpolType inType_;
    //@}

  };

} // end of namespace

#endif // FILE_NEWMARKDAMP
