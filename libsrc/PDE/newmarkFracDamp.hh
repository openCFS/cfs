#ifndef FILE_NEWMARKFRACDAMP_2001
#define FILE_NEWMARKFRACDAMP_2001

#include <General/environment.hh>
#include <DataInOut/ParamHandling/ConfFile.hh>
#include <Domain/grid.hh>

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
  */
  NewmarkFracDamp (std::string apdename, BaseSystem * algebraicsystem, NodeEQN * ptEQN, 
				   Grid * aptgrid, BasePDE * aptBasePDE, 
				   StdVector<std::string> subdomainList,
				   StdVector<std::string> adampingList,
				   Integer afrac_memory, InterpolType ainType, Boolean isaxi);

  //! deconstructor
  virtual ~NewmarkFracDamp();
  
  //! initilization
  virtual void Init(Double * matrix_factors, Double dt);

  //! perform predictor step
  virtual void Predictor(Vector<Double>& solold);

  //! perform corrector step
  virtual void Corrector(Vector<Double>& solnew);

  //! perform an update to RHS
  virtual void UpdateRHS();

  //! compute parameters for multiplication
  void CalcParameters(Double dt);

  //! get beta coefficient from Newmark time stepping scheme
  Double GetNewmarkBeta()
  { return beta_;};

private:

  //! get element solution, needed for assembling RHS in fractional damping model
  void GetElemSolution (const Vector<Double>& sol, 
						Vector<Double>& elemsol, 
						const StdVector<Integer> & connectPDE);

  //! compute Weights for Gruenwald-Letnikov formula
  void GLWeights(Integer memory, Double y);

  //! compute Weights for Luise Blanks frac diff spline collocation formula
  void BlankWeights(Integer memory, Double y, Boolean full);

  //! print solMemoryVal_ in .info file
  void PrintSolMemoryVal(Integer actSD);

  std::string pdename_;

  Double alpha_, gamma_, beta_;     //!< integration parameters
  Double a0_,a1_,a2_,a3_,a4_;       //!< coefficients from NewmarkFracDamp method

  Vector<Double> solpred_, solderiv1pred_; //!< predictors
  Grid * ptgrid_;
  BasePDE * ptBasePDE_;

  Integer laststepcalc_;  //!< last calculated time step
  Integer calclimit_;     //!< number of timesteps with which frac deriv is calculated

  //DampingType dampType_; //!< describes used damping model for whole domain
  StdVector<std::string> dampingList_; //!< damping type strings for all regions
  StdVector<std::string> subdoms_;     //!< all names of subdomains

  // For fractional damping model
  Double alpha0_, y_;         //!< parameter of damping model
  std::vector<Double> coeff_; //!< weights of BDF formula
  Integer fracMemory_;        //!< number of stored solution values
  Vector<Double> *solMemory_; //!< storing of solution values
  std::vector<InterpolType> solMemoryVal_; //!< describes storing in solmemory_
  InterpolType inType_;       //!< type of interpolation of solution values used


  //
  Boolean isaxi_;
  
};

} // end of namespace

#endif // FILE_NEWMARKDAMP
