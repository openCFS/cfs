#ifndef FILE_BASEPARAMETER_CLA
#define FILE_BASEPARAMETER_CLA

namespace CoupledField
{

//! Parameter class for lagebraic equation solvers
class BaseParameter
{
public:
  //! Constructor
  BaseParameter();

  //! Destructor
  virtual ~BaseParameter();

  //! Set all parameters of the solver to predefined values
  void Set();

  ///
  void SetAccuracy(Double aeps)
    {eps = aeps;}

  ///
  void SetMaxNumIter(Integer amaxnumiter)
    {maxnumiter = amaxnumiter;}

  ///
  void SetPrecond(Integer aprecond)
    {precond = aprecond;}

  ///
  void SetSolver(Integer asolver)
    {solver = asolver;}

  ///
  void SetDampIter(Double adampiter)
    {dampiter = adampiter;};

  ///
  void SetEpsMat(Double aepsmat)
    {epsmat = aepsmat;};

  ///
  void SetAlpha(Double aalpha)
    {alpha = aalpha;};

  ///
  void SetCoarseSystem(Integer acoarsesystem)
    {coarsesystem = acoarsesystem;};
  
  ///
  void SetCycle(Integer acycle)
    {cycle = acycle;};

  ///
  void SetSmoothType(Integer asmoothingtype)
    {smoothingtype = asmoothingtype;};

  ///
  void SetSmoothStepFor(Integer asmoothingfor)
    {smoothingfor = asmoothingfor;};

  ///
  void SetSmoothStepBack(Integer asmoothingback)
    {smoothingback = asmoothingback;};

  ///
  void SetSpectralMatrix()
    {spemat = 1;};

  ///
  void SetPenalty(Double apenalty)
    {penalty = apenalty;};

  ///
  void SetOutBackMemory()
    {outback = 1;};

  //! Get the preconditioner-type (IC, ILU, etc.)
  Integer GetPrecond() const {return precond;};

  //! Get the solver type (CG, etc.)
  Integer GetSolver() const {return solver;};

  ///
  Integer GetCoarseSystem() const {return coarsesystem;};

  ///
  Integer GetCoarseSolver() const {return coarsesolver;};

  ///
  Integer GetMaxNumIter() const {return maxnumiter;};

  ///
  Integer GetCoarsening() const {return coarsening;};

  ///
  Integer GetMaxNeighbour() const {return maxneighbour;};

  ///
  Integer GetTransfer() const {return transfer;};

  ///
  Integer GetNorm() const {return norm;};

  ///
  Integer GetCycle() const {return cycle;};

  ///
  Integer GetSmoothingStepFor() const {return smoothingfor;};
  
  ///
  Integer GetSmoothingStepBack() const {return smoothingback;};

  ///
  Integer GetSmoothingType() const {return smoothingtype;};

  ///
  Integer GetAuxiliaryMatrix() const {return auxmat;};

  ///
  Integer GetSpectralMatrix() const {return spemat;};

  ///
  Double GetEps() const {return eps;};

  ///
  Double GetEpsInt() const {return epsint;};

  ///
  Double GetEpsMat() const {return epsmat;};

  ///
  Double GetEpsMach() const {return epsmach;};

  ///
  Double GetDampIter() const {return dampiter;};

  ///
  Double GetDampSmooth() const {return dampsmooth;};

  ///
  Double GetAlpha() const {return alpha;};

  ///
  Double GetBeta() const {return beta;};

  ///
  Double GetPenalty() const {return penalty;};

  ///
  Integer GetOutBackMemory() const {return outback;};
  
private:
  ///
  Double eps, epsint, epsmat, epsmach;
  Double dampiter, dampsmooth;
  Double alpha, beta;
  Double penalty;

  Integer precond, solver;
  Integer coarsesystem, coarsesolver;
  Integer maxnumiter;

  Integer coarsening, maxneighbour, transfer;
  Integer norm, cycle, smoothingfor, smoothingback, smoothingtype, spemat, auxmat;
  Integer outback;
};

}

#endif // FILE_BASEPARAMETER_CLA
