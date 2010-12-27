// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
      \param algebraicsystem pointer to algebraic system used by PDE

      \param apdeId Id of PDE who called NewmarkFracDamp
      \param ptEQN
      \param aptgrid
      \param aptStdPDE pointer of class from which NewmarkFracDamp is initiated
      \param asubdomainList list of subdomains
      \param adampingList list damping description for subdomains
    */
    NewmarkFracDamp( BaseSystem * algebraicsystem,
                     const PdeIdType apdeId,
                     shared_ptr<EqnMap> eqnMap,
                     Grid * aptgrid,
                     StdPDE * aptStdPDE,
                     shared_ptr<ResultInfo> aresult,
                     StdVector<RegionIdType> asubdomainList,
                     std::map<RegionIdType,DampingType> adampingList);
  
    //! destructor
    virtual ~NewmarkFracDamp();
  
    //! initilization
    //! \param rhsSize size of right hand side vector
    void Init( Double dt, UInt rhsSize );
    
    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS();

    //! store current value to memory
    void AdvanceTimestep(Vector<Double>& solnew);

    //! compute parameters for multiplication
    void CalcParameters(Double dt);

    //! get beta coefficient from Newmark time stepping scheme
    Double GetNewmarkBeta()
    { return beta_;};

  private:

    //! get element solution, for assembling RHS in fractional damping model
    void GetElemSolution ( const Vector<Double>& sol, 
                           Vector<Double>& elemsol, 
                           const StdVector<Integer> eqnNumbers );

    //! compute Weights for Gruenwald-Letnikov formula
    void GLWeights(UInt memory, Double y);

    //! compute Weights for Luise Blanks frac diff spline collocation formula
    void BlankWeights(UInt memory, Double y, bool full);

    //! integrate interpolation of past function values in weight vector
    void CompressWeights();

    //! print solMemoryVal_ into string 'msg'
    void PrintSolMemoryVal(std::string & msg);

    //! name of the pde
    std::string pdename_;

    //! algsys identifier of the pde
    PdeIdType pdeId_;

    //@{
    //! integration parameters
    Double alpha_, gamma_, beta_;
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
    shared_ptr<EqnMap> eqnMap_;

    //! pointer to results description object
    shared_ptr<ResultInfo> result_;

    //! time step
    UInt actStep_;
   
    //! number of terms over which BDF is calculated
    UInt numValues_;

    //! number of truely stored values
    UInt numTrueValues_;

    //! damping type for all regions
    std::map<RegionIdType,DampingType> dampingList_;

    //! list of subdomains
    StdVector<RegionIdType> subdoms_;

    //! flag indicating axisymmetric model
    bool isaxi_;

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

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class NewmarkFracDamp
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! This concept does permit to have different regions with differnt
  //! damping methods.
  //! Although it would be desirable to eliminate the class NewmarkFracDamp
  //! to simplify direct coupling with acoustics.
  //! To improve the situation, we would need a generalized BaseForm
  //! class, which allows each integrator to store old values of the solution
  //! for a certain region. However, at the moment each integrator derived 
  //! from BaseForm only gets a pointer to a reference element and has 
  //! therefore no real memory of what it has calculated before.
  //! 

#endif
} // end of namespace

#endif // FILE_NEWMARKDAMP
