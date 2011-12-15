// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWMARKFRACDAMPMECH_2005
#define FILE_NEWMARKFRACDAMPMECH_2005

#include <map>
#include <string>
#include <vector>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "Utils/StdVector.hh"
#include "timestepping.hh"

namespace CoupledField {

  //! class for time stepping of hyperbolic PDE: method is NewmarkFracDampMech

class BaseFE;
class BaseMaterial;
class BaseSystem;
class EqnMap;
class Grid;
class StdPDE;
template <class TYPE> class Matrix;

  class NewmarkFracDampMech: public TimeStepping
  {
  public:
    //! constructor
    /*!
      \param algebraicsystem pointer to algebraic system used by PDE
      \param apdeId Id of PDE who called NewmarkFracDampMech
      \param ptEQN
      \param aptgrid
      \param aptStdPDE pointer of class from which NewmarkFracDampMech is initiated
      \param asubdomainList list of subdomains
      \param adampingList list damping description for subdomains
    */
    NewmarkFracDampMech( BaseSystem * algebraicsystem,
                         const PdeIdType apdeId,
                         shared_ptr<EqnMap> eqnMap,
                         Grid * aptgrid,
                         StdPDE * aptStdPDE, 
                         StdVector<RegionIdType> asubdomainList,
                         std::map<RegionIdType,DampingType> adampingList,
                         PtrParamNode systemNode );
  
    //! deconstructor
    virtual ~NewmarkFracDampMech();
  
    //! initilization
    //! \param rhsSize size of right hand side vector
    void Init( Double dt, UInt rhsSize );
    
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

    void CalcStress(BaseFE * aptelem, BaseMaterial* matDa, StdVector<Integer> connect_PDE, 
                    Matrix<Double> & ptCoord, Vector<Double> & stressVector,
                    Vector<Double> &displacementVector, Integer elemNr);

    void GetStressVector(Vector<Double> &stressVec,UInt elemNr,UInt memory);

    void InsertStressVector(Vector<Double> &stressVec, Integer elemNr, Integer memory);

    void GetAMat(Matrix<Double> & aMat);

    void GetAlphaMat(Matrix<Double> & alphaMat);

    void GetBetaMat(Matrix<Double> & betaMat, Double E, BaseMaterial* matData);

    UInt getStressDim();

    UInt getDim();

    /// returns dimension of D matrix
    virtual UInt getDimD(){return 3;};

    //! name of the pde
    std::string pdename_;

    //! subType of the pde
    std::string subType_;

    //! algsys identifier of the pde
    PdeIdType pdeId_;

    //@{
    //! integration parameters
    Double alpha_, gamma_, beta_;
    //@}

    //@{
    //! coefficients from NewmarkFracDampMech method
    Double a0_,a1_,a2_,a3_,a4_;
    Integer dofs_,numEQNs_;
    //@}
    
    Double dampAlpha_, dampBeta_;
    Double timeStepPowerFracDeriv_;

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

    std::string geomType_;    

    //@{ \name Fractional Damping Model

    //! weights of BDF formula
    std::vector<Double> coeff_;

    //! number of stored solution values
    UInt fracMemory_;

    //! storing of solution values
    Vector<Double> *solMemory_; 

    //! storing of stress values
    Vector<Double> *stressHistoryEl_; 

    //! describes storing in solmemory_
    std::vector<InterpolType> solMemoryVal_; 

    //! type of interpolation of solution values used
    InterpolType inType_;
    //@}
    Integer modulo_;
    Double elastModule_;

    //! node of parameter file of mechanic PDE
    PtrParamNode mechNode_;


  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class NewmarkFracDampMech
  //! 
  //! \purpose 
  //! Considering fractional damping for mechanical fields
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
  //! Although it would be desirable to eliminate the class NewmarkFracDampMech
  //! to simplify direct coupling with acoustics.
  //! To improve the situation, we would need a generalized BaseForm
  //! class, which allows each integrator to store old values of the solution
  //! for a certain region. However, at the moment each integrator derived 
  //! from BaseForm only gets a pointer to a reference element and has 
  //! therefore no real memory of what it has calculated before.
  //! 

#endif
} // end of namespace

#endif // FILE_NEWMARKDAMPMECH_2005
