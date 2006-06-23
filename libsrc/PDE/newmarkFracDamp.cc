#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "newmarkFracDamp.hh"
#include "Forms/massInt.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "basePDE.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/mathfunctions.hh"

namespace CoupledField {

  NewmarkFracDamp::NewmarkFracDamp( BaseSystem * algebraicsystem,
                                    const PdeIdType apdeId,
                                    shared_ptr<EqnMap> eqnMap,
                                    Grid * aptgrid,
                                    StdPDE * aptStdPDE,
                                    StdVector<RegionIdType> asubdomainList,
                                    std::map<RegionIdType,DampingType> adampingList) 
    :TimeStepping( algebraicsystem ){
	
    ENTER_FCN( "NewmarkFracDamp::NewmarkFracDamp" );
  
    pdename_     = aptStdPDE->GetName();
    pdeId_       = apdeId;
    ptgrid_      = aptgrid;
    ptStdPDE_    = aptStdPDE;
    eqnMap_      = eqnMap;
	
    subdoms_     = asubdomainList;
    dampingList_ = adampingList;
    fracMemory_  = aptStdPDE->GetFracMemory();

    bool found = false;
    std::map<RegionIdType,DampingType>::iterator it;
    for ( it = dampingList_.begin(); it != dampingList_.end(); it++ ) {
      
      if ( (*it).second == FRACTIONAL_GL ||
           (*it).second == FRACTIONAL_BLANK ) {
        found = true;
        break;
      }
    }
    
    if ( found == true ) {
      inType_ = NOTUSED;
    } else {
      inType_ = LIN1PT;
    }

    isaxi_       = aptStdPDE->GetIsaxi();
  
    alpha_ = 0.0;
    beta_  = 0.25;
    gamma_ = 0.5;

    //check if integration parameters are defined
    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    if(analysis != "paramIdent")
      Info->PrintF( pdename_, "NewmarkFracDamp: Using defaults for alpha, \
beta and gamma!\n" );

 
  }

  NewmarkFracDamp::~NewmarkFracDamp() {

    ENTER_FCN( "NewmarkFracDamp::~NewmarkFracDamp" );
  }

  void NewmarkFracDamp::Init( std::map<FEMatrixType,Double> & matrix_factors,
                              Double dt, UInt rhsSize ) {
	
    ENTER_FCN( "NewmarkFracDamp::Init" );

    rhsSize_ = rhsSize;
    dt_ = dt;
    CalcParameters(dt_);

    matrix_factors[STIFFNESS] = 1.0;
    matrix_factors[DAMPING] = 1.0*a4_; // needed for thermoviscous damping
    matrix_factors[CONVECTION] = 0.0;
    matrix_factors[MASS] = 1.0*a2_;


    // get the memory
    solderiv1_.Resize(rhsSize_);  
    solderiv1_.Init();
    solderiv2_.Resize(rhsSize_);  
    solderiv2_.Init();
  

    solpred_.Resize(rhsSize_); 
    solpred_.Init();
    solderiv1pred_.Resize(rhsSize_); 
    solderiv1pred_.Init();
  
    if ( fracMemory_ <= 1 )
      Error("Attenuation model needs frac_memory larger than 1",
            __FILE__,__LINE__);
    else {
      // get the memory
      solMemory_ = new Vector<Double>[fracMemory_];
      for (UInt i=0; i<fracMemory_; i++) {
        solMemory_[i].Resize(rhsSize_);  
        solMemory_[i].Init();
      }
    }
  }


  void NewmarkFracDamp::Predictor(Vector<Double>& solold) {

    ENTER_FCN( "NewmarkFracDamp::Predictor" );

    actStep_ = ptStdPDE_->GetTimeStepCounter();

    // determine number of terms over which BDF is calculated
    //   assumes first nstep = 1 (see transientdriver.cc)!
    numValues_ = solMemoryVal_.size();

    // determine number of past values stored in solMemory_
    numTrueValues_ = 0;
    for ( UInt i=0; i < numValues_; i++ ) {
      if ( solMemoryVal_[i] == trueVAL )
        numTrueValues_++;
    }

    solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
    solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
  }


  void NewmarkFracDamp::UpdateRHS() {

    ENTER_FCN( "NewmarkFracDamp::UpdateRHS" );

    // mass part
    Vector<Double> coeffMass;

    coeffMass = solpred_*a2_;
    algsys_->UpdateRHS(MASS,coeffMass.GetPointer());

    // damping part
    Matrix<Double>  elemmat,ptCoord;
    StdVector<UInt> connecth;
    Vector<Double>  rhsAssemble, rhsvec, elemsol;
    Double          density, compressibility, c0, alpha0, y, factor;
    std::map<RegionIdType, BaseMaterial*> mymaterialData;
    mymaterialData = ptStdPDE_->getPDEMaterialData();

    for ( UInt actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {

      if ( dampingList_[subdoms_[actSD]] == NONE ) {
        // no damping term has to be computed
      }
      else if ( dampingList_[subdoms_[actSD]] == THERMOVISCOUS || 
                dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
		
        Vector<Double> coeffDamp;
        coeffDamp = -solderiv1pred_ + solpred_*a4_;
        algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
      }
      else {
	mymaterialData[subdoms_[actSD]]->GetScalar(density,DENSITY,REAL);
	mymaterialData[subdoms_[actSD]]->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);
	mymaterialData[subdoms_[actSD]]->GetScalar(alpha0,ACOU_ALPHA,REAL);
	mymaterialData[subdoms_[actSD]]->GetScalar(y,FRACTIONAL_EXPONENT,REAL);

	c0 = sqrt(compressibility/density);

        // factor: pre factor in Newmark timestepping scheme
        //         times pre factor of damping term
        factor = a2_ * density * 2.0 * alpha0 / c0 / sin((y-1.0)*PI/2.0);

        // coeff_: weight factors in BDF
        if ( dampingList_[subdoms_[actSD]] == FRACTIONAL_GL ) {
          factor *= std::exp(-(y-1.0)*std::log(dt_));

          GLWeights(numValues_, y);
        }
        else if ( dampingList_[subdoms_[actSD]] == FRACTIONAL_GL_INT ) {
          factor *= std::exp(-(y-1.0)*std::log(dt_));

          GLWeights(numValues_, y);
          CompressWeights();
        }
        else if ( dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK ) {
          factor *= std::exp(-(y - 1.0)*std::log(dt_)) 
            * std::exp(-gammaln(1.0- (y - 1.0)) );

          BlankWeights(numValues_, y, true);
        }
        else if ( dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK_INT ) {
          factor *= std::exp(-(y - 1.0)*std::log(dt_)) 
            * std::exp(-gammaln(1.0- (y - 1.0)) );

          BlankWeights(numValues_, y, true);
          CompressWeights();
        }

        ElemList actSDList(ptgrid_ );
        actSDList.SetRegion( subdoms_[actSD] );
        
        EntityIterator it = actSDList.GetIterator();

        BaseForm * bilinear_mass = new MassInt(factor, isaxi_);
        for ( it.Begin(); !it.IsEnd(); it++ ) {

          connecth=it.GetElem()->connect;
          
          bilinear_mass->CalcElementMatrix(elemmat,it,it);
          
          // map connect to PDE node numbers
          StdVector<Integer> connect_PDE;
          eqnMap_->GetNodeEqn(connecth, connect_PDE);

#ifdef DEBUG
          UInt elemNum = it.GetElem()->elemNum;
		  
          // output matrix with which BDF is computed
          (*debug) << "Damping   matrix (fractional Damping) of Element " 
                   << elemNum << std::endl;
          (*debug) << elemmat << std::endl;
#endif
          // compute BDF
          GetElemSolution(solpred_, elemsol, connect_PDE);
          rhsvec = -elemsol * coeff_[0];

          for (UInt i=1; i<numTrueValues_; i++) {

            GetElemSolution(solMemory_[i-1], elemsol, connect_PDE);
            rhsvec += elemsol * coeff_[i];
          }

          rhsAssemble = elemmat * rhsvec;

#ifdef DEBUG
          (*debug) << "BDF vector of timestep " << actStep_ << std::endl
                   << rhsvec << std::endl;
#endif
          
          //assemble to RHS
          algsys_->SetElementRHS(&rhsAssemble[0], pdeId_, 
                                 connect_PDE.GetPointer(),
                                 connect_PDE.GetSize());
        }
        delete bilinear_mass;
      }
    }
  }


  void NewmarkFracDamp::Corrector(Vector<Double>& solnew)
  {
    ENTER_FCN( "NewmarkFracDamp::Corrector" );

    solderiv2_ = (solnew - solpred_) * a2_;
    solderiv1_ = solderiv1pred_ + solderiv2_*a3_;
  }


  void NewmarkFracDamp::AdvanceTimestep(Vector<Double>& solnew)
  {
    ENTER_FCN( "NewmarkFracDamp::AdvanceTimestep" );

    // store function values, reverse storing!!!!
    // solMemory_[0]=p_n
    //   solMemory_[1]=p_(n-1)
    //     ...
    //       solMemory_[fracMemory_ - 1]=p_(n-fracMemory_+1)

    // no interpolation
    if (inType_ == NOTUSED) {
	  
      for ( UInt i=fracMemory_-1; i>=1; i-- ) {
        solMemory_[i] = solMemory_[i-1];
      }
      solMemory_[0] = solnew - solpred_;

      if (actStep_ <= fracMemory_)
        solMemoryVal_.insert( solMemoryVal_.begin(),trueVAL);
      // when solMemoryVal_ reaches size fracMemory_ , all entries are trueVAL
      //  and stay the same for all time steps

    }
    // linear interpolation of one solution value 
    else if (inType_ == LIN1PT) {

      if (actStep_ <= fracMemory_) {
        for (UInt i=fracMemory_-1; i>=1; i--)
          solMemory_[i] = solMemory_[i-1];

        solMemoryVal_.insert( solMemoryVal_.begin(),trueVAL);
      }
      else if ( (actStep_ > fracMemory_) 
                && (actStep_ < 2*fracMemory_) ) {
        for (UInt i=(2*fracMemory_-actStep_-1); i>=1; i--)
          solMemory_[i] = solMemory_[i-1];
		
        solMemoryVal_[0] = trueVAL;
        solMemoryVal_.insert( (solMemoryVal_.begin()
                               + 2*fracMemory_-actStep_),LIN1PT);
      }
      else if (actStep_ >= 2*fracMemory_ ) {
        if ( (actStep_%2)==0 ) {
          for (UInt i=fracMemory_-1; i>=1; i--)
            solMemory_[i] = solMemory_[i-1];

          solMemoryVal_.pop_back();
          solMemoryVal_.pop_back();
          solMemoryVal_.insert(solMemoryVal_.begin(), trueVAL);
        }
        else if ( (actStep_%2)==1 )
          solMemoryVal_.insert( (solMemoryVal_.begin()+1), LIN1PT);
      }
      solMemory_[0] = solnew - solpred_;
    }
  }

  void NewmarkFracDamp::CalcParameters(Double dt)
  {
    ENTER_FCN( "NewmarkFracDamp::CalcParameters" );

    //for predictors
    a0_ = 0.5*(1-2.0*beta_)*dt_*dt_;
    a1_ = (1-gamma_)*dt_;

    //for correctors, matrices
    a2_ = 1.0/(beta_*dt_*dt_);
    a3_ = gamma_*dt_;

    //for RHS, matrices
    a4_ = gamma_ / (beta_*dt_);
  }


  void NewmarkFracDamp::
  GetElemSolution ( const Vector<Double>& sol, 
                    Vector<Double>& elemsol,
                    const StdVector<Integer> & connectPDE ) {

    ENTER_FCN( "NewmarkFracDamp::GetElemSolution" );
    elemsol.Resize(connectPDE.GetSize());
    for (UInt eqn=0; eqn<connectPDE.GetSize(); eqn++) 
      elemsol[eqn] = sol[connectPDE[eqn]-1];
  }

  void NewmarkFracDamp::GLWeights(UInt memory, Double y ) {

    ENTER_FCN( "NewmarkFracDamp::GLWeights" );
   
    // reserve memory for weights of BDF, order of derivative is y-1
    coeff_.resize(memory+1);
   
    coeff_[0] = 1.0; // Index 0

    for (UInt i=1; i <= memory; i++) // Index 1 .. memory
      coeff_[i] = coeff_[i-1] * (((Double)i - 1.0) - (y - 1.0))/(Double)i;
  }


  void NewmarkFracDamp::BlankWeights(UInt memory, Double y, bool full) {

    ENTER_FCN( "NewmarkFracDamp::BlankWeights" );
 
    Double pot;
    pot = 1.0- (y - 1.0);

    // reserve memory for weights of BDF, order of derivative is y-1
    coeff_.resize(memory+1);
      
    coeff_[0] = 1.0/ pot; // Index 0
 
    if (full) {
      for (UInt i=1; i<memory; i++) // Index 1 .. (memory-1)
        coeff_[i]= ( std::exp(pot* std::log((Double)i - 1.0))
                     - 2.0 * std::exp(pot* std::log((Double)i))
                     + std::exp(pot* std::log((Double)i + 1.0)) ) / pot;
      // Index memory
      coeff_[memory]= std::exp(-(y-1.0)*std::log((Double)fracMemory_)) + 
        ( std::exp(pot* std::log((Double)fracMemory_ - 1.0))
          - std::exp(pot* std::log((Double)fracMemory_)) ) / pot;
    }
    else {
      for (UInt i=1; i<=memory; i++) // Index 1 .. memory
        coeff_[i]= ( std::exp(pot* std::log((Double)i - 1.0))
                     - 2* std::exp(pot* std::log((Double)i))
                     + std::exp(pot* std::log((Double)i + 1.0)) ) / pot;
    }
  }

  void NewmarkFracDamp::CompressWeights() {

    ENTER_FCN( "NewmarkFracDamp::CompressWeights" );

    std::vector<Double> newCoeff;
    newCoeff.resize(numTrueValues_);
    //for (UInt i=0; i<numTrueValues_; i++) 
    //  newCoeff[i] = 0;
    newCoeff.assign(numTrueValues_+1, 0); // initialize with zeros

    // Index 0
    newCoeff[0] = coeff_[0];

    // Index 1 .. number of stored/interpolated values
    UInt noInt = 0;  // counts number of interpolated values
    for ( UInt i=1; i < coeff_.size(); i++) {

      if ( solMemoryVal_[i-1] == trueVAL ) {
        newCoeff[i - noInt] += coeff_[i];
      }
      else if (solMemoryVal_[i-1] == LIN1PT ) {
        newCoeff[i - noInt - 1] += 0.5 * coeff_[i];
        newCoeff[i - noInt] += 0.5 * coeff_[i];

        noInt++;
      }
    }
    // output weight factors coeff_
    // 	PrintSolMemoryVal();
    // 	std::cout << "Weight factors of BDF are:" << std::endl;
    // 	for(UInt i=0; i<coeff_.size(); i++)
    // 	  std::cout << coeff_[i] << "  ";
    // 	std::cout << std::endl << "Compressed Weight factors are:" << std::endl;
    // 	for(UInt i=0; i<newCoeff.size(); i++)
    // 	  std::cout << newCoeff[i] << "  ";
    // 	std::cout << std::endl;

    // assign compressed weight factors to coeff_
    coeff_ = newCoeff;
  }

  void NewmarkFracDamp::PrintSolMemoryVal() {

    ENTER_FCN( "NewmarkFracDamp::PrintSolMemoryVal" );

    std::string msg="solMemory_ is:";
    for (UInt i=0; i<solMemoryVal_.size(); i++) {
      if ( solMemoryVal_[i] == NOTUSED )
        msg += " N";
      else if ( solMemoryVal_[i] == trueVAL )
        msg += " T";
      else if ( solMemoryVal_[i] ==  LIN1PT )
        msg += " L";
    }
    msg += "\n\n";
    //Info->PrintF(pdename_, msg.c_str() );
    std::cout << msg;
  }

} // end of namespace
