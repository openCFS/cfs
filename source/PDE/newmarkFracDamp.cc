#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "newmarkFracDamp.hh"
#include "Forms/massInt.hh"
#include "DataInOut/WriteInfo.hh"
#include "basePDE.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/mathfunctions.hh"
#include "Domain/domain.hh"
#include "Driver/singleDriver.hh"
#include "General/exception.hh"

namespace CoupledField {

  NewmarkFracDamp::NewmarkFracDamp( BaseSystem * algebraicsystem,
                                    const PdeIdType apdeId,
                                    shared_ptr<EqnMap> eqnMap,
                                    Grid * aptgrid,
                                    StdPDE * aptStdPDE,
                                    shared_ptr<ResultInfo> aresult,
                                    StdVector<RegionIdType> asubdomainList,
                                    std::map<RegionIdType,DampingType> adampingList) 
    :TimeStepping( algebraicsystem ){
	
    ENTER_FCN( "NewmarkFracDamp::NewmarkFracDamp" );
  
    pdename_     = aptStdPDE->GetName();
    pdeId_       = apdeId;
    ptgrid_      = aptgrid;
    ptStdPDE_    = aptStdPDE;
    eqnMap_      = eqnMap;
    result_      = aresult;
	
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
    // std::cerr << "In NewmarkFracDamp isaxi_ has value of: " << isaxi_ << std::endl;

  
    alpha_ = 0.0;
    beta_  = 0.25;
    gamma_ = 0.5;

    //check if integration parameters are defined
    std::string analysis;
    Info->PrintF( pdename_, "NewmarkFracDamp: Using defaults for alpha, \
beta and gamma!\n" );

 
  }

  NewmarkFracDamp::~NewmarkFracDamp() {

    ENTER_FCN( "NewmarkFracDamp::~NewmarkFracDamp" );

    delete [] solMemory_;
  }

  void NewmarkFracDamp::Init( Double dt, UInt rhsSize ) {
	
    ENTER_FCN( "NewmarkFracDamp::Init" );

    rhsSize_ = rhsSize;
    dt_ = dt;
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = 1.0;
    matrix_factors_[DAMPING] = 1.0*a4_; // needed for thermoviscous damping
    matrix_factors_[CONVECTION] = 0.0;
    matrix_factors_[MASS] = 1.0*a2_;


    // get the memory
    if( !isDeriv1Set_ ) {
      solderiv1_.Resize(rhsSize_);  
      solderiv1_.Init();
    }

    if( !isDeriv2Set_ ) {
      solderiv2_.Resize(rhsSize_);  
      solderiv2_.Init();
    }
  

    solpred_.Resize(rhsSize_); 
    solpred_.Init();
    solderiv1pred_.Resize(rhsSize_); 
    solderiv1pred_.Init();
  
    if ( fracMemory_ <= 1 ) {
      EXCEPTION("Attenuation model needs frac_memory larger than 1" );
    } else {
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

    actStep_ = domain->GetSingleDriver()->GetActStep( pdename_ );

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
    StdVector<Integer> eqnNumbers;
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
        BaseForm * bilinear_mass = new MassInt(factor, 1, isaxi_);
        bilinear_mass->SetAnsatzFct( result_->fctType, result_->fctType );

        ElemList actSDList(ptgrid_ );
        actSDList.SetRegion( subdoms_[actSD] );
        EntityIterator it = actSDList.GetIterator();

        // iterate over all elements of subdomain
        for ( it.Begin(); !it.IsEnd(); it++ ) {

          bilinear_mass->CalcElementMatrix(elemmat,it,it);
          
          // map node numbers to equation numbers
          eqnMap_->GetEqns( eqnNumbers, *result_, it );

#ifdef DEBUG
          UInt elemNum = it.GetElem()->elemNum;
		  
          // output matrix with which BDF is computed
          (*debug) << "Damping   matrix (fractional Damping) of Element " 
                   << elemNum << std::endl;
          (*debug) << elemmat << std::endl;
#endif

          // compute BDF
          GetElemSolution( solpred_, elemsol, eqnNumbers );
          rhsvec = -elemsol * coeff_[0];

          for (UInt i=1; i<numTrueValues_; i++) {

            GetElemSolution(solMemory_[i-1], elemsol, eqnNumbers );
            rhsvec += elemsol * coeff_[i];
          }

          rhsAssemble = elemmat * rhsvec;

#ifdef DEBUG
          (*debug) << "BDF vector of timestep " << actStep_ << std::endl
                   << rhsvec << std::endl;
#endif
          
          //assemble to RHS
          algsys_->SetElementRHS(&rhsAssemble[0], pdeId_, 
                                 eqnNumbers.GetPointer(),
                                 eqnNumbers.GetSize());
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
                    const StdVector<Integer> eqnNumbers  ) {

    ENTER_FCN( "NewmarkFracDamp::GetElemSolution" );

    elemsol.Resize(eqnNumbers.GetSize());
    for (UInt eqn=0; eqn<eqnNumbers.GetSize(); eqn++) {
      if( eqnNumbers[eqn] != 0 ) {
        elemsol[eqn] = sol[eqnNumbers[eqn]-1];
      } else
        elemsol[eqn] = 0.0;
    }
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

    // output weight factors vector 'coeff_'
#ifdef DEBUG
    std::string mymessage;
    PrintSolMemoryVal(mymessage);
    (*debug) << mymessage << std::endl << std::endl;

    (*debug) << "Weight factors of BDF are:" << std::endl;
    for(UInt i=0; i<coeff_.size(); i++)
      (*debug) << coeff_[i] << "  ";
    (*debug) << std::endl << "Compressed Weight factors are:" << std::endl;
    for(UInt i=0; i<newCoeff.size(); i++)
      (*debug) << newCoeff[i] << "  ";
    (*debug) << std::endl;
#endif

    // assign compressed weight factors to coeff_
    coeff_ = newCoeff;
  }

  void NewmarkFracDamp::PrintSolMemoryVal(std::string & msg) {

    ENTER_FCN( "NewmarkFracDamp::PrintSolMemoryVal" );

    msg="BDF will be build from (T=true value, L=linear interpolated):\n";
    for (UInt i=0; i<solMemoryVal_.size(); i++) {
      if ( solMemoryVal_[i] == NOTUSED )
        msg += " N";
      else if ( solMemoryVal_[i] == trueVAL )
        msg += " T";
      else if ( solMemoryVal_[i] ==  LIN1PT )
        msg += " L";
    }

    //Info->PrintF(pdename_, msg.c_str() );
    //std::cerr << msg;
  }

} // end of namespace
