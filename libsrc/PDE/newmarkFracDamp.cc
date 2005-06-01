#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "Forms/forms_header.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "basePDE.hh"
#include "DataInOut/MaterialData.hh"
#include "Utils/mathfunctions.hh"

namespace CoupledField {

  NewmarkFracDamp::NewmarkFracDamp( BaseSystem * algebraicsystem,
                                    UInt rhsSize, 
                                    const PdeIdType apdeId,
                                    NodeEQN * ptEQN, Grid * aptgrid,
                                    StdPDE * aptStdPDE,
                                    StdVector<RegionIdType> asubdomainList,
                                    StdVector<DampingType> adampingList, 
                                    UInt afracMemory, 
                                    InterpolType ainType, Boolean isaxi)
    :TimeStepping( algebraicsystem, rhsSize )
  {
    ENTER_FCN( "NewmarkFracDamp::NewmarkFracDamp" );
  
    pdename_     = aptStdPDE->GetName();
    pdeId_       = apdeId;
    ptgrid_      = aptgrid;
    ptStdPDE_   =  aptStdPDE;
    ptEQN_ = ptEQN;

    subdoms_     = asubdomainList;
    dampingList_ = adampingList;
    fracMemory_  = afracMemory;
    inType_      = ainType;
    isaxi_       = isaxi;
  
    alpha_ = 0.0;
    beta_  = 0.25;
    gamma_ = 0.5;

    //check if integration parameters are defined
    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    if(analysis != "paramIdent")
      Info->PrintF( pdename_,
                    "Newmark: Using defaults for alpha, beta and gamma!\n" );


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

  NewmarkFracDamp::~NewmarkFracDamp()
  {
    ENTER_FCN( "NewmarkFracDamp::~NewmarkFracDamp" );

  }

  void NewmarkFracDamp::Init( std::map<FEMatrixType,Double> & matrix_factors,
                              Double dt ) 
  {
    ENTER_FCN( "NewmarkFracDamp::Init" );

    dt_ = dt;
    CalcParameters(dt_);

    matrix_factors[STIFFNESS] = 1.0;
    // matrix_factors[DAMPING] = 0.0;
    matrix_factors[DAMPING] = 1.0*a2_;
    matrix_factors[CONVECTION] = 0.0;
    matrix_factors[MASS] = 1.0*a2_;

  }


  void NewmarkFracDamp::Predictor(Vector<Double>& solold)
  {
    ENTER_FCN( "NewmarkFracDamp::Predictor" );

    laststepcalc_ = ptStdPDE_->GetTimeStepCounter();

    // determine number of terms over which BDF is calculated
    //   assumes first nstep = 1 (see transientdriver.cc)!
    calclimit_ = solMemoryVal_.size(); 

    solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
    solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
  }


  void NewmarkFracDamp::UpdateRHS()
  {
    ENTER_FCN( "NewmarkFracDamp::UpdateRHS" );

    Vector<Double> coeffMass, coeffDamp;

    // mass part
    coeffMass = solpred_*a2_;
    algsys_->UpdateRHS(MASS,coeffMass.GetPointer());

    // additional term for fractional damping part
    Matrix<Double> elemmat;         // element matrix
    Matrix<Double> ptCoord;
    BaseFE         * ptElem;
    StdVector<UInt> connecth;
    Vector<Double> rhsAssemble, rhsvec, elemsol; 
  
    UInt noInt;
  
    Double density, compressibility, c0, alpha0, y, factor;

    MaterialData *mymaterialData;
    mymaterialData = ptStdPDE_->getPDEMaterialData();

    if ( dampingList_.GetSize() != subdoms_.GetSize() )
      Error("Mismatch between dampingList_ and subdoms_!", 
            __FILE__, __LINE__);  
    
    for ( UInt actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {

	  if ( dampingList_[actSD]==THERMOVISCOUS || 
		   dampingList_[actSD]==RAYLEIGH ) {
		
		Vector<Double> coeffDamp;
		
		coeffDamp = -solderiv1pred_ + solpred_*a4_;
		algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
		
		Info->PrintF(pdename_,"  THERMOVISCOUS/RAYLEIGH damping for region %d\n"
					 , actSD );
	  }
	  else if ( dampingList_[actSD] == FRACTIONAL_GL || 
				dampingList_[actSD] == FRACTIONAL_BLANK ) {

		density         = mymaterialData[actSD].GetDensity();
		compressibility = mymaterialData[actSD].GetCompressibility();
		c0 = sqrt(compressibility/density);
		alpha0 = mymaterialData[actSD].GetDampingAlfa();
		y      = mymaterialData[actSD].GetDampingBeta();

		// factor: pre factor in Newmark timestepping scheme
		// coeff_: weight factors in BDF
		if ( dampingList_[actSD] == FRACTIONAL_GL ) {
		  factor = density * 2.0 * alpha0 / c0 / sin((y-1.0)*PI/2.0);
		  factor *= exp(-(y-1.0)*log(dt_));
		  factor *= a2_;
		  GLWeights(calclimit_, y);
		}
		else if ( dampingList_[actSD] == FRACTIONAL_BLANK ) {
		  factor =  density * 2.0 * alpha0 / c0 / sin((y-1.0)*PI/2.0);
		  factor *= exp(-(y-1.0)*log(dt_)) * exp(-gammaln(1.0- (y- 1.0)) );
		  factor *= a2_;
		  BlankWeights(calclimit_, y, TRUE);
		}

		// output weight factors coeff_
		//std::cout << "Weight factors of BDF are:" << coeff_ << std::endl;

		// get elements belonging to actual subdomain   
		StdVector<Elem*> elemssd;
		ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
        
		for (UInt el=0; el < elemssd.GetSize(); el++) {
		  // pointer to element
		  ptElem=elemssd[el]->ptElem;
		  BaseForm * bilinear_mass = new MassInt(ptElem, factor, isaxi_);
          
		  connecth=elemssd[el]->connect;
		  ptgrid_->GetElemNodesCoord(ptCoord,connecth);
          
		  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
          
		  // map connect to PDE node numbers
		  StdVector<Integer> connect_PDE;
		  ptEQN_->Node2EQN(connecth, connect_PDE);

#ifdef DEBUG
		  // output matrix with which BDF is computed
		  (*debug) << "Mass Matrix, damping part RHS, of Element " 
				   << el << std::endl;
		  (*debug) << elemmat << std::endl;
#endif
		  // compute BDF
		  GetElemSolution(solpred_, elemsol, connect_PDE);
		  rhsvec = -elemsol * coeff_[0];
		  noInt = 0; // counts how many values have allready been interpolated
		  for (UInt i=1; i<=calclimit_; i++) {

			if ( solMemoryVal_[i-1] == TRUEVAL ) {
			  GetElemSolution(solMemory_[i-1-noInt], elemsol, connect_PDE);
			  rhsvec += elemsol * coeff_[i];
			}
			else if ( solMemoryVal_[i-1] == LIN1PT ) {
			  noInt++;
			  GetElemSolution(solMemory_[i-1-noInt], elemsol, connect_PDE);
			  rhsvec += elemsol * 0.5 * coeff_[i];
			  GetElemSolution(solMemory_[i-1-noInt+1], elemsol, connect_PDE);
			  rhsvec += elemsol * 0.5 * coeff_[i];
			}
		  }
		  rhsAssemble = elemmat * rhsvec;

#ifdef DEBUG
		  (*debug) <<  "BDF vector of timestep " << laststepcalc_ << std::endl;
		  (*debug) << rhsvec << std::endl;
#endif
          
		  //assemble to RHS
		  algsys_->SetElementRHS(&rhsAssemble[0], pdeId_, 
								 connect_PDE.GetPointer(),
								 connect_PDE.GetSize());
          
		  delete bilinear_mass;
		}
		Info->PrintF(pdename_, "  FRACTIONAL damping for region  %d\n",actSD);
		Info->PrintF(pdename_, "    BDF over %d calculated\n", calclimit_ );
		PrintSolMemoryVal();
	  }
	}
  }


  void NewmarkFracDamp::Corrector(Vector<Double>& solnew)
  {
    ENTER_FCN( "NewmarkFracDamp::Corrector" );

    solderiv2_ = (solnew - solpred_) * a2_;
    solderiv1_ = solderiv1pred_ + solderiv2_*a3_;

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

      if (laststepcalc_ <= fracMemory_)
        solMemoryVal_.insert( solMemoryVal_.begin(),TRUEVAL);
      // when solMemoryVal_ reaches size fracMemory_ , all entries are TRUEVAL
      //  and stay the same for all time steps

    }
    // linear interpolation of one solution value 
    else if (inType_ == LIN1PT) {

      if (laststepcalc_ <= fracMemory_) {
        for (UInt i=fracMemory_-1; i>=1; i--)
          solMemory_[i] = solMemory_[i-1];

        solMemoryVal_.insert( solMemoryVal_.begin(),TRUEVAL);
      }
      else if ( (laststepcalc_ > fracMemory_) 
                && (laststepcalc_ < 2*fracMemory_) ) {
        for (UInt i=(2*fracMemory_-laststepcalc_-1); i>=1; i--)
          solMemory_[i] = solMemory_[i-1];
		
        solMemoryVal_[0] = TRUEVAL;
        solMemoryVal_.insert( (solMemoryVal_.begin()
                               + 2*fracMemory_-laststepcalc_),LIN1PT);
	  }
      else if (laststepcalc_ >= 2*fracMemory_ ) {
        if ( (laststepcalc_%2)==0 ) {
          for (UInt i=fracMemory_-1; i>=1; i--)
            solMemory_[i] = solMemory_[i-1];

          solMemoryVal_.pop_back();
          solMemoryVal_.pop_back();
          solMemoryVal_.insert(solMemoryVal_.begin(), TRUEVAL);
        }
        else if ( (laststepcalc_%2)==1 )
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

    //for RHS, Matrices
    a4_ = gamma_ / (beta_*dt_);
  }


  void NewmarkFracDamp::
  GetElemSolution ( const Vector<Double>& sol, 
                    Vector<Double>& elemsol,
                    const StdVector<Integer> & connectPDE )
  {
    ENTER_FCN( "NewmarkFracDamp::GetElemSolution" );
   
    elemsol.Resize(connectPDE.GetSize());
    for (UInt eqn=0; eqn<connectPDE.GetSize(); eqn++) 
      elemsol[eqn] = sol[connectPDE[eqn]-1];
  }

  void NewmarkFracDamp::GLWeights(UInt memory, Double y )
  {
    ENTER_FCN( "NewmarkFracDamp::GLWeights" );
   
    // reserve memory for weights of BDF, order of derivative is y-1
    coeff_.resize(memory+1);
   
    coeff_[0] = 1.0; // Index 0

    for (UInt i=1; i <= memory; i++) // Index 1 .. memory
      coeff_[i] = coeff_[i-1] * ((i-1) - (y-1))/i;
  }


  void NewmarkFracDamp::BlankWeights(UInt memory, Double y, Boolean full)
  {
    ENTER_FCN( "NewmarkFracDamp::BlankWeights" );
 
    Double pot;
    pot = 1.0- (y - 1.0);

    // reserve memory for weights of BDF, order of derivative is y-1
    coeff_.resize(memory+1);
      
    coeff_[0] = 1.0/ pot; // Index 0
 
    if (full) {
      for (UInt i=1; i<memory; i++) // Index 1 .. (memory-1)
        coeff_[i]= ( exp(pot* log(i-1))
                     - 2* exp(pot* log(i))
                     + exp(pot* log(i+1)) ) / pot;
      // Index memory
      coeff_[memory]= exp(-(y-1.0)*log(fracMemory_)) + 
        ( exp(pot* log(fracMemory_ -1))
          - exp(pot* log(fracMemory_)) ) / pot;
    }
    else {
      for (UInt i=1; i<=memory; i++) // Index 1 .. memory
        coeff_[i]= ( exp(pot* log(i-1))
                     - 2* exp(pot* log(i))
                     + exp(pot* log(i+1)) ) / pot;
    }
  }

  void NewmarkFracDamp::PrintSolMemoryVal()
  {
    ENTER_FCN( "NewmarkFracDamp::PrintSolMemoryVal" );

    std::string msg="solMemory_ is:";
    for (UInt i=0; i<solMemoryVal_.size(); i++) {
      if ( solMemoryVal_[i] == NOTUSED )
        msg += " N";
      else if ( solMemoryVal_[i] == TRUEVAL )
        msg += " T";
      else if ( solMemoryVal_[i] ==  LIN1PT )
        msg += " L";
    }
    msg += "\n\n";
    Info->PrintF(pdename_, msg.c_str() );
  }

} // end of namespace
