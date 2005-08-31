#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "Forms/forms_header.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "basePDE.hh"
#include "DataInOut/MaterialData.hh"
#include "Utils/mathfunctions.hh"

namespace CoupledField {

  NewmarkFracDampMech::
  NewmarkFracDampMech( BaseSystem * algebraicsystem,
                       UInt rhsSize, 
                       const PdeIdType apdeId,
                       NodeEQN * ptEQN,
                       Grid * aptgrid,
                       StdPDE * aptStdPDE,
                       StdVector<RegionIdType> asubdomainList,
                       std::map<RegionIdType,DampingType> adampingList) 
    :TimeStepping( algebraicsystem, rhsSize ){
    
    ENTER_FCN( "NewmarkFracDampMech::NewmarkFracDampMech" );
    
    pdename_     = aptStdPDE->GetName();
    pdeId_       = apdeId;
    ptgrid_      = aptgrid;
    ptStdPDE_    = aptStdPDE;
    ptEQN_       = ptEQN;

    params->Get( "subType", subType_, pdename_ );
	
    subdoms_     = asubdomainList;
    dampingList_ = adampingList;
    fracMemory_  = aptStdPDE->GetFracMemory();
    inType_ = NOTUSED;
    geomType_    = "3d"; //ptStdPDE_->GetSubType();
    isaxi_       = ptStdPDE_->GetIsaxi();

    UInt numElems = ptgrid_->GetNumElems(subdoms_);
    UInt dimStressVector = getStressDim();
  
    alpha_ = 0.0;
    beta_  = 0.25;
    gamma_ = 0.5;

    //check if integration parameters are defined
    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    if(analysis != "paramIdent")
      Info->PrintF( pdename_, "NewmarkFracDampMech: Using defaults for alpha, \
                               beta and gamma!\n" );
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
      Error("Damping model needs frac_memory larger than 1",
            __FILE__,__LINE__);
    else {
      // get the memory
      solMemory_  = new Vector<Double>[fracMemory_];
      stressHistoryEl_ = new Vector<Double>[fracMemory_];

      for (UInt i=0; i<fracMemory_; i++) {
        solMemory_[i].Resize(rhsSize_);  
        solMemory_[i].Init();
        stressHistoryEl_[i].Resize(numElems*dimStressVector);
        stressHistoryEl_[i].Init();
      }
    }
  }


  NewmarkFracDampMech::~NewmarkFracDampMech() {

    ENTER_FCN( "NewmarkFracDampMech::~NewmarkFracDampMech" );
  }

  void NewmarkFracDampMech::Init( std::map<FEMatrixType,Double> & matrix_factors,
                              Double dt ) {
	
    ENTER_FCN( "NewmarkFracDampMech::Init" );

    dt_ = dt;
    CalcParameters(dt_);

    //elastModule_ = 1.0;  
    //elastModule_ = 658.2;  
    StdVector<Double> elastModuleList_;
    params->GetList( "ElastModule", elastModuleList_, "mechanic", "damping" );
    elastModule_ = elastModuleList_[0];

    StdVector<Double> timeSlotList_;
    params->GetList( "timeSlot", timeSlotList_, "mechanic", "damping" );
    Double timeSlot = timeSlotList_[0];
  
    Double temp;
    temp   = (timeSlot/GetTimeStep())/fracMemory_;
  
    if (temp < 1.0)
      std::cerr << "The fracMemory value is to big for the given timeSlot and the timeStepSize" << std::endl;
  
    modulo_ = Integer(temp);
    
    matrix_factors[STIFFNESS] = 1.0;
    matrix_factors[DAMPING] = 1.0*a4_; // needed for thermoviscous damping
    matrix_factors[CONVECTION] = 0.0;
    matrix_factors[MASS] = 1.0*a2_;
  }


  void NewmarkFracDampMech::Predictor(Vector<Double>& solold) {

    ENTER_FCN( "NewmarkFracDampMech::Predictor" );

    actStep_ = ptStdPDE_->GetTimeStepCounter();

    // determine number of terms over which BDF is calculated
    //   assumes first nstep = 1 (see transientdriver.cc)!
    numValues_ = solMemoryVal_.size();

    // determine number of past values stored in solMemory_
    numTrueValues_ = 0;
    for ( UInt i=0; i < numValues_; i++ ) {
      if ( solMemoryVal_[i] == TRUEVAL )
	numTrueValues_++;
    }

    solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
    solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
  }


  void NewmarkFracDampMech::UpdateRHS() {

    ENTER_FCN( "NewmarkFracDampMech::UpdateRHS" );
    
    // mass part
    Vector<Double> coeffMass;

    coeffMass = solpred_*a2_;
    algsys_->UpdateRHS(MASS,coeffMass.GetPointer());

    // damping part
    Matrix<Double>  elemmat,ptCoord;
    BaseFE          * ptElem;
    StdVector<UInt> connecth;
    Vector<Double>  rhsAssemble, rhsvec, elemsol;
    Vector<Double>  fracDerivStressVec_, resultStressVector,stressVector;

    fracDerivStressVec_.Resize(getDim());
    stressVector.Resize(getDim());

    MaterialData    *mymaterialData;
    mymaterialData = ptStdPDE_->getPDEMaterialData();

    dampAlpha_ = mymaterialData[0].GetDampingAlfa();
    dampBeta_ = mymaterialData[0].GetDampingBeta();   

    StdVector<Double> fracDerivList_;
    params->GetList( "fracDeriv", fracDerivList_, "mechanic", "damping" );
    Double fracDeriv_ = fracDerivList_[0];

    Double timeStep = GetTimeStep();
    timeStepPowerFracDeriv_ = std::pow(timeStep,-fracDeriv_);


    if ( dampingList_.size() != subdoms_.GetSize() )
      Error("Mismatch between dampingList_ and subdoms_!", __FILE__, __LINE__);

    std::string model;
    StdVector<std::string> modelList_;
    params->GetList( "model", modelList_, "mechanic", "damping" );
    model =  modelList_[0]; 

    for ( UInt actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {
      if ( dampingList_[subdoms_[actSD]] == NONE ) {
        // no damping term has to be computed
      }
      else if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
	
        Vector<Double> coeffDamp;
        coeffDamp = -solderiv1pred_ + solpred_*a4_;
        algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
      }
      else {
        if ( dampingList_[subdoms_[actSD]]== FRACTIONAL_GL ) {
          //factor *= std::exp(-(y-1.0)*std::log(dt_));
          //GLWeights(numValues_, y);
          GLWeights(numValues_ + 1, fracDeriv_); // calclimit_+1 because the coeffcients in the loop start with i+1
        }
            
        // get elements belonging to actual subdomain   
        StdVector<Elem*> elemssd;
        ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
        
        BDBInt * rhsViscoMat = new LinViscoElastInt(mymaterialData[actSD],
                                                    subType_,
                                                    "MatDepRHSMatrix",
                                                    GetTimeStep() );

        BDInt * bdInt = new BDInt(mymaterialData[actSD],
                                  subType_, 
                                  GetTimeStep());

        for (UInt actEl=0; actEl < elemssd.GetSize(); actEl++) {
          ptElem=elemssd[actEl]->ptElem;
          const UInt nrNodes  = ptElem->GetNumNodes();
          dofs_ = ptElem->GetDim();
          
          resultStressVector.Resize(nrNodes * dofs_);
          rhsvec.Resize(nrNodes * dofs_);
          rhsvec.Init();
              
          rhsViscoMat->SetElemPtr(ptElem);
          bdInt->SetElemPtr(ptElem);
              
          connecth=elemssd[actEl]->connect;
          ptgrid_->GetElemNodesCoord(ptCoord,connecth);
          
          rhsViscoMat->CalcElementMatrix(ptCoord, elemmat);
          
          // map connect to PDE node numbers
          StdVector<Integer> connect_PDE;
          ptEQN_->Node2EQN(connecth, connect_PDE);

#ifdef DEBUG
          //UInt elemNum = elemssd[actEl]->elemNum;
          // output matrix with which BDF is computed
          (*debug) << "fractional Damping matrix of Element" << std::endl;
          (*debug) << elemmat << std::endl;
          (*debug) << "actStep_=" << actStep_ << std::endl;
          (*debug) << "numValues_=" << numValues_ << std::endl;
#endif
              
          for (UInt i=1; i<=numValues_; i++) {
            if ( solMemoryVal_[i-1] == TRUEVAL )
              {
                GetElemSolution(solMemory_[i-1], elemsol, connect_PDE);
                rhsvec += elemsol * coeff_[i];	     
#ifdef DEBUG
                (*debug) << "elemsolOld" << std::endl;
                (*debug) << elemsol << std::endl;
#endif
              }
            else { // see NewmarkFracDamp
            }		
          }
              
          rhsAssemble = -elemmat  * (rhsvec * timeStepPowerFracDeriv_);
          
          // compute stress history from displacement history with the fractional constitutive equation
          // later perhabs  this method should be implemented in a own class
          
          // calculates the stress history from the displacement history	  
          //CalcStressHistoryForElement( connect_PDE, ptElem->GetNumNodes(),ptCoord,mymaterialData[actSD],ptElem,el);
          // calculate linear stresses for stress history
          //	  CalcStressHist( ptCoord, mymaterialData[actSD],connect_PDE,ptElem);
          
          fracDerivStressVec_.Init();
          resultStressVector.Init();
          
          for (UInt i=1; i<numValues_; i++) {	     
            GetStressVector(stressVector,actEl,i-1);
            fracDerivStressVec_ += stressVector * coeff_[i+1];
          }
              
          bdInt->calcElementVector(ptCoord,resultStressVector,fracDerivStressVec_);
          
          //                   double t1,t2;
          //                   t1=0;
          //                   t2=0;
              
          //                   for (UInt i=0;i<rhsAssemble.GetSize();i++)
          //                     {
          //                       t1 += rhsAssemble[i];
          //                       t2 += resultStressVector[i];
          //                     }
              
#ifdef DEBUG
          (*debug) << "actStep_=" << actStep_ << std::endl 
                   << "rhsAssemble=" << std::endl << rhsAssemble << std::endl
                   << "resultStressVector=" << std::endl << resultStressVector << std::endl;
#endif
          if(model== "3param")
            rhsAssemble  = rhsAssemble +  resultStressVector;
          else if (model=="KelvinVoigt")
            rhsAssemble  = rhsAssemble; // +  resultStressVector;
          else
            std::cerr << "unknown model for fractional damping" << std::endl;
#ifdef DEBUG
          //(*debug) <<  "rhs vector of timestep " << actStep_ << std::endl;
          //(*debug) << rhsvec << std::endl;
#endif
          //assemble to RHS
          algsys_->SetElementRHS(&rhsAssemble[0], pdeId_, connect_PDE.GetPointer(),
                                 connect_PDE.GetSize());
              
        }
      }
    }
  }
  

  void NewmarkFracDampMech::Corrector(Vector<Double>& solnew)
  {
    ENTER_FCN( "NewmarkFracDampMech::Corrector" );

    solderiv2_ = (solnew - solpred_) * a2_;
    solderiv1_ = solderiv1pred_ + solderiv2_*a3_;

    Integer numEQNs = ptEQN_->GetNumEQNs();
    StdVector<UInt> connecth, connect_PDE;
    Vector<Double> stressVec, displacementVec;
    Matrix<Double> ptCoord;

    stressVec.Resize(getDim());
    displacementVec.Resize(numEQNs * dofs_);  

    MaterialData *mymaterialData;
    BaseFE         * ptElem;
    mymaterialData = ptStdPDE_->getPDEMaterialData();

    if( (actStep_ % modulo_) == 0) {
      for ( Integer i=fracMemory_-1; i>=1; i-- ) {
        solMemory_[i] = solMemory_[i-1];
        stressHistoryEl_[i] = stressHistoryEl_[i-1];
      }
      solMemory_[0] = solnew; //- solpred_;
      
      if ((actStep_ / modulo_)  <=  fracMemory_) {
        solMemoryVal_.insert( solMemoryVal_.begin(),TRUEVAL);
        // when solMemoryVal_ reaches size fracMemory_ , all entries are TRUEVAL
        //  and stay the same for all time for
      }
      
      for ( UInt actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {
        StdVector<Elem*> elemssd;
        ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
        //ptgrid_->GetElemSD(elemssd,subdoms_[actSD]);
          
        for (UInt el=0; el < elemssd.GetSize(); el++) {
          ptElem=elemssd[el]->ptElem;
          connecth=elemssd[el]->connect;
          ptgrid_->GetElemNodesCoord(ptCoord,connecth);
          StdVector<Integer> connect_PDE;
          ptEQN_->Node2EQN(connecth, connect_PDE);
          
          GetElemSolution(solMemory_[0], displacementVec, connect_PDE);
          
          CalcStress( ptElem, mymaterialData[actSD], connect_PDE, ptCoord,  stressVec, displacementVec, el);
          
          InsertStressVector(stressVec,el,0);		
        }
      }
    }
  }

  void NewmarkFracDampMech::CalcParameters(Double dt)
  {
    ENTER_FCN( "NewmarkFracDampMech::CalcParameters" );

    //for predictors
    a0_ = 0.5*(1-2.0*beta_)*dt_*dt_;
    a1_ = (1-gamma_)*dt_;

    //for correctors, matrices
    a2_ = 1.0/(beta_*dt_*dt_);
    a3_ = gamma_*dt_;

    //for RHS, matrices
    a4_ = gamma_ / (beta_*dt_);
  }


  void NewmarkFracDampMech::GetElemSolution ( const Vector<Double>& sol, 
                    Vector<Double>& elemsol,
                    const StdVector<Integer> & connectPDE ) {

    ENTER_FCN( "NewmarkFracDampMech::GetElemSolution" );
   
    elemsol.Resize(connectPDE.GetSize());

    for (UInt eqn=0; eqn<connectPDE.GetSize(); eqn++) {
      if (connectPDE[eqn]>=1)
        elemsol[eqn] = sol[connectPDE[eqn]-1];
    }
  }

  void NewmarkFracDampMech::GLWeights(UInt memory, Double y ) {

    ENTER_FCN( "NewmarkFracDampMech::GLWeights" );
   
    // reserve memory for weights of BDF, order of derivative is y-1
    coeff_.resize(memory+1);
    coeff_[0] = 1.0; // Index 0

#ifdef DEBUG
    (*debug) << "coeff_" <<std::endl;
    (*debug) << coeff_[0]<<std::endl;
#endif

    for (UInt i=1; i <= memory; i++) { // Index 1 .. memory 
      coeff_[i] = coeff_[i-1] * (i-1-y)/(i);

#ifdef DEBUG
      (*debug) << coeff_[i] << "   " << coeff_[i-1] * (i-1-y)/(i) <<std::endl;
#endif

    }
  }

void NewmarkFracDampMech::CalcStress(BaseFE * aptelem, MaterialData & matDa, StdVector<Integer> connect_PDE, 
				      Matrix<Double> & ptCoord, Vector<Double> & stressVector,
				      Vector<Double> &displacementVector, Integer elemNr)
{
  ENTER_FCN( "NewmarkFracDampMech::CalcStress" );
  
  Matrix<Double> bMat;
  Matrix<Double> dMat;
  Matrix<Double> alphaMat;
  Matrix<Double> betaMat;
  Matrix<Double> aMat;


  double timeStep = GetTimeStep();

  Integer dimStressVec = getDim();
  //Integer numEQNs = ptEQN_->GetNumEQNs();

  const Integer nrNodes  = aptelem->GetNumNodes();

  BDBInt * viscoIntegrator = new LinViscoElastInt(aptelem,matDa,subType_,
    							"modifiedStiffness",timeStep );
  Vector<Double> fracDerivStress, fracDerivDisplacement;
  Vector<Double> term1, term2,term3;

  fracDerivStress.Resize(dimStressVec);
  stressVector.Resize(dimStressVec);
  fracDerivDisplacement.Resize(nrNodes * dofs_);
  term1.Resize(dimStressVec);
  term2.Resize(dimStressVec);
  term3.Resize(dimStressVec);

  fracDerivStress.Init();
  fracDerivDisplacement.Init();

  viscoIntegrator->GetDMat(dMat);
  viscoIntegrator->GetBMat(bMat,ptCoord);
  GetBetaMat(betaMat,elastModule_, matDa);
  GetAlphaMat(alphaMat);
  GetAMat(aMat);

  //computation of the term1 (actual displacement * material, damping factors)
  term1 = dMat * bMat * displacementVector;

  //computation of term2 (fractional Derivative of displacement)
  // solMemory[0] is the displacement value of the previous time step -> loop from 0
  for(UInt k=0; k< numValues_;k++) {
    GetElemSolution(solMemory_[k], displacementVector, connect_PDE);	
    fracDerivDisplacement += displacementVector *  coeff_[k+1];   
  }

  term2 = aMat * betaMat * bMat * fracDerivDisplacement;

  //computation of term3 (fractional derivative of stress)
  for(UInt k=0; k< numValues_;k++) {
    GetStressVector(stressVector, elemNr,k);	
    fracDerivStress += stressVector *  coeff_[k+1];   
  }

  term3 =    aMat * alphaMat * fracDerivStress;

  double t1,t2,t3;
    t1=0;
    t2=0;
    t3=0;

    for (UInt i=0;i<term1.GetSize();i++) {
      t1 += term1[i];
      t2 += term2[i];
      t3 += term3[i];
    }
    
    stressVector = term1 + term2 - term3;  
}

  void NewmarkFracDampMech::GetStressVector(Vector<Double> &stressVec,UInt elemNr,UInt memory)
  {

    UInt dimStressVec = stressVec.GetSize();
    for(UInt i = 0; i < dimStressVec;i++) {
      stressVec[i]=stressHistoryEl_[memory][(elemNr*dimStressVec)+i];
    }
  }

void NewmarkFracDampMech::InsertStressVector(Vector<Double> &stressVec, Integer elemNr, Integer memory)
{
  Integer dimStressVec = stressVec.GetSize();

  for(Integer i = 0; i < dimStressVec;i++) {
    stressHistoryEl_[memory][(elemNr*dimStressVec)+i] = stressVec[i];
  }
}



void NewmarkFracDampMech::GetAMat(Matrix<Double> & aMat)
   { 
    ENTER_FCN( "NewmarkFracDampMech::GetAMat" );

     double val = 0.0;
     
     // compute matrix A,  same entries, 
     aMat.Resize(getDim());
     // set entries on the diagonal
      val = (timeStepPowerFracDeriv_)  * dampAlpha_ ;
	
     // set the emtries on the diagonal matrix, to get the inverse, 
     // the value on the diagonal are 1/val

     for(UInt i=0;i<getDimD();i++) {
       aMat.SetEntry(i,i,1/(1+val));
       //aMat.SetEntry(i,i,1/(val));
     }  
   }


void NewmarkFracDampMech::GetAlphaMat(Matrix<Double> & alphaMat)
   { 
    ENTER_FCN( "NewmarkFracDampMech::GetAlphaMat" );

     double val = 0.0;
     
     // compute matrix A,  same entries, 
     alphaMat.Resize(getDim());
     // set entries on the diagonal
      val = (timeStepPowerFracDeriv_)  * dampAlpha_ ;
	
     // set the emtries on the diagonal matrix, to get the inverse, 
     // the value on the diagonal are 1/val
     for(UInt i=0;i<getDimD();i++) {
       alphaMat.SetEntry(i,i,val);
     }  
   }


void NewmarkFracDampMech::GetBetaMat(Matrix<Double>& betaMat,  Double E, MaterialData & matData)
   { 
    ENTER_FCN( "NewmarkFracDampMech::GetBetaMat" );

     double val = 0.0;
     
     if(subType_ == "axi") {	 
       BDBInt * matMat = new mechAxiInt(matData);
       matMat->GetDMat(betaMat);
     }
       
     else if(subType_ == "planeStrain") {
       BDBInt * matMat = new mechPlainStrainInt(matData);
       matMat->GetDMat(betaMat);
     }
     else if(subType_ == "3d") {
       BDBInt * matMat = new mech3DInt(matData);
       matMat->GetDMat(betaMat);
     }     
     
     val = (dampBeta_/(E)) * timeStepPowerFracDeriv_;
     betaMat  = (betaMat *  val);
   }

  
  UInt NewmarkFracDampMech::getStressDim() {

    ENTER_FCN( "NewmarkFracDampMech::getStressDim" );

    if(subType_ == "axi") {
      return 4;
    }
    else if(subType_ == "planeStrain") {
      return 3;   	
    }
    else if(subType_ == "3d") {
      return 6;   	
    }
    else {
      Error("wrong subType(axi,planeStrain,3d) specified", __FILE__, __LINE__);  
      return 0;
    }
  }
  
  UInt NewmarkFracDampMech::getDim() {
    ENTER_FCN( "NewmarkFracDampMech::getDim" );
    
    if(subType_ == "axi") {
      return 4;
    }
    else if(subType_ == "planeStrain") {
      return 3;   	
    }
    else if(subType_ == "3d") {
      return 6;   	
    }
    else {
      Error("wrong subType(axi,planeStrein,3d) specified", __FILE__, __LINE__);  
      return 0;
    }
  }
} // end of namespace
