// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "Forms/bdInt.hh"
#include "Forms/bdbInt.hh"
#include "Forms/linViscoElastInt.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "basePDE.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/mathfunctions.hh"
#include "Domain/domain.hh"
#include "Driver/singleDriver.hh"
#include "General/exception.hh"
#include "OLAS/algsys/algebraicSys.hh"

namespace CoupledField {

  NewmarkFracDampMech::
  NewmarkFracDampMech( AlgebraicSys * algebraicsystem,
                       const FeFctIdType apdeId,
                       shared_ptr<EqnMap> eqnMap,
                       Grid * aptgrid,
                       StdPDE * aptStdPDE,
                       StdVector<RegionIdType> asubdomainList,
                       std::map<RegionIdType,DampingType> adampingList,
                       PtrParamNode systemNode ) 
    :TimeStepping( algebraicsystem){
    
    
    EXCEPTION( "This version is not working properly anymore and should "
               << "be reimplemented by Gerhard" );
 //    pdename_     = aptStdPDE->GetName();
//     pdeId_       = apdeId;
//     ptgrid_      = aptgrid;
//     ptStdPDE_    = aptStdPDE;
//     eqnMap_      = eqnMap;

//     // get subType of pde
//     UInt actSequenceStep = domain->GetSingleDriver()->GetActSequenceStep();
//     PtrParamNode actStepNode = 
//       param->Get("sequenceStep", "index", GenStr(actSequenceStep));
//     mechNode_ = actStepNode->Get("pdeList")->Get("mechanic");
//     subType_ = mechNode_->As<std::string>();
	
//     subdoms_     = asubdomainList;
//     dampingList_ = adampingList;
//     fracMemory_  = aptStdPDE->GetFracMemory();
//     inType_ = NOTUSED;
//     geomType_    = "3d"; //ptStdPDE_->GetSubType();
//     isaxi_       = ptStdPDE_->GetIsaxi();
  
//     alpha_ = 0.0;
//     beta_  = 0.25;
//     gamma_ = 0.5;

//     //check if integration parameters are defined
//     std::string analysis;
//     if( actStepNode->Get("analysis")->GetChild()->As<std::string>()  != "paramIdent" ) {
//       Info->PrintF( pdename_, "NewmarkFracDampMech: Using defaults for alpha, \
//                                beta and gamma!\n" );
//     }
  }


  NewmarkFracDampMech::~NewmarkFracDampMech() {

  }

  void NewmarkFracDampMech::Init( Double dt, UInt rhsSize ) {
	
    REFACTOR;
//    dt_ = dt;
//    rhsSize_ = rhsSize;
//    CalcParameters(dt_);
//    Vector<Double> dummyVec;
//    dummyVec.Resize(rhsSize_);
//    dummyVec.Init();
//
//    //elastModule_ = 1.0;  
//    //elastModule_ = 658.2;  
//    PtrParamNode firstRegionNode = (mechNode_->Get("regionList")->GetList("region"))[0];
//    firstRegionNode->Get("damping")->GetValue("ElastModul", elastModule_ );
//    
//    Double timeSlot;
//    firstRegionNode->Get("damping")->GetValue("timeSlot", timeSlot );
//  
//    Double temp;
//    temp   = (timeSlot/GetTimeStep())/fracMemory_;
//  
//    if (temp < 1.0)
//      std::cerr << "The fracMemory value is to big for the given timeSlot and the timeStepSize" << std::endl;
//  
//    modulo_ = Integer(temp);
//    
//    matrix_factors_[STIFFNESS] = 1.0;
//    matrix_factors_[DAMPING] = 1.0*a4_; // needed for thermoviscous damping
//    matrix_factors_[CONVECTION] = 0.0;
//    matrix_factors_[MASS] = 1.0*a2_;
//
//
//    // get the memory
//    if ( !is_Deriv_set(FIRST_DERIV) )
//    {
//      solDeriv_vec_[FIRST_DERIV] = dummyVec;
//    }
//    if ( !is_Deriv_set(SECOND_DERIV) )
//    {
//      solDeriv_vec_[SECOND_DERIV] = dummyVec;
//    }
//  
//
//    solpred_.Resize(rhsSize_); 
//    solpred_.Init();
//    solderiv1pred_.Resize(rhsSize_); 
//    solderiv1pred_.Init();
//  
//    if ( fracMemory_ <= 1 ) {
//      EXCEPTION("Damping model needs frac_memory larger than 1" );
//    } else {
//
//      //UInt numElems = ptgrid_->GetNumElems(subdoms_);
//      UInt numElems = eqnMap_->GetNumLocalElems();
//      UInt dimStressVector = getStressDim();
//
//      // get the memory
//      solMemory_  = new Vector<Double>[fracMemory_];
//      stressHistoryEl_ = new Vector<Double>[fracMemory_];
//
//      for (UInt i=0; i<fracMemory_; i++) {
//        solMemory_[i].Resize(rhsSize_);  
//        solMemory_[i].Init();
//        stressHistoryEl_[i].Resize(numElems*dimStressVector);
//        stressHistoryEl_[i].Init();
//      }
//    }
  }


  void NewmarkFracDampMech::Predictor(Vector<Double>& solold) {

    REFACTOR;
//    actStep_ = domain->GetSingleDriver()->GetActStep( pdename_ );
//
//    // determine number of terms over which BDF is calculated
//    //   assumes first nstep = 1 (see transientdriver.cc)!
//    numValues_ = solMemoryVal_.size();
//
//    // determine number of past values stored in solMemory_
//    numTrueValues_ = 0;
//    for ( UInt i=0; i < numValues_; i++ ) {
//      if ( solMemoryVal_[i] == trueVAL )
//        numTrueValues_++;
//    }
//    if( !omitFirstPredictor_) {
//      solpred_ = solold + solDeriv_vec_[FIRST_DERIV]*dt_ + solDeriv_vec_[SECOND_DERIV]*a0_;
//      solderiv1pred_ = solDeriv_vec_[FIRST_DERIV] + solDeriv_vec_[SECOND_DERIV]*a1_;
//    } else {
//      omitFirstPredictor_ = false;
//    }
  }


  void NewmarkFracDampMech::UpdateRHS() {

    REFACTOR;
//    // mass part
//    Vector<Double> coeffMass;
//
//    coeffMass = solpred_*a2_;
//    algsys_->UpdateRHS(MASS,coeffMass);
//
//    // damping part
//    Matrix<Double>  elemmat;
////    BaseFE          * ptElem;
//    StdVector<UInt> connecth;
//    Vector<Double>  rhsAssemble, rhsvec, elemsol;
//    Vector<Double>  fracDerivStressVec_, resultStressVector,stressVector;
//
//    fracDerivStressVec_.Resize(getDim());
//    stressVector.Resize(getDim());
//
//    std::map<RegionIdType, BaseMaterial*> mymaterialData;
//    mymaterialData = ptStdPDE_->getPDEMaterialData();
//
//    BaseMaterial * firstMat = mymaterialData.begin()->second;
//
//    firstMat->GetScalar(dampAlpha_,ACOU_ALPHA,Global::REAL);
//    firstMat->GetScalar(dampBeta_,FRACTIONAL_EXPONENT,Global::REAL);
//
//    PtrParamNode firstRegionNode = (mechNode_->GetList("region"))[0];
//    Double fracDeriv_;
//    firstRegionNode->Get("damping")->GetValue("fracDeriv", fracDeriv_ );
//
//    Double timeStep = GetTimeStep();
//    timeStepPowerFracDeriv_ = std::pow(timeStep,-fracDeriv_);
//
//
//    std::string model;
//    firstRegionNode->Get("damping")->GetValue( "model", model );
//
//    for ( UInt actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {
//      if ( dampingList_[subdoms_[actSD]] == NONE ) {
//        // no damping term has to be computed
//      }
//      else if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
//	
//        Vector<Double> coeffDamp;
//        coeffDamp = -solderiv1pred_ + solpred_*a4_;
//        algsys_->UpdateRHS(DAMPING,coeffDamp);
//      }
//      else {
//        if ( dampingList_[subdoms_[actSD]]== FRACTIONAL_GL ) {
//          //factor *= std::exp(-(y-1.0)*std::log(dt_));
//          //GLWeights(numValues_, y);
//          GLWeights(numValues_ + 1, fracDeriv_); // calclimit_+1 because the coeffcients in the loop start with i+1
//        }
//            
//         
//	//transform the type
//	SubTensorType type;
//	String2Enum(subType_,type);
//
//        BDBInt * rhsViscoMat = new LinViscoElastInt(mymaterialData[subdoms_[actSD]],
//                                                    type,
//                                                    "MatDepRHSMatrix",
//                                                    GetTimeStep() );
//
//        BDInt * bdInt = new BDInt(mymaterialData[subdoms_[actSD]],
//                                  subType_, 
//                                  GetTimeStep());
//
//
//        ElemList actSDList(ptgrid_ );
//        actSDList.SetRegion( subdoms_[actSD] );
//        
//        EntityIterator it = actSDList.GetIterator();
//        for ( it.Begin(); !it.IsEnd(); it++ ) {
//
//          //ptElem=it.GetElem()->ptElem;
//          const UInt nrNodes  = (it.GetElem()->connect).GetSize();
//          dofs_ = Elem::shapes[it.GetElem()->type].dim;
//          //dofs_ = ptElem->GetDim();
//          
//          resultStressVector.Resize(nrNodes * dofs_);
//          rhsvec.Resize(nrNodes * dofs_);
//          rhsvec.Init();
//              
//          connecth=it.GetElem()->connect;
//          
//          rhsViscoMat->CalcElementMatrix(elemmat,it,it);
//          
//          // map connect to PDE node numbers
//          StdVector<Integer> connect_PDE;
//          eqnMap_->GetNodeEqn(connecth, connect_PDE);
//
//          /*
//          // output matrix with which BDF is computed
//          (*debug) << "fractional Damping matrix of Element" << std::endl;
//          (*debug) << elemmat << std::endl;
//          (*debug) << "actStep_=" << actStep_ << std::endl;
//          (*debug) << "numValues_=" << numValues_ << std::endl;
//           */
//              
//          for (UInt i=1; i<=numValues_; i++) {
//            if ( solMemoryVal_[i-1] == trueVAL )
//              {
//                GetElemSolution(solMemory_[i-1], elemsol, connect_PDE);
//                rhsvec += elemsol * coeff_[i];	     
//
//                // (*debug) << "elemsolOld" << std::endl;
//                // (*debug) << elemsol << std::endl;
//              }
//            else { // see NewmarkFracDamp
//            }		
//          }
//              
//          rhsAssemble = -(elemmat  * rhsvec *timeStepPowerFracDeriv_);
//          
//          // compute stress history from displacement history with the fractional constitutive equation
//          // later perhabs  this method should be implemented in a own class
//          
//          // calculates the stress history from the displacement history	  
//          //CalcStressHistoryForElement( connect_PDE, ptElem->GetNumNodes(),ptCoord,mymaterialData[actSD],ptElem,el);
//          // calculate linear stresses for stress history
//          //	  CalcStressHist( ptCoord, mymaterialData[actSD],connect_PDE,ptElem);
//          
//          fracDerivStressVec_.Init();
//          resultStressVector.Init();
//          
//          for (UInt i=1; i<numValues_; i++) {	     
//            Integer localElemNum = eqnMap_->Mesh2PdeElem(it.GetElem()->elemNum); 
//            GetStressVector(stressVector,localElemNum,i-1);
//            fracDerivStressVec_ += stressVector * coeff_[i+1];
//          }
//              
//          bdInt->calcElementVector(resultStressVector,it,fracDerivStressVec_);
//          
//          //                   double t1,t2;
//          //                   t1=0;
//          //                   t2=0;
//              
//          //                   for (UInt i=0;i<rhsAssemble.GetSize();i++)
//          //                     {
//          //                       t1 += rhsAssemble[i];
//          //                       t2 += resultStressVector[i];
//          //                     }
//
//          /*
//          (*debug) << "actStep_=" << actStep_ << std::endl 
//                   << "rhsAssemble=" << std::endl << rhsAssemble << std::endl
//                   << "resultStressVector=" << std::endl << resultStressVector << std::endl;
//          */         
//          if(model== "3param")
//            rhsAssemble  = rhsAssemble +  resultStressVector;
//          else if (model=="KelvinVoigt")
//          {
//            assert(false);
//            rhsAssemble  = rhsAssemble; // +  resultStressVector;
//          }
//          else
//            std::cerr << "unknown model for fractional damping" << std::endl;
//          //(*debug) <<  "rhs vector of timestep " << actStep_ << std::endl;
//          //(*debug) << rhsvec << std::endl;
//          //assemble to RHS
//          algsys_->SetElementRHS( rhsAssemble, pdeId_, connect_PDE );
//        }
//      }
//    }
  }
  

  void NewmarkFracDampMech::Corrector(Vector<Double>& solnew)
  {
    EXCEPTION("Reimplement after refactoring!")
//
//    solDeriv_vec_[SECOND_DERIV] = (solnew - solpred_) * a2_;
//    solDeriv_vec_[FIRST_DERIV] = solderiv1pred_ + solDeriv_vec_[SECOND_DERIV]*a3_;
//
//    Integer numEQNs = eqnMap_->GetNumEqns();
//    StdVector<UInt> connecth, connect_PDE;
//    Vector<Double> stressVec, displacementVec;
//    Matrix<Double> ptCoord;
//
//    stressVec.Resize(getDim());
//    stressVec.Init();
//    displacementVec.Resize(numEQNs * dofs_);  
//
//    std::map<RegionIdType, BaseMaterial*> mymaterialData;
//    BaseFE         * ptElem;
//
//    mymaterialData = ptStdPDE_->getPDEMaterialData();
//
//    if( (actStep_ % modulo_) == 0) {
//      for ( Integer i=fracMemory_-1; i>=1; i-- ) {
//        solMemory_[i] = solMemory_[i-1];
//        stressHistoryEl_[i] = stressHistoryEl_[i-1];
//      }
//      solMemory_[0] = solnew; //- solpred_;
//      
//      if ((actStep_ / modulo_)  <=  fracMemory_) {
//        solMemoryVal_.insert( solMemoryVal_.begin(),trueVAL);
//        // when solMemoryVal_ reaches size fracMemory_ , all entries are trueVAL
//        //  and stay the same for all time for
//      }
//      
//      for ( UInt actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {
//        StdVector<Elem*> elemssd;
//        ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
//        //ptgrid_->GetElemSD(elemssd,subdoms_[actSD]);
//          
//        for (UInt el=0; el < elemssd.GetSize(); el++) {
//          ptElem=elemssd[el]->ptElem;
//          connecth=elemssd[el]->connect;
//          ptgrid_->GetElemNodesCoord(ptCoord,connecth);
//          StdVector<Integer> connect_PDE;
//          eqnMap_->GetNodeEqn(connecth, connect_PDE);
//          
//          GetElemSolution(solMemory_[0], displacementVec, connect_PDE);
//          
//          CalcStress( ptElem, mymaterialData[subdoms_[actSD]], connect_PDE, ptCoord,  
//		      stressVec, displacementVec, elemssd[el]->elemNum);
//          Integer localElemNum = eqnMap_->Mesh2PdeNode( elemssd[el]->elemNum );
//          InsertStressVector(stressVec,localElemNum,0);		
//        }
//      }
//    }
  }

  void NewmarkFracDampMech::CalcParameters(Double dt)
  {
     REFACTOR;
//
//    //for predictors
//    a0_ = 0.5*(1-2.0*beta_)*dt_*dt_;
//    a1_ = (1-gamma_)*dt_;
//
//    //for correctors, matrices
//    a2_ = 1.0/(beta_*dt_*dt_);
//    a3_ = gamma_*dt_;
//
//    //for RHS, matrices
//    a4_ = gamma_ / (beta_*dt_);
  }


  void NewmarkFracDampMech::GetElemSolution ( const Vector<Double>& sol, 
                                              Vector<Double>& elemsol,
                                              const StdVector<Integer> & connectPDE ) {
    REFACTOR;
//   
//    elemsol.Resize(connectPDE.GetSize());
//    elemsol.Init();
//
//    for (UInt eqn=0; eqn<connectPDE.GetSize(); eqn++) {
//      if (connectPDE[eqn]>=1)
//        elemsol[eqn] = sol[connectPDE[eqn]-1];
//    }
  }

  void NewmarkFracDampMech::GLWeights(UInt memory, Double y ) {

   REFACTOR;
//    // reserve memory for weights of BDF, order of derivative is y-1
//    coeff_.resize(memory+1);
//    coeff_[0] = 1.0; // Index 0
//
//    // (*debug) << "coeff_" <<std::endl;
//    // (*debug) << coeff_[0]<<std::endl;
//
//    for (UInt i=1; i <= memory; i++) { // Index 1 .. memory 
//      coeff_[i] = coeff_[i-1] * (i-1-y)/(i);
//
//      // (*debug) << coeff_[i] << "   " << coeff_[i-1] * (i-1-y)/(i) <<std::endl;
//
//    }
  }

void NewmarkFracDampMech::CalcStress(BaseFE * aptelem, 
				     BaseMaterial* matDa, 
				     StdVector<Integer> connect_PDE, 
				     Matrix<Double> & ptCoord, 
				     Vector<Double> & stressVector,
				      Vector<Double> &displacementVector, 
				     Integer elemNr)
{
  
  Matrix<Double> bMat;
  Matrix<Double> dMat;
  Matrix<Double> alphaMat;
  Matrix<Double> betaMat;
  Matrix<Double> aMat;


  //double timeStep = GetTimeStep();

  //Integer dimStressVec = getDim();
  //Integer numEQNs = ptEQN_->GetNumEQNs();

  //const Integer nrNodes  = aptelem->GetNumNodes();

  //transform the type
  SubTensorType type;
  String2Enum(subType_,type);

  EXCEPTION( "Not working at the moment (-> contact Andreas)" );
  
//   BDBInt * viscoIntegrator = new LinViscoElastInt(matDa, type,
//                                                 "modifiedStiffness",timeStep );

//   Vector<Double> fracDerivStress, fracDerivDisplacement;
//   Vector<Double> term1, term2,term3;

//   fracDerivStress.Resize(dimStressVec);
//   stressVector.Resize(dimStressVec);
//   fracDerivDisplacement.Resize(nrNodes * dofs_);
//   term1.Resize(dimStressVec);
//   term2.Resize(dimStressVec);
//   term3.Resize(dimStressVec);

//   fracDerivStress.Init();
//   fracDerivDisplacement.Init();

//   viscoIntegrator->GetDMat(dMat);
//   viscoIntegrator->GetBMat(bMat,ptCoord);
//   GetBetaMat(betaMat,elastModule_, matDa);
//   GetAlphaMat(alphaMat);
//   GetAMat(aMat);

//   //computation of the term1 (actual displacement * material, damping factors)
//   //term1 = (dMat * bMat) * displacementVector;
//   Matrix<Double> helpMat;
//   helpMat = dMat * bMat;
//   term1 = helpMat * displacementVector;

//   //computation of term2 (fractional Derivative of displacement)
//   // solMemory[0] is the displacement value of the previous time step -> loop from 0
//   for(UInt k=0; k< numValues_;k++) {
//     GetElemSolution(solMemory_[k], displacementVector, connect_PDE);	
//     fracDerivDisplacement += displacementVector *  coeff_[k+1];   
//   }

//   //term2 = aMat * betaMat * bMat * fracDerivDisplacement;
//   Matrix<Double> helpMat2;
//   helpMat = betaMat * bMat;
//   helpMat2 = aMat * helpMat;
//   term2 = helpMat2 * fracDerivDisplacement;

//   //computation of term3 (fractional derivative of stress)
//   for(UInt k=0; k< numValues_;k++) {
//     Integer locElemNum = eqnMap_->Mesh2PdeNode( elemNr );
//     GetStressVector(stressVector, localElemNum,k);	
//     fracDerivStress += stressVector *  coeff_[k+1];   
//   }

//   //term3 =    aMat * alphaMat * fracDerivStress;
//   helpMat = aMat * alphaMat;
//   term3 = helpMat * fracDerivStress;


//   double t1,t2,t3;
//     t1=0;
//     t2=0;
//     t3=0;

//     for (UInt i=0;i<term1.GetSize();i++) {
//       t1 += term1[i];
//       t2 += term2[i];
//       t3 += term3[i];
//     }
    
//     stressVector = term1 + term2 - term3;  
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

     double val = 0.0;
     
     // compute matrix A,  same entries, 
     aMat.Resize(getDim(),true);
     aMat.Init();
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

     double val = 0.0;
     
     // compute matrix A,  same entries, 
     alphaMat.Resize(getDim(),true);
     alphaMat.Init();
     // set entries on the diagonal
      val = (timeStepPowerFracDeriv_)  * dampAlpha_ ;
	
     // set the emtries on the diagonal matrix, to get the inverse, 
     // the value on the diagonal are 1/val
     for(UInt i=0;i<getDimD();i++) {
       alphaMat.SetEntry(i,i,val);
     }  
   }


void NewmarkFracDampMech::GetBetaMat(Matrix<Double>& betaMat,  Double E, 
				     BaseMaterial* matData)
   { 

     double val = 0.0;

     //transform the type
    SubTensorType tensorType;
    String2Enum(subType_,tensorType);

    BDBInt * matMat = new linElastInt(matData, tensorType);
    matMat->calcDMat(betaMat);

    val = (dampBeta_/(E)) * timeStepPowerFracDeriv_;
    betaMat  = (betaMat *  val);
   }

  
  UInt NewmarkFracDampMech::getStressDim() {


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
      EXCEPTION("wrong subType(axi,planeStrain,3d) specified" );
      return 0;
    }
  }
  
  UInt NewmarkFracDampMech::getDim() {
    
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
      EXCEPTION("wrong subType(axi,planeStrein,3d) specified" );
      return 0;
    }
  }
} // end of namespace
