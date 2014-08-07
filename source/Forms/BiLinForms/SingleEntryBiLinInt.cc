#include "SingleEntryBiLinInt.hh"
#include "Domain/Domain.hh"


namespace CoupledField {


SingleEntryBiLinInt::SingleEntryBiLinInt( UInt numDofs, PtrCoefFct& val )
    : BiLinearForm() {

    name_ = "SingleEntryBiLinInt";
    numDofs_ = numDofs;
    
    // check, if we have a constant expression coefficient function
    if( val->GetDependency() == CoefFunction::GENERAL) {
      EXCEPTION("SingleEntryInt only works with constant coefficients");
    }
    val_ = val;

  }

SingleEntryBiLinInt::SingleEntryBiLinInt(  UInt numDofs, const std::string& val, 
                                           UInt dof, MathParser* mp ) 
   : BiLinearForm() {
  
  name_ = "SingleEntryBiLinInt";
  numDofs_ = numDofs;
  
  // assemble coefficient function
  StdVector<std::string> realVals(numDofs_);
  realVals.Init("0.0");
  realVals[dof] = val;
  val_ = CoefFunction::Generate(mp, Global::REAL, realVals);
  
  // check, if we have a constant expression coefficient function
  if( val_->GetDependency() == CoefFunction::GENERAL) {
    EXCEPTION("SingleEntryInt only works with constant coefficients");
  }
}

SingleEntryBiLinInt::SingleEntryBiLinInt( UInt numDofs, const std::string& real, 
                                          const std::string& imag, UInt dof,
                                          MathParser* mp ) 
: BiLinearForm() {
  
  name_ = "SingleEntryBiLinInt";
  numDofs_ = numDofs;
  
  // assemble coefficient function
   StdVector<std::string> realVals(numDofs_), imagVals(numDofs_);
   realVals.Init("0.0");
   realVals[dof] = real;
   
   imagVals.Init("0.0");
   imagVals[dof] = imag;
   val_ = CoefFunction::Generate(mp, Global::REAL, realVals, imagVals);
   
   // check, if we have a constant expression coefficient function
   if( val_->GetDependency() == CoefFunction::GENERAL) {
     EXCEPTION("SingleEntryInt only works with constant coefficients");
   }
}


SingleEntryBiLinInt::SingleEntryBiLinInt::~SingleEntryBiLinInt() {
}
  

void SingleEntryBiLinInt::CalcElementMatrix( Matrix<Double>& stiffMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2) {
      
  // we use just a dummy local point, as we assume constant
  // expression coefficient function
  LocPointMapped lpm;
  
  if( val_->GetDimType() == CoefFunction::SCALAR) {
    stiffMat.Resize(1, 1);
    val_->GetScalar(stiffMat[0][0], lpm);
  } else  if( val_->GetDimType() == CoefFunction::VECTOR) {
    Vector<Double> elemVec;
    val_->GetVector(elemVec, lpm);
    stiffMat.Resize(numDofs_, numDofs_);
    stiffMat.Init();
    for( UInt i = 0; i < numDofs_; ++ i ) {
      stiffMat[i][i] = elemVec[i];
    }
  } else {
    EXCEPTION( "SingleEntryInt only works for SCALAR and VECTOR" );
  }
}

void SingleEntryBiLinInt::CalcElementMatrix( Matrix<Complex>& stiffMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2) {
  // we use just a dummy local point, as we assume constant
  // expression coefficient function
  LocPointMapped lpm;
  
  if( val_->GetDimType() == CoefFunction::SCALAR) {
    stiffMat.Resize(1, 1);
    val_->GetScalar(stiffMat[0][0], lpm);
  } else  if( val_->GetDimType() == CoefFunction::VECTOR) {
    Vector<Complex> elemVec;
    val_->GetVector(elemVec, lpm);
    stiffMat.Resize(numDofs_, numDofs_);
    stiffMat.Init();
    for( UInt i = 0; i < numDofs_; ++ i ) {
      stiffMat[i][i] = elemVec[i];
    }
  } else {
    EXCEPTION( "SingleEntryInt only works for SCALAR and VECTOR" );
  }
}


//  void SingleEntryBiLinInt::CalcElemVector( Vector<Double>& elemVec,
//                                       EntityIterator& ent1) {
//    
//    // we use just a dummy local point, as we assume constant
//    // expression coefficient function
//    LocPointMapped lpm; 
//    if( val_->GetDimType() == CoefFunction::SCALAR) {
//      elemVec.Resize(1);
//      val_->GetScalar(elemVec[0], lpm);
//    } else  if( val_->GetDimType() == CoefFunction::VECTOR) {
//      val_->GetVector(elemVec, lpm);
//    } else {
//      EXCEPTION( "SingleEntryInt only works for SCALAR and VECTOR" );
//    }
//  }
//  
//  void SingleEntryBiLinInt::CalcElemVector( Vector<Complex>& elemVec,
//                                       EntityIterator& ent1) {
//    
//    // we use just a dummy local point, as we assume constant
//    // expression coefficient function
//    LocPointMapped lpm; 
//    if( val_->GetDimType() == CoefFunction::SCALAR) {
//      elemVec.Resize(1);
//      val_->GetScalar(elemVec[0], lpm);
//    } else  if( val_->GetDimType() == CoefFunction::VECTOR) {
//      val_->GetVector(elemVec, lpm);
//    } else {
//      EXCEPTION( "SingleEntryInt only works for SCALAR and VECTOR" );
//    }
//  }

  
}
