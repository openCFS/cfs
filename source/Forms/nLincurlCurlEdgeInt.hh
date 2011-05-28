// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINCURLCURLEDGE
#define FILE_NLINCURLCURLEDGE

#include <fstream>
#include "curlCurlEdgeInt.hh"

namespace CoupledField {

  /// Class for calculation of nonlinear curl curl of edge elements
class nLinCurlCurlEdgeInt : public CurlCurlEdgeInt
{
public:
  /// Constructor
  nLinCurlCurlEdgeInt( BaseMaterial* matData, bool coordUpdate = false );

  /// 
  virtual ~nLinCurlCurlEdgeInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
  
  //! sets type of nonlinear algorithm
  virtual void SetNonLinMethod(NonLinMethodType atype);
  
  //! (De)Activate logging of average flux density
  void SetLogging(bool doLogging ) { logging_ = doLogging;} 

private:
  Vector<Double> magPot_;    //!< magnetic vector potential at nodes
  NonLinMethodType nonLinType_;  //!< type of nonlinear algorithm

  //! contains the current orthotropic reluctivities
  Vector<Double> currReluctivityVec_;  

  //! contains the current derivatives oforthotropic reluctivities
  Vector<Double> currDerivReluctivityVec_;
  
  // Temporary section
  //! Contains element numbers for which we want to log average flux density
  StdVector<UInt> fluxElems_;
  
};

} // end of namespace

#endif
