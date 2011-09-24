// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CURLCURLEDGE
#define FILE_CURLCURLEDGE

#include "baseForm.hh"

namespace CoupledField
{

// forward class declaration
class FeHCurl;

  /// Class for calculation of curl curl of edge elements
class CurlCurlEdgeInt : public BaseForm
{
public:
  /// Constructor
  CurlCurlEdgeInt( BaseMaterial* matData, bool coordUpdate = false );

  /// 
  virtual ~CurlCurlEdgeInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );

  /// Calculate curl derivative matrix
  void CalcBMat( Matrix<Double>& bMat, 
                 LocPointMapped& lp, FeHCurl* ptFE );

  /// calculates the curl, if the global derivates are already given in shapeDeriv
  /*!
    \param curl (output) Matrix with curl of edge shape functions
    \f[ \left( \begin{array}{c} E_{z,y} - E_{y,z}\\
    E_{x,z} - E_{z,x} \\
    E_{x,z} - E_{z,x} \\
    \end{array}\right) \f]
    \param shapeDeriv (input) Vector of matrices holding global derivates of edge shape functions
  */
  void CalcEdgeCurl(Matrix<Double>& curl, 
                    const StdVector< Matrix<Double> >& shapeDeriv);
  //! Calcualte the curl of an element solution
  void ApplyBMat( Vector<Double>& retVec,  
                  LocPointMapped& lp, BaseFE* ptFE, 
                  const Vector<Double>& solVec );
  
  //! Calcualte the curl of an element solution
    void ApplyBMat( Vector<Complex>& retVec,  
                    LocPointMapped& lp, BaseFE* ptFE, 
                    const Vector<Complex>& solVec );
    
    //! (De)Activate logging of average flux density (only for non-linear case)
    void SetLogging(bool doLogging ) { logging_ = doLogging;}

protected:
  UInt nrDofs_;

  //! multiplicative value for curl-curl operator 
  Double matVal_;

  //! true, if we have an orthotropic material
  bool isOrthotropic_;

  //! contains the orthotropic reluctivities
  Vector<Double> reluctivityVec_;

  //! Flag, if logging should be activated
  bool logging_;
};

}

#endif
