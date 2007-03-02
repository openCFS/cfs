#ifndef FILE_CURLCURLEDGE
#define FILE_CURLCURLEDGE

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation of curl curl of edge elements
class CurlCurlEdgeInt : public BaseForm
{
public:
  /// Constructor
  CurlCurlEdgeInt(Double laplVal);

  /// 
  virtual ~CurlCurlEdgeInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
  
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
                    const StdVector< Matrix<Double>* >& shapeDeriv);
  

private:
  /// multiplicative value for curl curl integration 
  Double reluctivity_;
};

}

#endif
