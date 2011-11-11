#ifndef CFS_INTERPOLATE_HH
#define CFS_INTERPOLATE_HH

#include <map>
#include <string>

#include "General/Environment.hh"
#include "MatVec/Vector.hh"

namespace CoupledField {
  
  //! Class for performin one-dimensional interpolation
  class Interpolate1D {
    
  public:
    
    //! Define available interpolation methods
    typedef enum {N_NEIGHBOUR = 0, 
                  LINEAR = 1, 
                  CUBIC_SPLINE= 2 } InterpolType;
    
    //! Calculate 1D interpolation using data specified in file

    //! Uses the (x,y) values specified in a file to interpolate the y-value
    //! for a given x-value \param entry. The method can be specified 
    //! according to the integer representation of InterpolType.
    //! \param fileName Name of the file which contains white-space separated
    //!                 (x,y) values
    //! \param entry x-value of the value to be interpolated
    //! \param method Integer representation of the interpolation method
    //!               (ref. InterpolType)
    //! \return Interpolated value
    static Double Interpolate( const char* fileName, 
                               Double xEntry, 
                               Double method );

  private:
    
    //! Read in data from file
    static void ReadFile( const char* fileName,
                          Vector<Double>& xVals,
                          Vector<Double>& yVals );

    //! Calculate second derivative for cubic spline interpolation
    static void Spline( const Vector<Double>& x, const Vector<Double>& y,
                        Double yp0, Double ypn,
                        Vector<Double>& y2 );
    
    //! Map filename <-> x-value vector
    static std::map<std::string, Vector<Double> > xVals_;

    //! Map filename <-> y-value vector
    static std::map<std::string, Vector<Double> > yVals_;
  };


} // end of namespace
#endif
