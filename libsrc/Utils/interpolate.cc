#include "interpolate.hh"

namespace CoupledField {
  
  // define static maps
  std::map<std::string, Vector<Double> > Interpolate1D::xVals_ = 
    std::map<std::string, Vector<Double> >();
  std::map<std::string, Vector<Double> > Interpolate1D::yVals_ = 
    std::map<std::string, Vector<Double> >();
  
  Double Interpolate1D::Interpolate( const char* fileName, 
                                     double entry,
                                     double method ) {
    // check if file was already read in
    if( xVals_.find(std::string(fileName)) == xVals_.end() ) {
      Vector<Double> xValsTemp, yValsTemp;
      ReadFile( fileName, xValsTemp, yValsTemp );
      xVals_[fileName] = xValsTemp;
      yVals_[fileName] = yValsTemp;
    }

    // get interpolation data
    Vector<Double> const & xVals = xVals_[fileName];
    Vector<Double> const & yVals = yVals_[fileName];    
    

    // Convert type of interpolation


    InterpolType type = InterpolType(method);
    
    
    // Perform interpolation
    Double val = 0.0;
    
    switch ( type ) {
      
    case N_NEIGHBOUR:
      std::cerr << "--> Nearest Neighbour Interpolation\n";
      break;
      
    case LINEAR:
      std::cerr << "--> Linear Interpolation\n";
      break;
      
    case CUBIC_SPLINE:
      std::cerr << "--> Cubic Spline Interpolation\n";
      break;
    
    default:
      *error << "Interpolation type " << type << " not known!";
      Error( __FILE__, __LINE__ );
    }

    return val;
    
  }
  
  void Interpolate1D::ReadFile( const char* fileName, 
                                Vector<Double>& xVals,
                                Vector<Double>& yVals ) {
    
    // check if file 'fileName' exists
    
    // open file
    
    // loop over all lines and split it up in (x,y) values
    
    // store values in vectors xVals and yVals
  }

}
