// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "interpolate.hh"

namespace CoupledField {
  
  // define static maps
  std::map<std::string, Vector<Double> > Interpolate1D::xVals_ = 
    std::map<std::string, Vector<Double> >();
  std::map<std::string, Vector<Double> > Interpolate1D::yVals_ = 
    std::map<std::string, Vector<Double> >();
  
  Double Interpolate1D::Interpolate( const char* fileName, 
                                     double xEntry,
                                     double method ) {

    // check if file was already read in
    if( xVals_.find(std::string(fileName)) == xVals_.end() ) {
      Vector<Double> xValsTemp, yValsTemp;
      ReadFile( fileName, xValsTemp, yValsTemp );
      xVals_[fileName] = xValsTemp;
      yVals_[fileName] = yValsTemp;
    }

    // get interpolation data,
    //  nescessary, when file was already read in
    Vector<Double> const & xVals = xVals_[fileName];
    Vector<Double> const & yVals = yVals_[fileName];

    // Convert type of interpolation
    InterpolType type = InterpolType(method);
    
    Double yValue = 0.0;

    // get index of last element
    const UInt kend = xVals.GetSize() - 1;

    // if coordinate is out of bounds, return boundary value (i.e.first or last) 
    if ( xEntry > xVals[kend] ) {
      yValue = yVals[kend];
    }
    else if ( xEntry < xVals[0] ) {
      yValue = yVals[0];
    }
    else {

      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      //  klo and khi bracket the input value of xEntry
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (xVals[k] > xEntry)
          khi=k;
        else
          klo=k;
      }

      // size of x interval
      Double dxVal=xVals[khi] - xVals[klo];

      // The x-values must be distinct!
      if (dxVal == 0.0) {
        EXCEPTION("You cannot have two equal x values!" );
      }

      // relative distance of xEntry to x-Value bounds
      Double a = ( xVals[khi] - xEntry )/dxVal;
      Double b = ( xEntry - xVals[klo] )/dxVal;

      // vector containing second derivative for cubic spline interpolation
      Vector<Double> y2Vals;

      // Perform interpolation    
      switch ( type ) {
        
      case N_NEIGHBOUR:
        if ( a <= 0.5 )
          yValue = yVals[khi];
        else
          yValue = yVals[klo];
        break;
        
      case LINEAR:
        yValue = a * yVals[klo] + b * yVals[khi];
        break;

      case CUBIC_SPLINE:
        Spline( xVals, yVals, 0.0, 0.0, y2Vals); // zero slope at endpoints
        yValue = a*yVals[klo] + b*yVals[khi] 
          + ( (a*a*a-a)*y2Vals[klo] + (b*b*b-b)*y2Vals[khi] )*(dxVal*dxVal)/6.0;
        break;

      default:
        EXCEPTION( "Interpolation type " << type << " not known!" );
      }
    }

    return yValue;
  }

  
  void Interpolate1D::ReadFile( const char* fileName, 
                                Vector<Double>& xVals,
                                Vector<Double>& yVals ) {
    
    // open file
    std::ifstream sampleData;
    sampleData.open( fileName );


    // check if file 'fileName' exists
    if ( !sampleData ) {
      EXCEPTION( "Failed to open file '" << fileName
                 << "' containing data of sample function!" );
    }
    sampleData.clear(); // clear flags

    // we don't trust .eof()
    sampleData.seekg(0,std::ios::end);
    std::string::size_type pos = 0, pos_end = sampleData.tellg(),
      line_end_pos , line_start_pos = 0;

    // start from the beginning
    sampleData.seekg(0,std::ios::beg);
    std::string buf;
    Double      x, y;
    // initialize vectors xVals,yVals
    xVals.Resize(0);
    yVals.Resize(0);
    
    // loop over all lines and split it up in (x,y) values
    while( pos <= pos_end ) {
      buf = "";
      line_start_pos = sampleData.tellg();
      std::getline(sampleData,buf,'\n');
      line_end_pos = sampleData.tellg();
        
      // big choice of signs for comment's
      if (buf.length() != 0 &&
          buf[0] != '#' &&
          buf[0] != '%' && 
          buf[0] != '!') {

        sampleData.seekg(line_start_pos); // rewind
        sampleData >> x >> y;                    

        // store values in vectors xVals and yVals
        xVals.Push_back(x);
        yVals.Push_back(y);

        sampleData.ignore(100,'\n');
      }

      pos = sampleData.tellg();  // and, where we are ?    

      std::string errMsg;          
      if( pos != line_end_pos) {
        EXCEPTION( "The sample data file '" << fileName
                   << "' is not correctly formatted.\n"
                   << "Please correct it!" );
      }
    }

    sampleData.close();
  }


  void Interpolate1D::Spline( const Vector<Double>& x, const Vector<Double>& y,
                              Double yp0, Double ypn,
                              Vector<Double>& y2) {
    
    // Taken from "Numerical Recipes in C", W.H. Press et al.
    
    // Given arrays x[0..n] and y[0..n] containing a tabulated function,
    // i.e., yi = f(xi), with x0 < x1 < .. . < xN, and given values yp0 and 
    // ypn for the first derivative of the interpolating function at 
    // points 0 and n, respectively, this routine returns an array y2[0..n] 
    // that contains the second derivatives of the interpolating function 
    // at the tabulated points xi. If yp0 and/or ypn are equal to 1 × 10^30 
    // or larger, the routine is signaled to set the corresponding boundary
    // condition for a natural spline, with zero second derivative on that
    // boundary.

    UInt n = x.GetSize() - 1; // last index
    Vector<Double> u;
    u.Resize(x.GetSize() );
    u.Init(0.0);
    y2.Resize(x.GetSize() );
    y2.Init(0.0);

    // The lower boundary condition is set either to be natural
    // or else to have a specified first derivative.
    if (yp0 > 0.99e30) {
      y2[0] = 0.0;
      u[0] = 0.0;
    }
    else {
      y2[0] = -0.5;
      u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp0);
    }

    // This is the decomposition loop of the tridiagonal algorithm.
    //  y2 and u are used for temporary storage of the decomposed factors.
    Double p,sig;
    for ( UInt i=1; i<=n-1; i++) {

      sig   = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p     = sig*y2[i-1]+2.0;
      y2[i] = (sig-1.0)/p;
      u[i]  = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]  = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
    }

    // The upper boundary condition is set either to be natural
    // or else to have a specified first derivative.
    Double un, qn;
    if (ypn > 0.99e30) {
      qn = 0.0;
      un = 0.0;
    }
    else {
      qn = 0.5;
      un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);


    // This is the backsubstitution loop of the tridiagonal algorithm.
    for ( Integer k=n-1; k>=0; k-- ) {
      y2[k] = y2[k]*y2[k+1]+u[k];
    }

  }

}
