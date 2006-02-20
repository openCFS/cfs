#include <fstream>
#include <iostream>
#include <math.h>

#include "tools.hh"
#include "Matrix/matrix.hh"
#include "Elements/elements_header.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"


namespace CoupledField {
 
  // =======================
  //   Issue Error Message
  // =======================
  void Error( const Char * Text, const Char * const filename,
              const UInt numline ) {
    Info->Error( Text, filename, numline );
  }

  // ====================================================
  //   Issue Error Message (using global string stream)
  // ====================================================
  void Error( const Char *const filename, const UInt numline ) {

    // Obtain error message and clear string stream
    std::string errMsg = "";
    if ( error != NULL ) {
      errMsg = error->str();
      error->str( "" );
    }
    else {
      errMsg = "<Caller did not provide error description>";
    }

    // Delegate work to WriteInfo::Error()
    Info->Error( errMsg.c_str(), filename, numline );
  }

  // =========================
  //   Issue Warning Message
  // =========================
  void Warning( const Char *Text, const Char * const filename,
                const UInt numline ) {
    Info->Warning( Text, filename, numline );
    //std::cerr << "\033[31mWARNING:\033[0m " << Text;
  }

  // ======================================================
  //   Issue Warning Message (using global string stream)
  // ======================================================
  void Warning( const Char *const filename, const UInt numline ) {

    // Obtain error message and clear string stream
    std::string warnMsg = "";
    if ( warning != NULL ) {
      warnMsg = warning->str();
      warning->str( "" );
    }
    else {
      warnMsg = "<Caller did not provide error description>";
    }

    // Delegate work to WriteInfo::Warning()
    Info->Warning( warnMsg.c_str(), filename, numline );
  }

  // =========================================
  //   Split a string into a list of strings
  // =========================================
  void SplitStringList( std::string list, StdVector<std::string> &strVec,
                        Char delimiter ) {

    UInt lastDelim = 0;
    strVec.Clear();
    UInt i=0;
  
    // ignore all leading spaces
    while ( i < list.length() && list[i] == ' ' ) {
      i++;
      lastDelim++;
    }

    // get the n-1 entries of the list
    for ( i = 0; i < list.length(); i++ ) {
      if ( list[i] == delimiter ) {
        strVec.Push_back( std::string( list, lastDelim, i-lastDelim ));
        i++;

        // ignore next spaces
        while ( i < list.length() && list[i] == ' ' ) {
          i++;
        }

        lastDelim = i;
      }
    }

    // get the n-th entry
    i = lastDelim;
    while ( i < list.length() && list[i] != ' ' ) {
      i++;
    }

    strVec.Push_back( std::string( list, lastDelim, i-lastDelim ));

  }


  template<UInt dim>
  Point<dim> & Point<dim>::operator+(const Point<dim>&t) {
    UInt i;
    for (i=0; i<dim; i++)
      p[i]+=t.p[i];     
    return *this;
  }


  template<UInt dim>
  void PrintPoint(Point<dim> point, std::ostream * out) {
    UInt i;
    for(i=0; i<dim; i++) 
      (*out)<< "   " << point[i];
  }

  template<UInt dim>
  Point<dim> & Point<dim>::operator=(const Point<dim>&t) {
    UInt i;
    for (i=0; i<dim; i++)
      p[i]=t.p[i];
    return *this;
  }

  template<UInt dim>
  Point<dim> & Point<dim>::operator-(const Point<dim>&t) {
    UInt i;
    for (i=0; i<dim; i++)
      p[i]-=t.p[i];     
    return *this;
  }

  Double Sin(const Double x) {
    return sin(x);
  }

  Double Cos(const Double x) {
    return cos(x);
  }

  pfn1var FncReader(const std::string namefnc) {
    if (namefnc == "sin") return Sin;
    else {
      if (namefnc == "cos") {
        return Cos;
      }
      else {
        (*error) << "Function '" << namefnc
                 << "' is not predefined in CFS++";
        Error( __FILE__, __LINE__ );
      }
    }
  }

  template<UInt dim>
  Double dist(Point<dim> a,Point<dim> b) {
    Double preSqrt=0;
    UInt i;
    for (i=0; i<dim; i++)
      preSqrt+=sqr(a[i]-b[i]);
    return sqrt(preSqrt);
  }

  Double dist_Mat(Matrix<Double> a) {
    Double preSqrt=0;
    UInt i;
    UInt k=a.GetSizeRow();
    // std::cout<<"tools.cc:size of matrix: "<<k<<std::endl;
    for (i=0; i<k; i++)
      preSqrt+=sqr(a[i][0]-a[i][1]);
    return sqrt(preSqrt);
  }

#ifdef __GNUC__
  template class Point<2>;
  template class Point<3>;

  template void PrintPoint(Point<2>, std::ostream *);
  template void PrintPoint(Point<3>, std::ostream *);

  template Double dist(Point<2>,Point<2>);
  template Double dist(Point<3>,Point<3>);
#endif


  /// a-->b
  void calcNormal2Line(Vector<Double> & normal,Point<2> a,Point<2> b) {
    Double distance=dist(a,b);
    normal[0]=(b[1]-a[1])/distance;
    normal[1]=(a[0]-b[0])/distance;
  }

  /// a-->b
  void calcNormal2Line_Mat(Vector<Double> & normal,Matrix<Double> a)
  {
    Double distance=dist_Mat(a);
    // std::cout<<"distance :"<<distance<<std::endl;
  
    normal[0]=(a[1][1]-a[1][0])/distance;
    normal[1]=(a[0][1]-a[0][0])/distance;
  }


  /// a-->b-->c. no fix orientation.
  void calcNormal2Surface(Vector<Double> & normal,Point<3> a,Point<3> b,
                          Point<3> c) {
    Point<3> t,s;
    Double L2_normal; 
    s=(a-b);
    t=(c-b);
    normal[0]=t[1]*s[2]-t[2]*s[1];
    normal[1]=t[2]*s[0]-t[0]*s[2];  
    normal[2]=t[0]*s[1]-t[1]*s[0];
    L2_normal=sqrt(sqr(normal[0])+sqr(normal[1])+sqr(normal[2]));
    normal[0]=normal[0]/L2_normal;
    normal[1]=normal[1]/L2_normal;
    normal[2]=normal[2]/L2_normal;  
  }

  /// a-->b-->c. no fix orientation.
  void calcNormal2Surface_Mat(Vector<Double> & normal,Matrix<Double> ptCoord){
    Point<3> t,s,a,b,c;

    for (UInt i=0;i<3;i++) {
      a[i]=ptCoord[i][0];
      b[i]=ptCoord[i][1];      
      c[i]=ptCoord[i][2];
    }
  
    Double L2_normal; 
    s=(a-b);
    t=(c-b);
    normal[0]=t[1]*s[2]-t[2]*s[1];
    normal[1]=t[2]*s[0]-t[0]*s[2];  
    normal[2]=t[0]*s[1]-t[1]*s[0];
    L2_normal=sqrt(sqr(normal[0])+sqr(normal[1])+sqr(normal[2]));
    normal[0]=normal[0]/L2_normal;
    normal[1]=normal[1]/L2_normal;
    normal[2]=normal[2]/L2_normal;  
  }

  char * c_string(const std::string & s) {
    char * p = new char[s.length()+1];
    s.copy(p, std::string::npos);
    p[s.length()]=0;
    return p;
  }

  UInt defineRefinements(const Double tolElem, const Double tolTotal,
                         const UInt noOfChilds){
    Double tmp = log( tolElem/tolTotal ) / log( (Double) noOfChilds );
    return (UInt)tmp + 1;
  }

  Double CalcArea(Elem * ptE, Grid * ptgrid, const UInt level){
    Double         area = 0;
    BaseFE         * ptelem = ptE->ptElem;
    const StdVector<UInt> & connect = ptE->connect;
    Matrix<Double> ptCoord;
    UInt        nrIntPnts = ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    UInt        i;
    Double         jacDet;

    ptgrid->GetElemNodesCoord(ptCoord,connect);

    for (i=0; i<nrIntPnts; i++) {
      jacDet = ptelem->CalcJacobianDetAtIp(i+1, ptCoord);
      area +=jacDet*intWeights[i];
    }

    return area;
  }


  // **************** Conversion Function *****************
  
  Double String2Double( const std::string & val) {

   
    Double retVal = 0.0;
    char *endp;
    retVal = strtod(val.c_str(), &endp);
    if (! (val.c_str() != endp && *endp == '\0') ) {
      (*error) << "Could not convert '" << val << "' into a double!";
      Error( __FILE__, __LINE__ );
    }
    return retVal;
    
  }

  Integer String2Int( const std::string & val) {
    Integer retVal = 0;
    char *endp;
    retVal = strtol(val.c_str(), &endp, 0);
    if (! (val.c_str() != endp && *endp == '\0') ) {
      (*error) << "Could not convert '" << val 
               << "' into an integer!";
      Error( __FILE__, __LINE__ );
    }
    return retVal;
  }

  UInt String2UInt( const std::string & val) {

    UInt retVal=0;
    char *endp;
    retVal = strtoul(val.c_str(), &endp, 0);
    if (! (val.c_str() != endp && *endp == '\0') ) {
      (*error) << "Could not convert '" << val 
               << "' into an unsigned integer!";
      Error( __FILE__, __LINE__ );
    }
    return retVal;
  }

   //! Converts a string vector into a double vector
  void String2Double( StdVector<Double> & retVal, 
                      const StdVector<std::string> & val ) { 
    retVal.Resize( val.GetSize() );
    
    for ( UInt i = 0; i < val.GetSize(); i++ ) {
      if ( std::sscanf ( val[i].c_str(), "%lf", &retVal[i]) != 1 ) {
        (*error) << "Could not convert '" << val 
                 << "' into a double!";
        Error( __FILE__, __LINE__ );
      }
      
    }
  }
    

  //! Converts a double vector into a string vector
  void Double2String( StdVector<std::string> & retVal, 
                      const StdVector<Double> & val ) {

    retVal.Resize( val.GetSize() );
    std::ostringstream mystream;

    for ( UInt i = 0; i < val.GetSize(); i++ ) {
      mystream << val[i];
      retVal[i] = mystream.str();
      mystream.str("");
    }
  }
  
  //! Converts a string vector into an integer vector
  void String2Int( StdVector<Integer> & retVal, 
                   const StdVector<std::string> & val ) {
    retVal.Resize( val.GetSize() );
    
    for ( UInt i = 0; i < val.GetSize(); i++ ) {
      if ( std::sscanf ( val[i].c_str(), "%d", &retVal[i]) != 1 ) {
        (*error) << "Could not convert '" << val 
                 << "' into an integer!";
        Error( __FILE__, __LINE__ );
      }
      
    }
  }
  
  //! Converts an integer vector into a string vector
  void Int2String( StdVector<std::string> & retVal, 
                   const StdVector<Integer> & val ) {
    
    retVal.Resize( val.GetSize() );
    std::ostringstream mystream;
    
    for ( UInt i = 0; i < val.GetSize(); i++ ) {
      mystream << val[i];
      retVal[i] = mystream.str();
      mystream.str("");
    }
  }
  
  //! Converts a string vector into an unsigned integer vector
  void String2UInt( StdVector<UInt> & retVal, 
                    const StdVector<std::string> & val ) {
    retVal.Resize( val.GetSize() );
    
    for ( UInt i = 0; i < val.GetSize(); i++ ) {
      if ( std::sscanf ( val[i].c_str(), "%u", &retVal[i]) != 1 ) {
        (*error) << "Could not convert '" << val 
                 << "' into an unsigned integer!";
        Error( __FILE__, __LINE__ );
      }
    }
  }
  
  //! Converts an unsigned integer vector into a string vector
  void UInt2String( StdVector<std::string> & retVal, 
                    const StdVector<UInt> & val ) {
    retVal.Resize( val.GetSize() );
    std::ostringstream mystream;
    
    for ( UInt i = 0; i < val.GetSize(); i++ ) {
      mystream << val[i];
      retVal[i] = mystream.str();
      mystream.str("");
    }
  }

  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Point<2>
#pragma instantiate Point<3>
#pragma instantiate void PrintPoint(Point<2>, std::ostream *)
#pragma instantiate void PrintPoint(Point<3>, std::ostream *)
#endif
}// namespace CoupledField
