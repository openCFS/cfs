#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "tools.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Exception.hh"
#include "Domain/Domain.hh"
#include "Utils/mathParser/mathParser.hh"

using boost::char_separator;
using boost::tokenizer;
using boost::bad_lexical_cast;

DECLARE_LOG(tools)
DEFINE_LOG(tools, "tools")

namespace CoupledField {

  // =========================================
  //   Split a string into a list of strings
  // =========================================
  void SplitStringList( const std::string &list, StdVector<std::string> &strVec,
                        const char delimiter ) {

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

  std::string ConvertToFilename(std::string org)
  {
    std::string result = org;
    boost::replace_all(result, ":", "_");
    boost::replace_all(result, ",", "_");
    boost::replace_all(result, "(", "");
    boost::replace_all(result, ")", "");
    return result;
  }

  void SplitStringListWhitespace(const std::string &s, StdVector<std::string> &strVec)
  {
    char_separator<char> sep(" ,\t\n");
    tokenizer<char_separator<char> > tok(s, sep);
    for(tokenizer<char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg)
      strVec.Push_back(*beg);
  }


  // =========================================================================
  //  COMPLEX CONVERSION
  // =========================================================================
  //! Convert (real,imag) => amplitude
  Double RealImagToAmpl( const Complex& in ) {
    return sqrt(in.real() * in.real() + in.imag() * in.imag());
  }

  //! Convert (real,imag) => phase
  Double RealImagToPhase( const Complex& in ) {
    return (fabs(in.imag()) > 1e-16) ?                   
        atan2(in.imag(),in.real() )*180/M_PI : 
        ( in.real() < Double(0.0) ) ? Double(180) : Double(0) ; 
  }


  //! Convert (ampl,phase) => (real,imag)
  Complex AmplPhaseToComplex( Double val, Double phase ) {
    return Complex( val * cos( phase / 180 * M_PI ),
                    val * sin( phase / 180 * M_PI ) ); 
  }

  //! Convert (ampl,phase) => real
  Double AmplPhaseToReal( Double val, Double phase ) {
    return val * cos( phase / 180 * M_PI );
  }

  //! Convert (ampl,phase) => imag
  Double AmplPhaseToImag( Double val, Double phase ) {
    return val * sin( phase / 180 * M_PI );
  }


  //! Convert (ampl,phase) => real (strings)
  std::string AmplPhaseToReal( const std::string& val, 
                               const std::string& phase,
							   bool isInRad ) {
    if( phase == "0.0" || phase == "0" ) {
      return "( " + val + " ) ";
    } else {
      if (isInRad)
        return "( (" + val + ") * cos( " + phase + " ) )";
      else
        return "( (" + val + ") * cos( " + phase + " / 180 * pi ) )";
    }
  }

  //! Convert (ampl,phase) => imag (strings)
  std::string AmplPhaseToImag( const std::string& val, 
                               const std::string& phase,
							   bool isInRad ) {
    if( phase == "0.0" || phase == "0" ) {
      return "( 0.0 )";
    } else {
      if (isInRad)
        return "( (" + val + ") * sin( " + phase + " ) )";
      else
        return "( (" + val + ") * sin( " + phase + " / 180 * pi ) )";
    }
  }

  // ------ vector versions -----
  //! Convert (ampl,phase) => (real,imag) (strings vectors)
  void AmplPhaseToRealImag( const StdVector<std::string>& val, 
                            const StdVector<std::string>& phase,
                            StdVector<std::string>& real,
                            StdVector<std::string>& imag ) {
    assert(val.GetSize() == phase.GetSize());
    real.Resize(val.GetSize());
    imag.Resize(val.GetSize());
    for( UInt i = 0; i < val.GetSize(); ++i ) {
      real[i] = AmplPhaseToReal(val[i], phase[i]);
      imag[i] = AmplPhaseToImag(val[i], phase[i]);
    }    
  }

  // -----------------------------------------------------------------------

  Double dist_Mat(const Matrix<Double> &a) {
    Double preSqrt=0;
    UInt i;
    const UInt k=a.GetNumRows();
    // std::cout<<"tools.cc:size of matrix: "<<k<<std::endl;
    for (i=0; i<k; i++)
      preSqrt+= (a[i][0]-a[i][1]) * (a[i][0]-a[i][1]);
    return sqrt(preSqrt);
  }


  std::string ToValidXML(const std::string& input)
  {
    std::string out = input;
    std::replace_if(out.begin(), out.end(), boost::is_any_of(", |"), '_');
    return out;
  }


  Double NormL2(const Double* data, const UInt size)
  {
    Double result = 0.0;
    for(UInt i = 0; i < size; i++)
      result += data[i] * data[i];

    return sqrt(result);
  }

  Double NormL2(const Double* data1, const Double* data2, const UInt size)
  {
    Double result = 0.0;
    for(UInt i = 0; i < size; i++)
      result += (data1[i] - data2[i]) * (data1[i] - data2[i]);

    return sqrt(result);
  }

  Double NormL2(const SingleVector* data, const SingleVector* data2)
  {
    assert(data != NULL && data2 != NULL);

    if(data->GetSize() != data2->GetSize()) EXCEPTION("incompatible sizes");
    if(data->GetEntryType() != data2->GetEntryType()) EXCEPTION("incompatible entry types");

    if(data->GetEntryType() == BaseMatrix::COMPLEX)
      return dynamic_cast<const Vector<Complex>& >(*data).NormL2(dynamic_cast<const Vector<Complex>& >(*data2));
    else
      return dynamic_cast<const Vector<Double>& >(*data).NormL2(dynamic_cast<const Vector<Double>& >(*data2));
  }


  template <class TYPE>
  std::string ToString(const TYPE* data, unsigned int size)
  {
    std::ostringstream os;

    os << "[";
    for(unsigned int i = 0; i < size; i += 1)
      os << data[i] << ((i < size-1) ? ", ": "]");

    return os.str();
  }

  template <class TYPE>
  std::string ToString(const StdVector<Vector<TYPE> >& data, bool new_line)
  {
    std::ostringstream os;

    for(unsigned int i = 0; i < data.GetSize(); i++)
      os << i << "=" << "[" << data[i].ToString() << "]" << (new_line ? "\n": " ");

    return os.str();
  }

  template <class TYPE>
  std::string ToString(const StdVector<StdVector<TYPE> >& data, bool new_line)
  {
    std::ostringstream os;

    for(unsigned int i = 0; i < data.GetSize(); i++)
      os << i << "=" << "[" << data[i].ToString() << "]" << (new_line ? "\n": " ");

    return os.str();
  }


  template <>
  std::string ToString<std::complex<Double> >(const std::complex<Double>* data, unsigned int size)
  {
    std::ostringstream os;

    // prepare for copy and paste to matlab
    os << "[";
    for(unsigned int i = 0; i < size; i += 1)
      os << data[i].real() << " + " << data[i].imag() << "i" << ((i < size-1) ? ", ": "]");

    return os.str();
  }


  Double Average(const Double* data, unsigned int size)
  {
    Double sum = 0.0;
    for(unsigned int i = 0; i < size; i++)
      sum += *(data + i);

    return sum / size;
  }

  Double StandardDeviation(const Double* data, unsigned int size)
  {
    // expected value
    Double e = Average(data, size);

    // variance
    Double v = 0.0;
    for(unsigned int i = 0; i < size; i++)
      v += (*(data + i) - e) * (*(data + i) - e);

    return v / size;

    return sqrt(v);
  }


  void Assign(Matrix<Double>& target, const Matrix<Double>& other, const Double factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Complex>& other, const Complex factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Double>& other, const Complex factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Double>& other, const Double factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }
  
  void Assign(Vector<Double>& target, const Vector<Double>& other, const Double factor)
   {
     UInt n = other.GetSize();
     target.Resize(n);
     for(UInt i = 0; i < n; i++)
         target[i] = factor * other[i];
   }
  
 
  void Assign(Vector<Complex>& target, const Vector<Complex>& other, const Double factor)
  {
    UInt n = other.GetSize();
    target.Resize(n);
    for(UInt i = 0; i < n; i++)
      target[i] = factor * other[i];
  }

  unsigned int SearchMinMax(const Matrix<double>& mat, unsigned int row, bool minimum, double* val, EigenInfo* info)
  {
    // make sure we have an info to operate on
    EigenInfo tmp;
    if(info == NULL)
      info = &tmp;

    info->col = 0;
    info->max = mat[row][0];
    info->min = mat[row][0];

    for(unsigned int c = 1, n = mat.GetNumCols(); c < n; c++)
    {
      double val = mat[row][c];
      if(val < info->min)
      {
        info->min = val;
        if(minimum)
          info->col = c;
      }
      if(val > info->max)
      {
        info->max = val;
        if(!minimum)
          info->col = c;
      }
    }

    if(val != NULL)
      *val = minimum ? info->min : info->max;

    return info->col;
  }

  void Conj(Matrix<Complex>& mat)
  {

    for(UInt r = 0; r < mat.GetNumRows(); r++)
      for(UInt c = 0; c < mat.GetNumCols(); c++)
        mat[r][c] = std::conj(mat[r][c]);
  }

  int MemoryUsage(bool peak)
  {
    // if the file dows not exist we return 0
    std::ifstream file("/proc/self/status", std::ifstream::in);

    std::string data;
    while(file)
    {
      file >> data;
      if(data == (peak ? "VmPeak:" : "VmSize:"))
      {
        file >> data; // read next value
        try
        {
          return lexical_cast<int>(data);
        }
        catch(boost::bad_lexical_cast &)
        {
          return 0;
        }
      }
    }

    return 0;
  }



  Double SmoothMax(Double left, Double right, Double beta)
  {
    assert(beta > 0 || beta == -1.0);
    if(beta == -1.0) return fmax(left, right);

    // the continuous Kreisselmeier and Steinhauser max approximation taken
    // from Sigmund; Morphology-based black and white filters for topology optimization; 2007
    // x = log ( sum(exp(beta * x_i)) / sum 1) / beta

    return log(0.5 * (exp(left * beta) + exp(right * beta))) / beta;
  }

  Double SmoothMax(const StdVector<Double>& values, Double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(values.GetSize() > 0);

    Double res = 0.0;

    if(beta == -1.0)
    {
      res = values[0];
      for(unsigned int i = 1, n = values.GetSize(); i < n; i++)
        res = fmax(res, values[i]);
    }
    else
    {
      // see  SmoothMax(double left, double right, double beta)
      // x = log ( sum(exp(beta * x_i)) / sum 1) / beta
      Double sum = 0.0;
      for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
        sum += exp(values[i] * beta);

      res = log(sum / (Double) values.GetSize()) / beta;
    }

    // LOG_DBG3(tools) << "SmoothMax v=" << values.ToString() << " beta=" << beta << " -> " << res;
    // std::cout << "SmoothMax v=" << values.ToString() << " beta=" << beta << " -> " << res << std::endl;
    return res;
  }

  Double SmoothMin(Double left, Double right, Double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(right > 0 && left > 0);

    if(beta == -1.0)
      return fmin(left, right);

    // see comment in CalcMaxApproximation()
    return 1.0 - log(0.5 * (exp((1.0 - left) * beta) + exp((1.0 - right) * beta))) / beta;
  }

  Double SmoothMin(const StdVector<Double>& values, Double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(values.GetSize() > 0);
    // see  SmoothMax(Double left, Double right, Double beta)

    if(beta == -1.0)
    {
      Double t = values[0];
      for(unsigned int i = 1, n = values.GetSize(); i < n; i++)
        t = fmin(t, values[i]);
      return t;
    }

    Double sum = 0.0;
    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
      sum += exp((1.0 - values[i]) * beta);

    return 1.0 - log(sum / (Double) values.GetSize()) / beta;
  }


  Double DerivSmoothMax(Double left, Double right, Double beta, int derive)
  {
    assert(derive == -1 || derive == 1); // left or right
    assert(beta > 0);

    Double exp_left  = exp(left * beta);
    Double exp_right = exp(right * beta);

    if(derive == -1)
      return exp_left / (exp_left + exp_right);
    else
      return exp_right / (exp_left + exp_right);
  }

  Double DerivSmoothMax(const StdVector<Double>& values, Double beta, unsigned int derive)
  {
    assert(beta > 0);

    Double my_exp = -1.0;
    Double sum = 0.0;

    if(values.GetSize() == 1)
      return 1.0;

    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
    {
      Double v = exp(values[i] * beta);
      sum += v;
      if(derive == i) my_exp = v; // not derive ist not a window based index
    }
    assert(my_exp != -1.0);

    return my_exp / sum;
  }

  Double DerivSmoothMin(Double left, Double right, Double beta, int derive)
  {
    assert(derive == -1 || derive == 1); // left or right
    Double exp_left  = exp((1.0 - left) * beta);
    Double exp_right = exp((1.0 - right) * beta);

    if(derive == -1)
      return exp_left / (exp_left + exp_right);
    else
      return exp_right / (exp_left + exp_right);
  }

  Double DerivSmoothMin(const StdVector<Double>& values, Double beta, unsigned int derive)
  {
    Double my_exp = -1.0;
    Double sum = 0.0;

    if(values.GetSize() == 1)
      return 1.0;

    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
    {
      Double v = exp((1.0 - values[i]) * beta);
      sum += v;
      if(derive == i) my_exp = v;
    }
    assert(my_exp != -1.0);

    return my_exp / sum;
  }

  Double SmoothAbs(Double x, Double eps)
  {
    assert(eps >= 0);
    return sqrt(x*x + eps*eps) - eps;
  }

  Double DerivSmoothAbs(Double x, Double eps)
  {
    assert(eps >= 0);
    assert(fabs(x) + eps > 0);
    return x / sqrt(x*x + eps*eps);
  }


  // explicit template instantiation
  template std::string ToString<Double>(const Double* data, unsigned int size);
  template std::string ToString<int>(const int* data, unsigned int size);
  template std::string ToString<unsigned int>(const unsigned int* data, unsigned int size);
  template std::string ToString<std::complex<Double> >(const std::complex<Double>* data, unsigned int size);

  template std::string ToString<Double>(const StdVector<Vector<Double> >& data, bool new_line);
  template std::string ToString<unsigned int>(const StdVector<StdVector<unsigned int> >& data, bool new_line);

}// namespace CoupledField
