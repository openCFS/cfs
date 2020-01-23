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
    return (std::abs(in.imag()) > 1e-16) ?                   
        std::atan2(in.imag(),in.real() )*180/M_PI : 
        ( in.real() < 0.0 ) ? 180 : 0 ; 
  }


  //! Convert (ampl,phase) => (real,imag)
  Complex AmplPhaseToComplex( Double val, Double phase ) {
    return Complex( val * std::cos( phase / 180 * M_PI ),
                    val * std::sin( phase / 180 * M_PI ) ); 
  }

  //! Convert (ampl,phase) => real
  Double AmplPhaseToReal( Double val, Double phase ) {
    return val * std::cos( phase / 180 * M_PI );
  }

  //! Convert (ampl,phase) => imag
  Double AmplPhaseToImag( Double val, Double phase ) {
    return val * std::sin( phase / 180 * M_PI );
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

    return std::sqrt(result);
  }

  Double NormL2(const Double* data1, const Double* data2, const UInt size)
  {
    Double result = 0.0;
    for(UInt i = 0; i < size; i++)
      result += (data1[i] - data2[i]) * (data1[i] - data2[i]);

    return std::sqrt(result);
  }

  double NormL2(const SingleVector* data, const SingleVector* data2)
  {
    assert(data != NULL && data2 != NULL);

    if(data->GetSize() != data2->GetSize()) EXCEPTION("incompatible sizes");
    if(data->GetEntryType() != data2->GetEntryType()) EXCEPTION("incompatible entry types");

    if(data->GetEntryType() == BaseMatrix::COMPLEX)
      return dynamic_cast<const Vector<Complex>& >(*data).NormL2(dynamic_cast<const Vector<Complex>& >(*data2));
    else
      return dynamic_cast<const Vector<double>& >(*data).NormL2(dynamic_cast<const Vector<double>& >(*data2));
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
  std::string ToString<std::complex<double> >(const std::complex<double>* data, unsigned int size)
  {
    std::ostringstream os;

    // prepare for copy and paste to matlab
    os << "[";
    for(unsigned int i = 0; i < size; i += 1)
      os << data[i].real() << " + " << data[i].imag() << "i" << ((i < size-1) ? ", ": "]");

    return os.str();
  }


  double Average(const double* data, unsigned int size)
  {
    double sum = 0.0;
    for(unsigned int i = 0; i < size; i++)
      sum += *(data + i);

    return sum / size;
  }

  double StandardDeviation(const double* data, unsigned int size)
  {
    // expected value
    double e = Average(data, size);

    // variance
    double v = 0.0;
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



  double SmoothMax(double left, double right, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    if(beta == -1.0) return std::max(left, right);

    // the continuous Kreisselmeier and Steinhauser max approximation taken
    // from Sigmund; Morphology-based black and white filters for topology optimization; 2007
    // x = log ( sum(exp(beta * x_i)) / sum 1) / beta

    return std::log(0.5 * (std::exp(left * beta) + std::exp(right * beta))) / beta;
  }

  double SmoothMax(const StdVector<double>& values, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(values.GetSize() > 0);

    double res = 0.0;

    if(beta == -1.0)
    {
      res = values[0];
      for(unsigned int i = 1, n = values.GetSize(); i < n; i++)
        res = std::max(res, values[i]);
    }
    else
    {
      // see  SmoothMax(double left, double right, double beta)
      // x = log ( sum(exp(beta * x_i)) / sum 1) / beta
      double sum = 0.0;
      for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
        sum += std::exp(values[i] * beta);

      res = std::log(sum / (double) values.GetSize()) / beta;
    }

    // LOG_DBG3(tools) << "SmoothMax v=" << values.ToString() << " beta=" << beta << " -> " << res;
    // std::cout << "SmoothMax v=" << values.ToString() << " beta=" << beta << " -> " << res << std::endl;
    return res;
  }

  double SmoothMin(double left, double right, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(right > 0 && left > 0);

    if(beta == -1.0)
      return std::min(left, right);

    // see comment in CalcMaxApproximation()
    return 1.0 - std::log(0.5 * (std::exp((1.0 - left) * beta) + std::exp((1.0 - right) * beta))) / beta;
  }

  double SmoothMin(const StdVector<double>& values, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(values.GetSize() > 0);
    // see  SmoothMax(double left, double right, double beta)

    if(beta == -1.0)
    {
      double t = values[0];
      for(unsigned int i = 1, n = values.GetSize(); i < n; i++)
        t = std::min(t, values[i]);
      return t;
    }

    double sum = 0.0;
    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
      sum += std::exp((1.0 - values[i]) * beta);

    return 1.0 - std::log(sum / (double) values.GetSize()) / beta;
  }


  double DerivSmoothMax(double left, double right, double beta, int derive)
  {
    assert(derive == -1 || derive == 1); // left or right
    assert(beta > 0);

    double exp_left  = std::exp(left * beta);
    double exp_right = std::exp(right * beta);

    if(derive == -1)
      return exp_left / (exp_left + exp_right);
    else
      return exp_right / (exp_left + exp_right);
  }

  double DerivSmoothMax(const StdVector<double>& values, double beta, unsigned int derive)
  {
    assert(beta > 0);

    double my_exp = -1.0;
    double sum = 0.0;

    if(values.GetSize() == 1)
      return 1.0;

    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
    {
      double v = std::exp(values[i] * beta);
      sum += v;
      if(derive == i) my_exp = v; // not derive ist not a window based index
    }
    assert(my_exp != -1.0);

    return my_exp / sum;
  }

  double DerivSmoothMin(double left, double right, double beta, int derive)
  {
    assert(derive == -1 || derive == 1); // left or right
    double exp_left  = std::exp((1.0 - left) * beta);
    double exp_right = std::exp((1.0 - right) * beta);

    if(derive == -1)
      return exp_left / (exp_left + exp_right);
    else
      return exp_right / (exp_left + exp_right);
  }

  double DerivSmoothMin(const StdVector<double>& values, double beta, unsigned int derive)
  {
    double my_exp = -1.0;
    double sum = 0.0;

    if(values.GetSize() == 1)
      return 1.0;

    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
    {
      double v = std::exp((1.0 - values[i]) * beta);
      sum += v;
      if(derive == i) my_exp = v;
    }
    assert(my_exp != -1.0);

    return my_exp / sum;
  }

  double SmoothAbs(double x, double eps)
  {
    assert(eps >= 0);
    return std::sqrt(x*x + eps*eps) - eps;
  }

  double DerivSmoothAbs(double x, double eps)
  {
    assert(eps >= 0);
    assert(abs(x) + eps > 0);
    return x / std::sqrt(x*x + eps*eps);
  }


  double MathParse(const std::string& expr)
  {
    // obtain handle
    MathParser* parser = domain->GetMathParser();
    MathParser::HandleType handle = parser->GetNewHandle(false);

    // Set expression and evaluate
    parser->SetExpr(handle, expr);
    double ret = parser->Eval(handle);

    // release handle
    parser->ReleaseHandle(handle);

    return ret;
  }

StdVector<Vector<double> > GetNewtonCotes()
{
  StdVector<Vector<double> > nc(11);
  nc[0].Resize(0);

  nc[1].Resize(2);
  nc[1][0] = 1.0/2;
  nc[1][1] = 1.0/2;

  nc[2].Resize(3);
  nc[2][0] = 1.0/6;
  nc[2][1] = 2.0/3;
  nc[2][2] = 1.0/6;

  nc[3].Resize(4);
  nc[3][0] = 1.0/8;
  nc[3][1] = 3.0/8;
  nc[3][2] = 3.0/8;
  nc[3][3] = 1.0/8;

  nc[4].Resize(5);
  nc[4][0] = 7.0/90;
  nc[4][1] = 16.0/45;
  nc[4][2] = 2.0/15;
  nc[4][3] = 16.0/45;
  nc[4][4] = 7.0/90;

  nc[5].Resize(6);
  nc[5][0] = 19.0/288;
  nc[5][1] = 25.0/96;
  nc[5][2] = 25.0/144;
  nc[5][3] = 25.0/144;
  nc[5][4] = 25.0/96;
  nc[5][5] = 19.0/288;

  nc[6].Resize(7);
  nc[6][0] = 41.0/840;
  nc[6][1] = 9.0/35;
  nc[6][2] = 9.0/280;
  nc[6][3] = 34.0/105;
  nc[6][4] = 9.0/280;
  nc[6][5] = 9.0/35;
  nc[6][6] = 41.0/840;

  nc[7].Resize(8);
  nc[7][0] = 751.0/17280;
  nc[7][1] = 3577.0/17280;
  nc[7][2] = 49.0/640;
  nc[7][3] = 2989.0/17280;
  nc[7][4] = 2989.0/17280;
  nc[7][5] = 49.0/640;
  nc[7][6] = 3577.0/17280;
  nc[7][7] = 751.0/17280;

  nc[8].Resize(9);
  nc[8][0] = 989.0/28350;
  nc[8][1] = 2944.0/14175;
  nc[8][2] = -464.0/14175;
  nc[8][3] = 5248.0/14175;
  nc[8][4] = -454.0/2835;
  nc[8][5] = 5248.0/14175;
  nc[8][6] = -464.0/14175;
  nc[8][7] = 2944.0/14175;
  nc[8][8] = 989.0/28350;

  nc[9].Resize(10);
  nc[9][0] = 2857.0/89600;
  nc[9][1] = 15741.0/89600;
  nc[9][2] = 27.0/2240;
  nc[9][3] = 1209.0/5600;
  nc[9][4] = 2889.0/44800;
  nc[9][5] = 2889.0/44800;
  nc[9][6] = 1209.0/5600;
  nc[9][7] = 27.0/2240;
  nc[9][8] = 15741.0/89600;
  nc[9][9] = 2857.0/89600;

  nc[10].Resize(11);
  nc[10][0] = 16067.0/598752;
  nc[10][1] = 26575.0/149688;
  nc[10][2] = -16175.0/199584;
  nc[10][3] = 5675.0/12474;
  nc[10][4] = -4825.0/11088;
  nc[10][5] = 17807.0/24948;
  nc[10][6] = -4825.0/11088;
  nc[10][7] = 5675.0/12474;
  nc[10][8] = -16175.0/199584;
  nc[10][9] = 26575.0/149688;
  nc[10][10] = 16067.0/598752;

  return nc;
}

Vector<double> Linspace(double start, double end, int n)
{
  Vector<double> res(n);
  for(int i = 0; i < n; n++)
    res[i] = start + i*(end-start)/(n-1);
  return res;
}

void Sub2Ind(Vector<unsigned int> size, StdVector<int> sub, unsigned int &ind)
{
  assert(size.GetSize() >= sub.GetSize());

  // cumulative product
  StdVector<int> cumprod = StdVector<int>(size.GetSize());
  cumprod[0] = size[0];
  for(unsigned int i = 1; i < size.GetSize(); ++i) {
    cumprod[i] = cumprod[i-1] * size[i];
  }

  int idx = sub[0] + 1; // zero based
  for(unsigned int i = 1; i < sub.GetSize(); ++i) {
    idx += sub[i] * cumprod[i-1];
  }
  ind = idx - 1; // zero based
}

void Ind2Sub(Vector<unsigned int> size, unsigned int ind, StdVector<int> &sub)
{
  sub.Resize(size.GetSize());

  // cumulative product
  StdVector<int> cumprod = StdVector<int>(size.GetSize());
  cumprod[0] = size[0];
  for(unsigned int i = 1; i < size.GetSize(); ++i) {
    cumprod[i] = cumprod[i-1] * size[i];
  }

  unsigned int idx = ind + 1; //zero based
  for(unsigned int i = sub.GetSize()-1; i > 0; i--) {
    int v = (idx-1) % cumprod[i-1] + 1;
    sub[i] = (idx - v) / cumprod[i-1];
    idx = v;
  }
  sub[0] = idx-1;
}

Vector<double> LogspaceBase(double start, double end, int n)
{
  Vector<double> res(n);
  for(int i = 0; i < n; i++)
    res[i] = std::pow(10, start + i*(end-start)/(n-1)); // pow(10,linspace(s,e,n))
  return res;
}

  // explicit template instantiation
  template std::string ToString<double>(const double* data, unsigned int size);
  template std::string ToString<int>(const int* data, unsigned int size);
  template std::string ToString<unsigned int>(const unsigned int* data, unsigned int size);
  // std::string ToString<std::complex<double> >(const std::complex<double>* data, unsigned int size) is expicitly instantiated

  template std::string ToString<double>(const StdVector<Vector<double> >& data, bool new_line);
  template std::string ToString<unsigned int>(const StdVector<StdVector<unsigned int> >& data, bool new_line);

}// namespace CoupledField
