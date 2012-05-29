#include "General/Exception.hh"
#include "BLASLAPACKInterface.hh"

namespace CoupledField {

  // ***********************************************************
  //   Conversion of C++ to F77 data type ( entryC -> entryF )
  // ***********************************************************
  void CC2F77( const double &v, F77real8 &val ) {
    val = (F77real8)v;
  }
  void CC2F77( const double &v, F77complex16 &val ) {
    EXCEPTION("Invalid conversion attempt double -> F77complex16");
  }
  void CC2F77( const std::complex<double> &v, F77real8 &val ) {
    EXCEPTION("Invalid conversion attempt complex<double> -> F77real8");
  }
  void CC2F77( const std::complex<double> &v, F77complex16 &val ) {
    val.real = (F77real8)v.real();
    val.imag = (F77real8)v.imag();
  }


  // ***********************************************************
  //   Conversion of F77 to C++ data type ( entryF -> entryC )
  // ***********************************************************
  void F772CC( const F77real8 &v, double &val ) {
    val = (double)v;
  }
  void F772CC( const F77real8 &v, std::complex<double> &val ) {
    EXCEPTION("Invalid conversion attempt F77real8 -> complex<double>");
  }
  void F772CC( const F77complex16 &v, double &val ) {
    EXCEPTION("Invalid conversion attempt F77complex16 -> double");
  }
  void F772CC( const F77complex16 &v, std::complex<double> &val ) {
    std::complex<double> aux(v.real,v.imag);
    val = aux;
  }

}
