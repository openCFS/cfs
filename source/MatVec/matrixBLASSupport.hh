#ifndef MATRIXBLASSUPPORT_HH_
#define MATRIXBLASSUPPORT_HH_

#include <def_use_blas.hh>

namespace CoupledField
{
#ifdef USE_BLAS

#define DGEMM dgemm_ //matrix matrix multiplication
#define ZGEMM zgemm_ //complex matrix matrix multiplication
#define DGEMV dgemv_ //matrix vector multiplication
#define ZGEMV zgemv_ //complex matrix vector multiplication


extern "C"{
  void DGEMM(char*, char*, int*, int*, int*, double*, double*, int* , double*, int*, double*, double*, int*);
	void ZGEMM(char*, char*, int*, int*, int*, std::complex<double>*, std::complex<double>*,int*, std::complex<double>*,int*,std::complex<double>*, std::complex<double>*,int*);
	void DGEMV(char*, int*, int*, double*, double*, int*, const double*, int*, double*, double*, int*);
	void ZGEMV(char*, int*, int*, std::complex<double>*, std::complex<double>*,int*, const std::complex<double>*, int*, std::complex<double>*, std::complex<double>*,int*);

}

}
#endif
#endif /* MATRIXBLASSUPPORT_HH_ */
