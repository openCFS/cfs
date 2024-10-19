#define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
//#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <Python.h>

#include "Vector.hh"

namespace CoupledField
{

template<typename T>
Vector<T>::Vector(pyObject* obj, bool decref)
{
  this->capacity_ = 0;
  this->size_ = 0;
  this->memBelongsToMe_= true;
  this->data_= NULL;

  Fill(obj, decref);
}

template<typename T>
void Vector<T>::Fill(pyObject* obj, bool decref)
{
  assert(obj != NULL);

  PyArrayObject* pao = (PyArrayObject*) obj;

  if(PyArray_NDIM(pao) != 1)
    EXCEPTION("numpy array is not a vector but has dimension " + std::to_string(PyArray_NDIM(pao)));

  Resize(PyArray_DIM(pao,0));

  for(unsigned int i = 0; i < size_; i++)
    data_[i] = *((T*) PyArray_GETPTR1(pao, i));

  if(decref)
    Py_DECREF(obj);
}

template<typename T>
void Vector<T>::Export(pyObject* obj)
{
  assert(obj != NULL);

  PyArrayObject* pao = (PyArrayObject*) obj;

  if(PyArray_NDIM(pao) != 1)
    EXCEPTION("numpy array is not a vector but has dimension " + std::to_string(PyArray_NDIM(pao)));

  if(PyArray_DIM(pao,0) != size_)
    EXCEPTION("numpy array has size " + std::to_string(PyArray_DIM(pao,0)) + " but cfs vector has " + std::to_string(size_));

  for(unsigned int i = 0; i < size_; i++)
    *((T*) PyArray_GETPTR1(pao,i)) = data_[i];
}

template class Vector<Double>;
template class Vector<Complex>;
template class Vector<Integer>;
template class Vector<UInt>;

} // end of name space
