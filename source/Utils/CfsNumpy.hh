#ifndef UTILS_CFSNUMPY_HH
#define UTILS_CFSNUMPY_HH

// Single numpy C-API entry point for all of cfs. By default each compiled file (translation unit) including <numpy/arrayobject.h>
// gets its own static, NULL PyArray_API table that only import_array() initializes - so table
// functions (PyArray_Check, array creation) crash from any TU that never ran import_array.
// PY_ARRAY_UNIQUE_SYMBOL turns that per-TU table into ONE shared symbol; every TU includes this
// header, so the single table filled by CfsImportNumpy() (in the owner TU that defines
// CFS_NUMPY_OWNER) is visible everywhere. Owner: Utils/PythonKernelFunctions.cc.
#define PY_ARRAY_UNIQUE_SYMBOL CFS_PyArray_API
#ifndef CFS_NUMPY_OWNER
  #define NO_IMPORT_ARRAY
#endif
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

namespace CoupledField
{
  /** Initialize the shared numpy C-API table. Call once after Py_Initialize() and before any numpy use. */
  void CfsImportNumpy();
}

#endif
