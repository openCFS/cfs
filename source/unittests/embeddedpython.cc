#define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
//#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <Python.h>


#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <def_use_embedded_python.hh>

#ifdef USE_EMBEDDED_PYTHON
  #include "MatVec/Vector.hh"
#endif

using namespace CoupledField;


static int callback_val=5;

/* expects a long of value 4 and returns callback_val which is 5 */
static pyObject* cfs_val(pyObject *self, pyObject *args)
{
  long c;
  PyArg_ParseTuple(args, "l", &c);
  assert(c == 4);

  return PyLong_FromLong(callback_val);
}


/** get a numpy array (dim = 1) and a long with the size of the array */
static pyObject* cfs_vec(pyObject *self, pyObject *args)
{
  PyArrayObject *vec1;
  PyArrayObject *vec2;
  long test;
  // in python 'verify' is a bool. Here it must be of type 'int'
  // e.g. see https://github.com/numpy/numpy/issues/12400
  int verify;
  std::cout << "cfs_vec: args= " << PyTuple_Size(args) << std::endl;
  PyArg_ParseTuple(args, "O!O!np", &PyArray_Type, &vec1, &PyArray_Type, &vec2, &test, &verify);
  BOOST_TEST(PyArray_NDIM(vec1) == 1);
  BOOST_TEST(PyArray_TYPE(vec1) == NPY_DOUBLE);
  BOOST_TEST(verify);
  int n = PyArray_DIM(vec1,0);
  BOOST_TEST(n == test);
  //  BOOST_TEST(PyArray_DIM((PyArrayObject*) obj,0) == n);

  return PyLong_FromLong(n);
}

static pyObject* cfs_mod_mat(pyObject *self, pyObject *args){
  // consider list of 2d arrays as and 3d array
  PyArrayObject *list;
  PyArg_ParseTuple(args, "O!", &PyArray_Type, &list);

  if(!list)
    PyErr_Print();

  // number of matrices inside list
  int n_mat = PyArray_DIM(list,0);
  // assume all matrices have same dimensions
  int rows = PyArray_DIM(list,1);
  int cols = PyArray_DIM(list,2);

  std::cout << "Got " << n_mat << " matrices with dim=(" << rows << "," << cols << ") from python" << std::endl;

  // loop over matrices
  for (int e = 0; e < n_mat; e++) {
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++){
        (*((double*) PyArray_GETPTR3(list,e,i,j))) = 11 + e;
      }
  }
  BOOST_TEST((*((double*) PyArray_GETPTR3(list,0,0,0))) == 11);
  BOOST_TEST((*((double*) PyArray_GETPTR3(list,1,0,0))) == 12);
  BOOST_TEST((*((double*) PyArray_GETPTR3(list,2,0,0))) == 13);
  return PyLong_FromLong(1);
}

static PyMethodDef cfs_methods[] = {
    {"val", cfs_val, METH_VARARGS, "Shall return 5"},
    {"vec", cfs_vec, METH_VARARGS, "get array and size, return size+1"},
    {"mod_mat", cfs_mod_mat, METH_VARARGS, "modify a numpy array of 2d matrices (numpy.array)"},
    {NULL, NULL, 0, NULL}
};

static PyModuleDef cfs_modules = {
    PyModuleDef_HEAD_INIT, "cfs", NULL, -1, cfs_methods, NULL, NULL, NULL, NULL
};

static pyObject* PyInit_cfs(void)
{
  // https://stackoverflow.com/questions/37943699/crash-when-calling-pyarg-parsetuple-on-a-numpy-array
  import_array();

  return PyModule_Create(&cfs_modules);
}


BOOST_AUTO_TEST_CASE(embedded_python)
{
  // needs to be done before Py_Initialize
  PyImport_AppendInittab("cfs", &PyInit_cfs);

  Py_Initialize();

  PyRun_SimpleString("import os; print('embedded python: runs in ' + os.getcwd())");

  // test system
  pyObject* version = PySys_GetObject("version");
  std::cout << "pyObject version=" << version << std::endl;
  if(!version)
    PyErr_Print();

  const char* c_str = PyUnicode_AsUTF8(version);
  std::cout << "version c_str=" << c_str << std::endl;

  Py_XDECREF(version);

  // assume we are in the build directory
  boost::filesystem::path test = boost::filesystem::path("../source/unittests/embeddedpython.py");
  std::cout << "test filename=" << test.filename() << std::endl;
  std::cout << "test path=" << test.parent_path()  << std::endl; // is "" in case of test = "embeddedpython.py"
  boost::filesystem::path path = boost::filesystem::absolute(test.parent_path()); // is pwd for ""
  std::cout << "test absolute path=" << path  << std::endl;
  std::cout << "test no extension=" << boost::filesystem::change_extension(test.filename(), "") << std::endl;



  // add it to the system path
  if(boost::filesystem::is_directory(path))
  {
    pyObject* sysPath = PySys_GetObject((char*) "path"); // must not decref after append
    PyList_Append(sysPath, PyUnicode_FromString(path.string().c_str()));
  }
  else
    std::cout << "WARNING: not running in build directory, make sure */source/unittest is in PYTHONPATH" << std::endl;

  // https://docs.python.org/3.8/extending/embedding.html
  pyObject* pModule = PyImport_ImportModule("embeddedpython");

  if(pModule != NULL)
  {
    // call python function
    pyObject* pFunc = pyObject_GetAttrString(pModule, "inc_a");
    assert(pFunc && PyCallable_Check(pFunc));
    pyObject* pValue = pyObject_CallObject(pFunc, NULL);
    BOOST_TEST(pValue);
    BOOST_TEST(PyLong_AsLong(pValue) == 1);
    std::cout << "inc_c returned: " << PyLong_AsLong(pValue) << std::endl;
    Py_XDECREF(pValue);
    Py_XDECREF(pFunc); // could be NULL

    // get numpy object from python
    // https://stackoverflow.com/questions/43437885/how-do-i-access-a-numpy-array-in-embedded-python-from-c
    PyArrayObject* array = (PyArrayObject*) pyObject_GetAttrString(pModule, "A");
    BOOST_TEST(array);
    std::cout << "dims: " << PyArray_NDIM(array) << ":" << PyArray_DIM(array,0) << "," << PyArray_DIM(array,1) << std::endl;
    BOOST_TEST(PyArray_NDIM(array) == 2);
    for(int i = 0; i < PyArray_DIM(array,0); i++)
      for(int j = 0; j < PyArray_DIM(array,1); j++) {
        (*((double*) PyArray_GETPTR2(array, i,j)))++;
         std::cout << i << ":" << j << " -> " << *((double*) PyArray_GETPTR2(array, i,j)) << std::endl;
      }
    Py_XDECREF(array);

    pyObject* pA = pyObject_GetAttrString(pModule, "print_A");
    assert(pA && PyCallable_Check(pA));
    pyObject_CallObject(pA, NULL);
    Py_DECREF(pA);

    // call python method with more than one return value
    pyObject* mv = pyObject_GetAttrString(pModule, "many_values");
    assert(mv && PyCallable_Check(mv));
    pyObject* arg = PyTuple_New(2);
    PyTuple_SetItem(arg, 0, PyLong_FromLong(6));
    PyTuple_SetItem(arg, 1, PyFloat_FromDouble(3.14));
    pyObject* ret = pyObject_CallObject(mv, arg);
    Py_XDECREF(arg); // I guess val1 and val1 are decremented here?!
    if(!ret)
      PyErr_Print();
    Py_XDECREF(arg);
    Py_XDECREF(mv);
    Py_XDECREF(ret);
    BOOST_TEST_REQUIRE(ret);
    BOOST_TEST_REQUIRE(PyTuple_Size(ret) == 2);

    // test Vector and numpy conversion
    pyObject* V = pyObject_GetAttrString(pModule, "V");
    if(!V)
      PyErr_Print();
    Vector<double> v(V,true);
    std::cout << v.ToString() << std::endl;

    // test for numpy type
    pyObject* a = (pyObject*) pyObject_GetAttrString(pModule, "a");

    if(!a)
      PyErr_Print();

    std::cout << "a as long is " << PyLong_AsLong(a) << std::endl;
    // std::cout << "isinstance : " << pyObject_IsInstance(a,(pyObject*) &PyArray_Type) << std::endl;
    // std::cout << "check a: " << PyArray_Check(a) << std::endl;
    Py_XDECREF(a);

    // test parsing boolean
    pyObject* derivative = pyObject_GetAttrString(pModule, "derivative");
    BOOST_TEST(!pyObject_IsTrue(derivative));

    // test dictionary
    pyObject* dict = PyDict_New();
    PyDict_SetItemString(dict, "franz", PyUnicode_FromString("franzl"));
    PyDict_SetItemString(dict, "hans", PyLong_FromLong(5));
    arg = PyTuple_New(1);
    PyTuple_SetItem(arg, 0, dict);
    pFunc = pyObject_GetAttrString(pModule, "print_dict");
    assert(pFunc && PyCallable_Check(pFunc));
    pValue = pyObject_CallObject(pFunc, arg);
    Py_XDECREF(arg);
    Py_XDECREF(pValue);
    Py_XDECREF(pFunc);


    // call python function with a tuple (not numpy), and empty tuple and a number
    pyObject* gv = pyObject_GetAttrString(pModule, "get_vec");
    assert(gv && PyCallable_Check(gv));
    pyObject* gva = PyTuple_New(3); // gvl + gve + 2
    pyObject* gvl = PyTuple_New(2); // the list object as first attribute
    PyTuple_SetItem(gvl, 0, PyFloat_FromDouble(0.1));
    PyTuple_SetItem(gvl, 1, PyFloat_FromDouble(0.2));

    pyObject* gve = PyTuple_New(0); // empty tuple

    PyTuple_SetItem(gva, 0, gvl);
    PyTuple_SetItem(gva, 1, gve);
    PyTuple_SetItem(gva, 2, PyLong_FromLong(2)); // gvl size

    pyObject* gvr = pyObject_CallObject(gv, gva);

    if(!gvr) PyErr_Print();
    Py_XDECREF(gvr); // I guess val1 and val1 are decremented here?!
    Py_XDECREF(gve);
    Py_XDECREF(gva); // do NOT Py_XDECREF(gvl); -> segfault when looping
    Py_XDECREF(gv);



    // PyArray_Zeros, PyArray_SimpleNew(), PyArray_SimpleNewFromData and PyArray_Resize don't work.  :(
    // pyObject* resized = PyArray_Resize(vec, &ad, 0, NPY_CORDER);
    //pyObject* new_vec = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data.GetPointer());


    Py_DECREF(pModule); // we know it is not NULL
  }
  else
  {
    PyErr_Print();
    std::cout << "make sure the current directory is in the PYTHONPATH!";
    assert(false);
  }

  Py_Finalize();
}
