// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "py-messenger.hh"

#include "Utils/StdVector.hh"
#include "General/environment.hh"
#include "pythonrun.h"
#include <boost/algorithm/string.hpp>

namespace CoupledField {
  

  // --- Declare static functions ---
  StdVector<std::string> PY_CFSMessenger::curParams_;

  // initialization of static structure
  static PyMethodDef CFSMethods[]  = {
    {"cfs",  PY_CFSMessenger::PY_CFSEval, METH_VARARGS,
     "Send messages to cfs."},
    {NULL, NULL, 0, NULL}    
  };

  PyMODINIT_FUNC initcfs() {
    (void) Py_InitModule("cfs", CFSMethods);
  }
  
  // end of static declaration block


  PY_CFSMessenger::PY_CFSMessenger() {
    
    // Initialize interpreter
    Py_Initialize();

    // Get global module and related dictionary
    PyObject* main_module =
      PyImport_AddModule("__main__");
    mainDict_ = PyModule_GetDict(main_module);
    
    // register cfs 'module' in python interpreter
    initcfs();

    // register (dummy) event procedures
    RegisterEvents();
    
  }
  
  
  PY_CFSMessenger::~PY_CFSMessenger() {

    // finalize script execution
    Py_Finalize();
  }
  
  void PY_CFSMessenger::ReadScriptFile ( const std::string & fileName ) {

    // open file
    FILE *scriptFile;
    scriptFile = fopen( fileName.c_str(), "r");

    // evaluate script file
    isEvaluating_ = true;
    PyObject* ret = NULL;
    ret = PyRun_File(scriptFile, fileName.c_str(), Py_file_input,
                     mainDict_, mainDict_);
    isEvaluating_ = false;

    // close file
    fclose( scriptFile );
    
    // check if error occured
    if( ret == NULL ) {
      std::cerr << "\n\n------------------------------------------------------\n";
      PyErr_Print();
      std::cerr << "\n------------------------------------------------------\n\n";
      PyErr_Print();
      std::string msg = "Error during reading of scriptfile! ";
      msg += "For details see previous error message.";
      EXCEPTION( msg.c_str() );
    }
  }
  
  bool PY_CFSMessenger::TriggerEvent( const EventType event, 
                                      const StdVector<std::string> & context) {
    
    // get name of event
    curEvent_ = eventNames_[event];

    std::string procName;
    procName += eventNames_[event] + "(";
    for ( UInt i=0; i<context.GetSize()-1; i++ ) {
      std::string actContext = context[i];
      boost::replace_all(actContext, "\n", "\\n");
      procName += "\"" + actContext + "\",  ";
    }
    if( context.GetSize() > 0 ) {
      std::string actContext = context.Last();
      boost::replace_all(actContext, "\n", "\\n");
      procName += "\"" + actContext + "\"";
    }
    procName += " )";

    isEvaluating_ = true;

    // evaluate script
    PyObject* ret = NULL;
    ret = PyRun_String( procName.c_str(), Py_eval_input,
                        mainDict_, mainDict_ );
    isEvaluating_ = false;
    
    if( ret == NULL ) {
      std::cerr << "\n\n------------------------------------------------------\n";
      PyErr_Print();
      std::cerr << "\n------------------------------------------------------\n\n";
      std::string error = "PYTHON error in function '";
      error+= eventNames_[event];
      error+= ". For details see previous error message.";
      //error += ExtractErrorInfo( type, pvalue, ptraceback );
      EXCEPTION( error.c_str() );
      return false;
    } else {
      return true;
    }
  }
  
  void PY_CFSMessenger::Warning( const char * msg, const char * const filename,
                                 const UInt numline) {
    
    std::stringstream warn;
    warn <<  "PYTHON warning in function '"  << curEvent_ << "':\n";
    for (UInt i = 0; i < curParams_.GetSize(); i++ ) {
      warn << curParams_[i] << " ";
    }
    warn << ": "<< msg;

    // After having generated the correct error string,
    // the bucket is passed back to the global error handler
    isEvaluating_ = false;
    ::Warning( warn.str().c_str(), filename, numline );
    isEvaluating_ = true;

  }
  
  void PY_CFSMessenger::Error( const char * msg, const char * const filename,
                               const UInt numline) {

    std::stringstream error;

    // Generate more accurate error message, if error occurs during
    // calling of a event procedure
    if (curEvent_ != std::string() ) {
      error <<  "PYTHON error in function '"  << curEvent_ << "':\n\n ";
      for (UInt i = 0; i < curParams_.GetSize(); i++ ) {
        error  << curParams_[i] << " ";
      }
      error << ": " << msg;
    } else {
      error << msg;
    }

    // After having generated the correct error string,
    // the bucket is passed back to the global error handler
    isEvaluating_ = false;
    EXCEPTION( error.str().c_str() );
  }
  
  
  PyObject* PY_CFSMessenger::PY_CFSEval( PyObject *self, PyObject *args) {

    bool success;

    // get first entry of args object
    PyObject * tuples;
    
    if( PyTuple_Size(args) == 1 ) {
      tuples = PyTuple_GetItem(args,0);
      
      // if object is not of type tuple, leave
      if( !PyTuple_Check(tuples) ) {
        EXCEPTION( "Argument to 'cfs' must be of type 'list'!");
        success = false;
      }
    } else {
      tuples = args;
    }

    // get number of entries in list
    UInt numArgs =  PyTuple_Size(tuples);
    

    //    std::cerr << "numArgs = " << numArgs << std::endl;
    // check if there are more than 1 entries present
    if (numArgs < 2 ) { 
      EXCEPTION( "cfs needs at least 2 arguments!" );
      success = false;
    } 
    
    // convert arguments into std::strings
    StdVector<std::string> stringArgs(numArgs+1);
    stringArgs[0] = "cfs";
    for( UInt i = 0; i < numArgs; i++) {
      const char * ch;
      if (!(ch = PyString_AsString(PyTuple_GetItem(tuples,i))) ) {
        return NULL;
      }
      //std::cerr << "successfully parsed string" << ch << std::endl;
      stringArgs[i+1] = std::string(ch);
    }

    // Save arguments to curent argument list
    curParams_ = stringArgs;

    // call language independent CFSEval function
    StdVector<std::string> retVal;
    success = CFSMessenger::CFSEval( stringArgs, retVal );

    if( success == true ) {
      PyObject * myList =  PyList_New( retVal.GetSize());
      for( UInt i = 0; i < retVal.GetSize(); i++ ) {
        PyList_SetItem( myList, i,PyString_FromString(retVal[i].c_str()));
      }
      return myList;
    } else {
      PyErr_SetString(	PyExc_Exception , strdup(errMsg_.c_str()));
      return NULL;
    }
  }
  
  void PY_CFSMessenger::RegisterEvents() {
   
    // First of all, generate mapping from event enumerations
    // to string representation
    eventNames_[CFS_Init] = "CFS_Init";
    eventNames_[CFS_PdeInit] = "CFS_PdeInit";
    eventNames_[CFS_ReadBCs] = "CFS_ReadBCs";
    eventNames_[CFS_AssembleMat] = "CFS_AssembleMat";
    eventNames_[CFS_AssembleRhs] = "CFS_AssembleRhs";
    eventNames_[CFS_SetBCs] = "CFS_SetBCs";
    eventNames_[CFS_CalcResults] = "CFS_CalcResults";
    eventNames_[CFS_Finish] = "CFS_Finish";
    
    // Check, if for all events the number
    // of parameters has been passed correctly
    if ( eventNumParams_.size() != eventNames_.size() ) {
      EXCEPTION("Size of eventNumParams_ and eventNames_ "
                << "mismatch!\n Please ensure, that for each "
                << "event the number of parameters is defined!");
    }
    
    // Register all possible events
    std::map<EventType,std::string>::iterator it;
    for ( it = eventNames_.begin(); it != eventNames_.end(); it++ ) {
      std::stringstream procName;
      procName << "def " << (*it).second;
      
      // Append dummy arguments for each parameter
      procName << " ( ";
      for ( UInt i = 0; i < eventNumParams_[(*it).first]-1; i++ ) {
        procName << "p" << i << " ,";
      }

      if( eventNumParams_[(*it).first] > 0 ) {
        procName << "p" << eventNumParams_[(*it).first];
      }
      
      procName << " ) : \n \t return\n";
      
       // register function ( = evaluate dummy proc body)
      PyRun_SimpleString( procName.str().c_str());
    }

    // In the end use the module
    PyRun_SimpleString( "from cfs import *");
  }
  
 //  std::string PY_CFSMessenger::ExtractErrorInfo( PyObject *exception,
//                                                  PyObject *value,
//                                                  PyObject *traceback ) {
//     PyErr_NormalizeException(&exception, &value, &traceback);
//     Py_INCREF(value);
//     PyObject *f = PyString_FromString("");
//     int err = 0;
//     std::string ret;

//   //   if (traceback && traceback != Py_None) {
// //       std::cerr << "Doing something nasty!";
// //       err = PyTraceBack_Print(traceback, f);
// //     }
// //     if( err != 0 ) {
// //       std::cerr << "->Getting the trace-back did not work!\n";
                                                                        
// //     }
//     std::cerr << "ret = " << ret << std::endl;
//     if ( PyObject_HasAttrString(value, "print_file_and_line")) {
//       std::cerr << "--> Has 'print_file_and_line'\n";
//       PyObject *message;
//       const char *filename, *text;
//       int lineno, offset;
//       if (!parse_syntax_error(value, &message, &filename,
//                               &lineno, &offset, &text))
//         PyErr_Clear();
//       else {
//         char buf[10];
//         ret += "  File \"";
//         if (filename == NULL)
//           ret += "<string>";
//         else
//           ret += filename;
//         ret += "\", line ";
//         ret += lineno;
//         ret += "\n";
//         //        if (text != NULL)
//         //  print_error_text(f, offset, text);
//         //Py_DECREF(value);
//         //value = message;
//         /* Can't be bothered to check all those
//            PyFile_WriteString() calls */
//         if (PyErr_Occurred())
//           err = -1;
//       }
//     }
//    if (err) {
//      /* Don't do anything else */
//    }
//    else if (PyClass_Check(exception)) {
//      PyClassObject* exc = (PyClassObject*)exception;
//      PyObject* className = exc->cl_name;
//      PyObject* moduleName =
//        PyDict_GetItemString(exc->cl_dict, "__module__");
     
//      if (moduleName == NULL)
//        err = PyFile_WriteString("<unknown>", f);
//      else {
//        char* modstr = PyString_AsString(moduleName);
//        if (modstr && strcmp(modstr, "exceptions"))
//          {
//            err = PyFile_WriteString(modstr, f);
//            err += PyFile_WriteString(".", f);
//          }
//      }
//      if (err == 0) {
//        if (className == NULL)
//          err = PyFile_WriteString("<unknown>", f);
//        else
//          err = PyFile_WriteObject(className, f,
//                                   Py_PRINT_RAW);
//      }
//    }
//    else
//      err = PyFile_WriteObject(exception, f, Py_PRINT_RAW);
//    if (err == 0) {
//      if (value != Py_None) {
//        PyObject *s = PyObject_Str(value);
//        /* only print colon if the str() of the
//           object is not the empty string
//        */
//        if (s == NULL)
//          err = -1;
//        else if (!PyString_Check(s) ||
//                 PyString_GET_SIZE(s) != 0)
//          err = PyFile_WriteString(": ", f);
//        if (err == 0)
//          err = PyFile_WriteObject(s, f, Py_PRINT_RAW);
//        Py_XDECREF(s);
//      }
//    }
//    if (err == 0)
//      err = PyFile_WriteString("\n", f);
   
//    Py_DECREF(value);
//    /* If an error happened here, don't show it.
//       XXX This is wrong, but too many callers rely on this behavior. */
//    if (err != 0)
//      PyErr_Clear();
   
   
   
   
//     //return std::string(PyString_AsString(f));
//     return ret;
   
   
   
//   }


//   int
//   PY_CFSMessenger::parse_syntax_error(PyObject *err, PyObject **message, const char **filename,
// 		   int *lineno, int *offset, const char **text)
// {
// 	long hold;
// 	PyObject *v;

// 	/* old style errors */
// 	if (PyTuple_Check(err))
// 		return PyArg_ParseTuple(err, "O(ziiz)", message, filename,
// 				        lineno, offset, text);

// 	/* new style errors.  `err' is an instance */

// 	if (! (v = PyObject_GetAttrString(err, "msg")))
// 		goto finally;
// 	*message = v;

// 	if (!(v = PyObject_GetAttrString(err, "filename")))
// 		goto finally;
// 	if (v == Py_None)
// 		*filename = NULL;
// 	else if (! (*filename = PyString_AsString(v)))
// 		goto finally;

// 	Py_DECREF(v);
// 	if (!(v = PyObject_GetAttrString(err, "lineno")))
// 		goto finally;
// 	hold = PyInt_AsLong(v);
// 	Py_DECREF(v);
// 	v = NULL;
// 	if (hold < 0 && PyErr_Occurred())
// 		goto finally;
// 	*lineno = (int)hold;

// 	if (!(v = PyObject_GetAttrString(err, "offset")))
// 		goto finally;
// 	if (v == Py_None) {
// 		*offset = -1;
// 		Py_DECREF(v);
// 		v = NULL;
// 	} else {
// 		hold = PyInt_AsLong(v);
// 		Py_DECREF(v);
// 		v = NULL;
// 		if (hold < 0 && PyErr_Occurred())
// 			goto finally;
// 		*offset = (int)hold;
// 	}

// 	if (!(v = PyObject_GetAttrString(err, "text")))
// 		goto finally;
// 	if (v == Py_None)
// 		*text = NULL;
// 	else if (! (*text = PyString_AsString(v)))
// 		goto finally;
// 	Py_DECREF(v);
// 	return 1;

// finally:
// 	Py_XDECREF(v);
// 	return 0;
// }
}
