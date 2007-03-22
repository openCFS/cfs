// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>

#include "utils/defs.hh"
#include "utils/environment.hh"
#include "utils/Exception.hh"

namespace OLAS {


  Exception::Exception(const char* msg, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      PostInit();
  }  

  Exception::Exception(std::string msg, const int line, const char* function, Segfault seg)
  {
      Init(NULL, line, function, seg);
      this->string_1 = msg;
      PostInit();      
  }  
  
  Exception::Exception(const char* msg, Double value, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->double_1 = value;
      PostInit();      
  }  

  Exception::Exception(const char* msg, UInt value, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->uint_1 = value;
      PostInit();      
  }  

  Exception::Exception(const char* msg, UInt value1, UInt value2, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->uint_1 = value1;
      this->uint_2 = value2;
      PostInit();      
  }  


  Exception::Exception(const char* msg, Integer value, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->int_1 = value;
      PostInit();      
  }  

  Exception::Exception(const char* msg, Integer value1, Integer value2, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->int_1 = value1;
      this->int_2 = value2;
      PostInit();      
  }  

  Exception::Exception(const char* msg, const std::string& value1, Integer value2, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->string_1 = value1;
      this->int_1    = value2;
      PostInit();      
  }  



  Exception::Exception(const char* msg, const std::string& value, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->string_1 = value;
      PostInit();      
  }  

  Exception::Exception(const char* msg, const std::string& value1,  const std::string& value2, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->string_1 = value1;
      this->string_2 = value2;
      PostInit();      
  }  

  Exception::Exception(const char* msg, const std::string& value1, const std::string& value2, const std::string& value3, const int line, const char* function, Segfault seg)
  {
      Init(msg, line, function, seg);
      this->string_1 = value1;
      this->string_2 = value2;
      this->string_3 = value3;
      PostInit();      
  }  


  Exception::~Exception()
  {
  }

  void Exception::Init(const char* msg, const int line, const char* function, Segfault seg)
  {
      this->msg       = msg != NULL ? msg : "";
      this->line      = line;
      this->function  = function != NULL ? function : "";
      this->segfault  = seg;

      // default values for "not set"
      this->double_1 = 12345.7;
      this->uint_1   = 1234567;
      this->uint_2   = 1234567;
      this->int_1    = 1234567;
      this->int_2    = 1234567;
      this->string_1 = "not_set";
      this->string_2 = "not_set";
      this->string_3 = "not_set";
  }
  
  void Exception::PostInit()
  {
      if(segfault == NO_SEGFAULT) return;

      // killme -todo, use command-line parameter and determin if we
      // are compilled with debug information
      bool debug = true;
      
      // only AUTO and SEGFAULT here
      if(segfault == AUTO && !debug) return;
      
      Dump();
      
      // force segfault to gain stacktrace
      int* ip = NULL;
      (*ip)++;
  }

  void Exception::ToStringCore(std::string& out)
  {
     std::stringstream ss;
     ss << "msg='" << msg << "'";
     if(double_1 != 12345.7)     ss << " value=" << double_1;
     if(uint_1   != 1234567)     ss << " value=" << uint_1;
     if(uint_2   != 1234567)     ss << " value=" << uint_2;
     if(int_1    != 1234567)     ss << " value=" << int_1;
     if(int_2    != 1234567)     ss << " value=" << int_2;
     if(string_1 != "not_set")   ss << " value='" << string_1 << "'";
     if(string_2 != "not_set")   ss << " value='" << string_2 << "'";
     if(string_3 != "not_set")   ss << " value='" << string_3 << "'";
     ss << std::endl;
     if(line != -1)              ss << " line=" << line;
     if(function.length() != 0)  ss << " function=" << function;
     out.assign(ss.str());
  }  

  void Exception::ToString(std::string& out)
  {
     ToStringCore(out); // temporary
     std::stringstream ss;
     ss << "Exception: " << out;
     out.assign(ss.str());
  }  

  std::string Exception::ToString()
  {
      std::string string;
      ToString(string);
      return string;
  }
        
  void Exception::Dump()
  {
     std::string out;
     ToString(out);
     std::cerr << out << std::endl; 
 }     
        

} // namespace
