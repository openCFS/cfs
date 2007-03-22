// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef EXCEPTIONS_HH_
#define EXCEPTIONS_HH_

#include <iostream>
#include <string>



namespace OLAS
{
 
  class Exception
  {
     public:
        /** The state of forced segfault to give a stack trace */
        typedef enum { NO_SEGFAULT, SEGFAULT, AUTO } Segfault;
     
        Exception(const char* msg, double value, const int line = -1, const char* function = NULL, Segfault = AUTO);     

        Exception(const char* msg, unsigned int value, const int line = -1, const char* function = NULL, Segfault = AUTO);        
     
        Exception(const char* msg, unsigned int value1, unsigned int value1, const int line = -1, const char* function = NULL, Segfault = AUTO);     
     
        Exception(const char* msg, int value, const int line = -1, const char* function = NULL, Segfault = AUTO);
        
        Exception(const char* msg, int value1, int value2, const int line = -1, const char* function = NULL, Segfault = AUTO);

        Exception(const char* msg, const std::string& value1, int value2, const int line = -1, const char* function = NULL, Segfault = AUTO);

        Exception(const char* msg, const std::string& value, const int line = -1, const char* function = NULL, Segfault = AUTO);     
     
        Exception(const char* msg, const std::string& value1, const std::string& value2,  const int line = -1, const char* function = NULL, Segfault = AUTO);
     
        Exception(const char* msg, const std::string& value1, const std::string& value2, const std::string& value3, const int line = -1, const char* function = NULL, Segfault = AUTO);     
     
        Exception(const char* msg, const int line = -1, const char* function = NULL, Segfault = AUTO);
        
        Exception(std::string msg, const int line = -1, const char* function = NULL, Segfault = AUTO);        
        
        ~Exception();
        
        std::string ToString();
        
        void ToString(std::string& out);

        /** makes a ToString() to std::cerr */
        void Dump();
       
        
     protected:
        void ToStringCore(std::string& out);

        double       double_1;
        unsigned int uint_1;   
        unsigned int uint_2;
        int          int_1;
        int          int_2;  
        std::string  string_1;
        std::string  string_2;
        std::string  string_3;
        std::string  msg;
        std::string  file;
        std::string  function;
        int line;
     private:
        /** the common constructor. call this first in the constructors as it sets values to defaults
         * which might be overwritten my parameters. */
        void Init(const char* msg, const int line, const char* function, Segfault segfault);
        
        /** The post-init to be called as last action in the constructors. Evaluates segfault */
        void PostInit();

        /** determines if a segfault should be forced to create a stack trace */ 
        Segfault segfault;           
  }; 
} // namespace



#endif /*EXCEPTIONS_HH_*/
