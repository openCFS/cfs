// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>

#include "exception.hh"

namespace CoupledField {

    bool Exception::segfault_ = false;

    /** this is the constructor for the EXCEPTION macro */
    Exception::Exception( const Exception* reason,
               const char* const fileName, 
               const unsigned int lineNum,
               const char* const message) throw ()
    {
        init(reason, fileName, lineNum, message);
    }

    /** this is the constructor to be called manually */
    Exception::Exception( const std::string& message,
               const char* const fileName, 
               const unsigned int lineNum) throw ()
    {
        init(NULL, fileName, lineNum, message.c_str());
    }

    void Exception::init(const Exception* reason,
               const char* const fileName, 
               const unsigned int lineNum,
               const char* const message) throw () 
    {
        this->fileName_ = fileName;
        this->lineNum_  = lineNum;
        this->msg_      = message;
         
        if(reason)
            reason_ = new Exception(*reason);
        else
            reason_ = NULL;

        std::ostringstream ostream;
        
        ostream << "\n" << msg_;

        if( fileName_ != "" && fileName_ != "NO_FILENAME") {
            ostream << "\nIn file '" << fileName_ << "'";
        }
        
        if( lineNum_ != 0 ) {
            ostream << " at line " << lineNum_ << std::endl;
        }
        
        ostream << "\n";
        
        if(reason_)
        {
            ostream << reason_->what();
        }

        what_ = ostream.str();

        // force segfault to gain stacktrace
        if(segfault_)
        {
            std::cerr << what_;
            
            int* ip = NULL;
            (*ip)++;
        }
    }
  
    Exception::Exception( const Exception& exc ) throw ()
    {
      reason_ = NULL;
      fileName_ = exc.fileName_;
      lineNum_  = exc.lineNum_;
      msg_ = exc.msg_;
      if( exc.reason_ ) {
        reason_= new Exception( *(exc.reason_) );
      }
      what_ = exc.what_;
    }
    
    Exception::~Exception() throw (){
        if(reason_)
            delete reason_;
    }


    const char * Exception::what() const throw ()
    {
        return what_.c_str();
    }
  
  
    std::string Exception::GetMsg() const {
        return msg_;
    }
 
 
    std::string Exception::GetFileName() const {
        return fileName_;
    }
  
  
    unsigned int Exception::GetLineNum() const {
        return lineNum_;
    }
                        
}
