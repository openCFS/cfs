#include <iostream>

#include "exception.hh"

namespace CoupledField {

    bool Exception::segfault_ = false;

    Exception::Exception( const Exception* reason,
               const char* const fileName, 
               const unsigned int lineNum,
               const char* const message) throw () :
        fileName_(fileName),
        lineNum_(lineNum),
        msg_(message)
    {
        if(reason)
            reason_ = new Exception(*reason);
        else
            reason_ = NULL;

        std::ostringstream ostream;
        
        ostream << msg_;

        if( fileName_ != "" ) {
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
        fileName_ = exc.fileName_;
        lineNum_  = exc.lineNum_;
        msg_ = exc.msg_;
        reason_= exc.reason_;
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
