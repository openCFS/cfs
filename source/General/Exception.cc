#include <iostream>

#include "Exception.hh"

namespace CoupledField {
    // Initialization of static members
    bool Exception::segfault_ = false;
    std::function<void(Exception&)> Exception::exCallback_;
    std::function<void(Exception&)> Exception::warnCallback_;
    

    /** this is the constructor for the EXCEPTION macro */
    Exception::Exception( const Exception* reason,
               const char* const fileName, 
               const unsigned int lineNum,
               const char* const message,
               SeverityType severity) throw ()
    {
        init(reason, fileName, lineNum, message, severity);
    }

    /** this is the constructor to be called manually */
    Exception::Exception( const std::string& message,
               const char* const fileName, 
               const unsigned int lineNum,
               SeverityType severity) throw ()
    {
        init(NULL, fileName, lineNum, message.c_str(), severity);
    }

    Exception::Exception(const std::string& message, const Exception& reason) throw ()
    {
      init(&reason, "NO_FILENAME", 0, message.c_str(), EXCEPTION);
    }


    void Exception::init(const Exception* reason,
               const char* const fileName, 
               const unsigned int lineNum,
               const char* const message,
               SeverityType severity) throw () 
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
            ostream << "\n\nIn file '" << fileName_ << "'";
        }
        
        if( lineNum_ != 0 ) {
            ostream << " at line " << lineNum_ << std::endl;
        }
        
        //ostream << "\n";
        
        if(reason_)
        {
            ostream << reason_->what();
        }

        what_ = ostream.str();

        switch(severity)
        {
          case EXCEPTION:
            
            // execute callback method
            if( exCallback_ )
              exCallback_(*this);
            
            // force segfault to gain stacktrace
            if( segfault_) {
              std::cerr << what_ << std::endl;
              std::cerr.flush();
            
              int* ip = NULL;
              (*ip)++;
            }
            break;
          case WARNING:

            // execute callback method
            if( warnCallback_ )
              warnCallback_(*this);
            
//            
//            // This is pretty hard-coded.
//            // In general we should think of a way
//            // to register some sort of call-back mechanism which
//            // gets called in case of en EXCEPTION of WARNING event.
//            // This way we do not need to include the ParamNode library 
//            // here
//            std::string msg = message;
//            boost::trim(msg);
//            std::cerr << "\n "
//                << fg_yellow << "WARNING:" << fg_reset << "\n "
//                << msg << endl;
//            
            // ahauck, 2010-03-15
            // The following section is commented out, as it introduces a 
            // dependency of the exception class to the ParamNode class.
            // This in addition creates a dependency of the h5tool to the
            // libparamh ..
            //            PtrParamNode out = info->Get(ParamNode::WARNING)->Get("warning", ParamNode::APPEND);
           //            
           //            out->Get("lineNum")->SetValue(lineNum);
           //            out->Get("fileName")->SetValue(fileName);
           //            out->Get("message")->SetValue(message);
            
            break;
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

    Exception::~Exception() {
        if(reason_)
            delete reason_;
    }


    void Exception::SetCallbackEx(std::function<void (Exception& x)> cb)
    {
    exCallback_ = cb;
    }
    
    void Exception::SetCallbackWarn(std::function<void (Exception& x)> cb)
    {
      warnCallback_ = cb;
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
