// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <vector>

//#include "DataInOut/simInput.hh"
#include "ApproxData.hh"

namespace CoupledField
{ 

  ApproxData::ApproxData(std::string nlFileName, MaterialType matType )
  {
    matType_ = matType;
    nlFileName_ = nlFileName;

  }


  void ApproxData::ReadNlinFunc(std::string fncName)  {
  
    std::ifstream datafile;
  
    datafile.open(fncName.c_str());
    if ( !datafile ) {
      std::string str = "Failed to open file '" + fncName +
        "' suspected to contain nonlinear data";
      EXCEPTION( str );
    }
  
    datafile.clear(); // clear flags
  
    // we don't trust .eof() =)
    datafile.seekg(0,std::ios::end);
    std::string::size_type pos = 0, line_start_pos = 0,
      pos_end = datafile.tellg();
    
    datafile.seekg(0,std::ios::beg); // start from the beginning
    std::string     buf;
    Double          xval,yval;
    std::vector<Double> xx, yy;
    
    while(pos <= pos_end)
      {     
        line_start_pos = datafile.tellg();
        std::getline(datafile,buf);
      
        // big choice of signs for comment's
        if ((buf[0] != '#' || buf[0] != '%' || buf[0] != '!') && buf.size() > 0) 
          {
            datafile.seekg(line_start_pos); // rewind
          
            datafile >> xval >> yval;                  
            yy.push_back(yval);
            xx.push_back(xval);
            datafile.ignore(100,'\n');

          }
      
        pos = datafile.tellg();  // and, where we are ?    
      }
  
    datafile.close();
    numMeas_ = xx.size();

    x_.Resize(numMeas_);
    y_.Resize(numMeas_);
      
    for ( UInt k=0; k<numMeas_; k++ ) {
      x_[k] = xx[k];
      y_[k] = yy[k];
    }
  }


  void ApproxData::PerformChecksOnInputData() {

    if  ( matType_ == MAG_PERMEABILITY ) {
      //all x=H have to be greater than zero
      for ( UInt k=0; k<numMeas_; k++ ) {
        if ( x_[k] < 0.0 ) {
          std::string str = "There are values for H in file '" + nlFileName_ +
            "', that are smaller than zero!";
          EXCEPTION( str );
        }
      }
      //all y=B have to be greater than zero
      for ( UInt k=0; k<numMeas_; k++ ) {
        if ( y_[k] < 0.0 ) {
          std::string str = "There are values for B in file '" + nlFileName_ +
            "', that are smaller than zero!";
          EXCEPTION( str );
        }
      }
      
      // starting point should be (0,0)
      if ( x_[0] != 0.0 ) {
        WARN("First measured point in BH curve read from file '"
             << nlFileName_ << "' is not zero. We will add such a point!");
        
        Vector<Double> xx = x_; 
        Vector<Double> yy = y_; 
        numMeas_ +=1;
        
        x_.Resize(numMeas_);
        y_.Resize(numMeas_);
        x_[0] = 0.0;
        y_[0] = 0.0;
        
        for ( UInt k=1; k<numMeas_; k++ ) {
          x_[k] = xx[k-1];
          y_[k] = yy[k-1];
        }
      }
      
      //now check, if we have monotonically increasing x- and y-values
      for ( UInt k=0; k<numMeas_-1; k++ ) {
        Double epsX = 1.0;
        Double epsY = 1e-3;
        if ( x_[k+1] - x_[k] < epsX ) {
          std::string str = "The H-values in file '" + nlFileName_ +
            "', are not monotonically increasing (epsH = 1)!";
          EXCEPTION( str );
        }
        if ( y_[k+1] - y_[k] < epsY ) {
          std::string str = "The B-values in file '" + nlFileName_ +
            "', are not monotonically increasing (epsB = 1e-3)!";
          EXCEPTION( str );
        }
      }
    }
  }

}

