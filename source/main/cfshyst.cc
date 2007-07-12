// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>
#include <vector>

#include "Utils/myclock.hh"
#include "Utils/tracing.hh"
#include "DataInOut/DefineFiles/definefiles.hh"
#include "DataInOut/timefunc.hh"
#include "DataInOut/WriteInfo.hh"

#include "Driver/driver_header.hh"
#include "Domain/domain.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "Domain/GridCFS/interface_gridcfs.hh"
#include "General/environment.hh"

#include "Utils/preisach.hh"

using namespace CoupledField;


int main(int argc, char *argv[]) {

  //  std::cout << "Read in the file" << std::endl;

  char *name = argv[argc-1];
  char *filename = new char[100];
  strcpy(filename, name);

  //  std::cout << "filename:" << filename << std::endl;

  std::ifstream timefile;

  timefile.open(filename);
  if (!timefile)  
    Error("Can't open file with data of time function");

  timefile.clear(); // clear flags

  // we don't trust .eof() =)
  timefile.seekg(0,std::ios::end);
  std::string::size_type pos = 0, pos_end = timefile.tellg(), line_end_pos = 0;

  timefile.seekg(0,std::ios::beg); // start from the beginning
  std::string     buf;
  Double          timeT, valT;
  
  Vector<Double> timeTF(500); 
  Vector<Double> valTF(500); 

  Integer idx =0;
  while(pos <= pos_end) {         
    buf = "";
    std::getline(timefile,buf,'\n');
    line_end_pos = timefile.tellg();
    
    // big choice of signs for comment's
    if (buf.length() != 0 &&
        buf[0] != '#' &&
        buf[0] != '%' && 
        buf[0] != '!') 
      {
        timefile.seekg(- buf.size() - 1,std::ios::cur); // rewind
        
        timefile >> timeT >> valT ;                  
        
        valTF[idx] = valT;
        timeTF[idx] = timeT;
        idx++;
        timefile.ignore(100,'\n');
          
      }
      
    pos = timefile.tellg();  // and, where we are ?    
    
    if( pos != line_end_pos)
      std::cerr << "Wrong format of input file" << std::endl;     
  }
  
  timefile.close();

  Integer actsize = idx;
  //   for (Integer i=0; i<actsize; i++) 
  //     std::cout << timeTF[i] << " " << valTF[i] << std::endl;

  std::vector<Integer> test;

  for (Integer k=0; k<10; k++)
    test.push_back(k);

  Integer numEl = 1;
  Double  xSat  = 1;
  Double  ySat  = 1;
  Double  xCoer = 1;
  bool isVirgin = true;

  Hysteresis * hyst;
  hyst = new Preisach(numEl, xSat, ySat, xCoer, isVirgin);

  Integer elNr = 1;
  Double * Pval = new Double[actsize];
  for (Integer i=0; i<actsize; i++) {
    Pval[i] = hyst->computeValue(valTF[i], elNr);
    std::cout <<  valTF[i] << "  " << (Pval[i])*ySat << std::endl;
    //    std::cout <<  valTF[i] << "  " << (Pval[i]-0.5)*ySat << std::endl;
  }
  std::cerr << "FINISHED" << std::endl;


}


