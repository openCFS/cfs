// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <string>
#include <fstream>
#include "math.h"

#include "freqfunc.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  FreqFunc::FreqFunc() {


    maxnumFF_  = 0;

    // query parameter file
    std::string nameFreqFile;
    param->Get("sequenceStep") 
      ->Get("analysis")->Get("harmonic")->Get("freqDataFile")
      ->GetValue("name", nameFreqFile );
    fnc_names_.Push_back( nameFreqFile );
    
  }

  void FreqFunc::ReadFreqFuncs(UInt nodeNumber)
  {

    maxnumFF_ =  fnc_names_.GetSize();

    valFF_    =  new std::list<Double>[maxnumFF_];
    phaseFF_    =  new std::list<Double>[maxnumFF_];
    freqFF_     =  new std::list<Double>[maxnumFF_]; 
    std::string errMsg;
    std::string nodalFilename;
    std::string pathAndFilename;

    // loop over freq-fncs
    for ( UInt i = 0; i < maxnumFF_; i++ ) {

      std::ifstream freqfile;

  
      nodalFilename=fnc_names_[i];
      nodalFilename.append( lexical_cast<std::string>( nodeNumber ) );
      pathAndFilename="nodalSrcs/" + nodalFilename;
      freqfile.open( pathAndFilename.c_str() );
      if ( !freqfile ) {
        EXCEPTION( "Failed to open file '" << nodalFilename
                 << "' in subdir 'nodalSrcs/' containing nodal freq. function");
      }

      freqfile.clear(); // clear flags

      // we don't trust .eof() =)
      freqfile.seekg(0,std::ios::end);
      std::string::size_type pos = 0, pos_end = freqfile.tellg(),
        line_end_pos , line_start_pos = 0;

      freqfile.seekg(0,std::ios::beg); // start from the beginning
      std::string     buf;
      Double          freqF, valF, phaseF;

      while( pos <= pos_end ) {         
        buf = "";
        line_start_pos = freqfile.tellg();
        std::getline(freqfile,buf,'\n');
        line_end_pos = freqfile.tellg();
        
        // big choice of signs for comment's
        if (buf.length() != 0 &&
            buf[0] != '#' &&
            buf[0] != '%' && 
            buf[0] != '!') {
          freqfile.seekg(line_start_pos); // rewind
          freqfile >> freqF >> valF >> phaseF;                    

          phaseFF_[i].push_back(phaseF);
          valFF_[i].push_back(valF);
          freqFF_[i].push_back(freqF);
          freqfile.ignore(100,'\n');
          
        }

        pos = freqfile.tellg();  // and, where we are ?    
          
        if( pos != line_end_pos) {
          EXCEPTION("The freq data file '" << fnc_names_[i]
                                           << "' is not correctly formatted.\n"
                                           << "Please correct it!");
        }

      }
      
      freqfile.close();
        
    } // loop over fncs
  }

  StdVector<Double> FreqFunc::NodalFreqFuncAtFreq(const Double freq,  
                                                  const std::string fncname, UInt nodeNumber)
  {
  
    Integer numfnc;
    StdVector<Double> Amp_Phase_atF;
    Amp_Phase_atF.Resize(2);

    //if name of freq function not defined, then a constant freq function with value
    //1.0 is assumed 
    if (fncname == "none")
      {
        WARN( "Using ampl=1.0 and phase=0.0 since fncname was not defined" );
        Amp_Phase_atF[0]=1.0;
        Amp_Phase_atF[1]=0.0;
        return Amp_Phase_atF;
      }
    
    //get correct freq function
    numfnc = fnc_names_.Find(fncname);

    if (numfnc == -1)
      {
      EXCEPTION("Freq Function '" << fncname
      << "' is not defined within 'freqDataFile' "
      << "element in section 'Harmonic'!");
      }


    // Read the freqfile for the specific node
    ReadFreqFuncs(nodeNumber);


    if (freq < *freqFF_[numfnc].begin() )
      EXCEPTION("Wrong freq in FreqFuncAtFreq");

    //if freq larger as defined in freq function, return the last value
    if (freq > freqFF_[numfnc].back())
      {
        Amp_Phase_atF[0]=valFF_[numfnc].back();
        Amp_Phase_atF[1]=phaseFF_[numfnc].back();
        return Amp_Phase_atF;
      }
    

    //loop over array with freq
    Double                              hlp=0 , hlp1, hlp2, hlp3, hlp_prev=0;
    Double                              v_prev =0, phase_prev =0;
    std::list<Double>::const_iterator   p,v,phase,p_prev;

    for ( p=freqFF_[numfnc].begin(), v=valFF_[numfnc].begin(), phase=phaseFF_[numfnc].begin(); 
          p!= freqFF_[numfnc].end(); p++,v++,phase++)
      {

        hlp_prev = hlp;
        hlp      = *p;

        if (freq == hlp)
          {
            Amp_Phase_atF[0]=*v;
            Amp_Phase_atF[1]=*phase;
            return Amp_Phase_atF;
          }
        if (freq < hlp) break;
     
        v_prev   = *v;
        phase_prev = *phase;
      }
 
    // linear interpolation
    hlp1 = hlp - hlp_prev;
    hlp2 = ((hlp-freq)/hlp1)*v_prev+((freq-hlp_prev)/hlp1)*(*v);
    hlp3 = ((hlp-freq)/hlp1)*phase_prev+((freq-hlp_prev)/hlp1)*(*phase);
    
    Amp_Phase_atF[0]=hlp2;
    Amp_Phase_atF[1]=hlp3;
    return Amp_Phase_atF;

  }

  FreqFunc::~FreqFunc()
  {
 
    if (maxnumFF_) 
      {
        UInt i;
        for(i=0; i<maxnumFF_; i++)
          {
            freqFF_[i].clear();
            valFF_[i].clear();
            phaseFF_[i].clear();
          }

        delete [] phaseFF_;
        delete [] valFF_;
        delete [] freqFF_;

      }
  }

} // end of namespace
