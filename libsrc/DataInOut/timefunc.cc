#include <iostream>
#include <string>
#include <fstream>

#include "timefunc.hh"
#include "ParamHandling/ConfFile.hh"
#include "WriteInfo.hh"
#include "ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

TimeFunc :: TimeFunc(FileType * aptFileType)
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::TimeFunc" << std::endl;
#endif

  ptFileType = aptFileType;
  maxnumTF_  = 0;

  timeFncDatFiles_=FALSE;
  std::string nametf;

#ifndef XMLPARAMS
  if (conf->ifgetliststr("time_data_files",fnc_names_))
      timeFncDatFiles_ = TRUE;

  else if (conf->ifget("time_data_file",nametf))
     Error("Use instead of time_data_file-command the time_data_files-command!");

  else if (conf->ifget("time_fnc",nametf))
      Error("time_func-command currently not supported");
    
  if (timeFncDatFiles_) ReadTimeFuncs();
#else
  params->GetList( "name", fnc_names_, "transient", "timeDataFile" );
#endif

}

void TimeFunc :: ReadTimeFuncs()
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::ReadTimeFuncs " << std::endl;
#endif

  maxnumTF_ =  fnc_names_.size();

  valTF_    =  new std::list<Double>[maxnumTF_];
  timeTF_     =  new std::list<Double>[maxnumTF_]; 



  // loop over time-fncs
  for (Integer i=0; i< maxnumTF_; i++)
    {
      std::ifstream timefile;

      timefile.open(fnc_names_[i].c_str());
      if (!timefile)  
	Error("Can't open file with data of time function");

      timefile.clear(); // clear flags

      // we don't trust .eof() =)
      timefile.seekg(0,std::ios::end);
      std::string::size_type pos = 0, pos_end = timefile.tellg();

      timefile.seekg(0,std::ios::beg); // start from the beginning
      std::string     buf;
      Double          timeT, valT;
      while(pos <= pos_end)
	{	  
	  getline(timefile,buf);

	  // big choice of signs for comment's
	  if (buf[0] != '#' || buf[0] != '%' || buf[0] != '!') 
	    {
	      timefile.seekg(- buf.size() - 1,std::ios::cur); // rewind

	      timefile >> timeT >> valT ;		     
 
	      valTF_[i].push_back(valT);
	      timeTF_[i].push_back(timeT);
	      timefile.ignore(100,'\n');
	  
	    }

	  pos = timefile.tellg();  // and, where we are ?    
	}
      
      timefile.close();

    } // loop over fncs
}

Double TimeFunc::TimeFuncAtTime(const Double time,  const std::string fncname)
{
#ifdef TRACE
  (*trace) << " entering TimeFunc::TimeFuncAtTime " << std::endl;
#endif 
  
  Integer     numfnc;
  Integer     i;
 
 //if name of time function not defined, than a constant time function with value
 //1.0 is assumed 
 if (fncname == "---not-defined--")
   return 1.0;

 //get correct time function
 numfnc = -1;
 for (i=0; i<fnc_names_.size(); i++)
   if (fnc_names_[i] == fncname) numfnc = i;

 if (numfnc == -1)
   {
     std::string fncstring = "Time Function " + fncname + " not defined within time_data_files-command";
     Error(c_string(fncstring));
   }

 if (time < *timeTF_[numfnc].begin() )
   Error("Wrong time in TimeFuncAtTime",__FILE__,__LINE__);

 //if time larger as defined in time function, return the last value
 if (time > timeTF_[numfnc].back()) 
   return valTF_[numfnc].back();

 //loop over array with time
 Double                              hlp , hlp1, hlp2, hlp_prev=0;
 Double                              v_prev =0;
 std::list<Double>::const_iterator   p,v,p_prev;

 for ( p=timeTF_[numfnc].begin(), v=valTF_[numfnc].begin(); p!= timeTF_[numfnc].end(); p++,v++)
   {

     hlp_prev = hlp;
     hlp      = *p;

     if (time == hlp) return *v;
     if (time < hlp) break;
     
     v_prev   = *v;
   }
 
 // linear interpolation
 hlp1 = hlp - hlp_prev;
 hlp2 = ((hlp-time)/hlp1)*v_prev+((time-hlp_prev)/hlp1)*(*v);
 return hlp2;
 
}

TimeFunc :: ~TimeFunc()
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::~TimeFunc" << std::endl;
#endif
 
  if (maxnumTF_) 
    {
      Integer i;
      for(i=0; i<maxnumTF_; i++)
	{
	  timeTF_[i].clear();
	  valTF_[i].clear();
	}

      delete [] valTF_;
      delete [] timeTF_;

    }
}

 void TimeFunc::Print(std::ostream * outfileDat) const
{
  (*outfileDat) << "------------- Print Time function ----------------" 
		<< std::endl;
  for (Integer i=0; i < maxnumTF_; i++)
  {
  (*outfileDat) << " Number of time function is " << i+1 << std::endl;
   (*outfileDat) << "\t" << "time" <<
                       "\t" << "value" << std::endl;
   std::list<Double>::const_iterator   t,v;
  for (t=valTF_[i].begin(),v=valTF_[i].end(); t!=valTF_[i].end(); t++,v++ )
      (*outfileDat) << "\t" << *t << "\t" << *v << std::endl;
  } // 
 
} // end of TimeFunc::Print

} // end of namespace
