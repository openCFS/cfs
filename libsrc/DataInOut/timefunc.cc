#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

#include "timefunc.hh"
#include "conffile.hh"
#include "WriteInfo.hh"

namespace CoupledField
{

TimeFunc :: TimeFunc(FileType * aptFileType)
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::TimeFunc" << std::endl;
#endif

  ptFileType=aptFileType;

  maxnumTimeFunc_ = 0;

  timeFncDatFile_=FALSE;
  timeFncDatFiles_=FALSE;
  std::string nametf;

  if (conf->ifgetliststr("time_data_files",fnc_names_))
      timeFncDatFiles_ = TRUE;

  else if (conf->ifget("time_data_file",nametf))
     Error("Use instead of time_data_file-command the time_data_files-command!");

  else if (conf->ifget("time_fnc",nametf))
    {
      Error("time_func-command currently not supported");

      if (timeFncDatFiles_==TRUE)
	Error("You can just use the time_data_file and the time_fnc-command simulteneously!");

      ptTimeFnc=FncReader(nametf);
      conf->get("time_arg",argTimeFnc);
      conf->get ("time_interval_fnc_a",intervalTF_a);
      conf->get ("time_interval_fnc_b",intervalTF_b);
    }             
    
  if (timeFncDatFiles_) ReadTimeFuncs();

}



void TimeFunc :: ReadTimeFunc(const std::string nametf)
{
#ifdef TRACE
  (*trace) << " entering TimeFunc :: ReadTimeFunc " << std::endl;
#endif

    std::ifstream timefile;
    timefile.open(nametf.c_str());
    if (!timefile)  
      Error("Can't open file with data of time function");

    std::string buffer;
    std::getline(timefile,buffer,'\n');

    // at present we are working only with 1 time function
    maxnumTimeFunc_ = 1; 

    maxvalTimeFunc = new Integer [maxnumTimeFunc_]; 
    
    maxvalTimeFunc[0] = 0;
    Double dummyTime, dummyVal;

    Double firstTime, firstVal;
    
    timefile >> firstTime >> firstVal;
    
    while(!timefile.eof())
    {
      timefile >> dummyTime >> dummyVal;
      maxvalTimeFunc[0]++;
    }

    if (firstVal == maxvalTimeFunc[0])
      {
	Info->Warning("Your are probably using an old time-file-style, which holds in the first line the function number and the number of values!! This line isn't needed any more.");
      }
    

    // start again reading from the top of the file
    timefile.clear();
    timefile.seekg (0, std::ios::beg);

//     timefile >> dummyTime >> dummyVal;
//     timefile >> dummyTime >> dummyVal;    
//     timefile.seekg (0, std::ios::beg);


    timeTimeFunc = new Double * [maxnumTimeFunc_];  // timeTF and valTF

    for (Integer i=0; i < maxnumTimeFunc_; i++)
      timeTimeFunc[i] = new Double[maxvalTimeFunc[i]];

    valTimeFunc = new Double * [maxnumTimeFunc_];
    for (Integer i=0; i < maxnumTimeFunc_; i++)
      valTimeFunc[i] = new Double[maxvalTimeFunc[i]];

    for (Integer i=0; i < maxnumTimeFunc_; i++)
      for (Integer j=0; j < maxvalTimeFunc[i]; j++)
	timefile >> timeTimeFunc[i][j] >> valTimeFunc[i][j];

    timefile.close();
}


void TimeFunc :: ReadTimeFuncs()
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::ReadTimeFuncs " << std::endl;
#endif

  maxnumTimeFunc_ = fnc_names_.size();

  maxvalTimeFunc = new Integer [maxnumTimeFunc_];
  timeTimeFunc   = new Double * [maxnumTimeFunc_]; 
  valTimeFunc    = new Double * [maxnumTimeFunc_];

  for (Integer i=0; i<fnc_names_.size(); i++)
    {
      std::ifstream timefile;

      timefile.open(fnc_names_[i].c_str());
      if (!timefile)  
	Error("Can't open file with data of time function");


      Double dummyTime, dummyVal;
      maxvalTimeFunc[i] = 0;

      Double firstVal;
      timefile >> dummyTime >> firstVal;
    
      while(!timefile.eof())
	{
	  timefile >> dummyTime >> dummyVal;
	  maxvalTimeFunc[0]++;
	}
      
      // +1 is due to the fact, that maxvalTimeFunc includes the first (header) line
      if (firstVal+1 == maxvalTimeFunc[0])
	{
	  std::string warnText = "Your are probably using an old time-file-style, \n";
	  warnText += "which holds in the first line the function number and the number of values!!\n";
	  warnText += "This line isn't needed any more.";

	  Info->Warning(warnText);
	}

      // start again reading from the top of the file
      timefile.clear();
      timefile.seekg (0, std::ios::beg);

      timeTimeFunc[i] = new Double[maxvalTimeFunc[i]];
      valTimeFunc[i]  = new Double[maxvalTimeFunc[i]];

      for (Integer j=0; j < maxvalTimeFunc[i]; j++)
	timefile >> timeTimeFunc[i][j] >> valTimeFunc[i][j];

      timefile.close();
    }
}


// this type needs a header with a function number and number of entries
void TimeFunc :: ReadTimeFuncOldType(const std::string nametf)
{
#ifdef TRACE
  (*trace) << " entering TimeFunc :: ReadTimeFunc " << std::endl;
#endif

    std::ifstream timefile;
    timefile.open(nametf.c_str());
    if (!timefile)  
      Error("Can't open file with data of time function");

    std::string buffer;
    std::getline(timefile,buffer,'\n');

    // at present we are working only with 1 time function
    maxnumTimeFunc_ = 1; 

    maxvalTimeFunc = new Integer [maxnumTimeFunc_];

    Integer ibuf;
    timefile >> ibuf >> maxvalTimeFunc[0];
 
    timeTimeFunc = new Double * [maxnumTimeFunc_];  // timeTF and valTF

    for (Integer i=0; i < maxnumTimeFunc_; i++)
      timeTimeFunc[i] = new Double[maxvalTimeFunc[i]];

    valTimeFunc = new Double * [maxnumTimeFunc_];
    for (Integer i=0; i < maxnumTimeFunc_; i++)
      valTimeFunc[i] = new Double[maxvalTimeFunc[i]];

    for (Integer i=0; i < maxnumTimeFunc_; i++)
      for (Integer j=0; j < maxvalTimeFunc[i]; j++)
	timefile >> timeTimeFunc[i][j] >> valTimeFunc[i][j];

    timefile.close();
}

// Double TimeFunc::TimeFuncWaveSt(const Double time)
// {
//   Double valtf=sin(2*3.1416*50000*time);
//   testtf << time << "    " << valtf << std::endl;
//   return valtf;
// }

Double TimeFunc::TimeFuncAtTime(const Double time,  const std::string fncname)
{
#ifdef TRACE
  (*trace) << " entering TimeFunc::TimeFuncAtTime " << std::endl;
#endif 
 
  if (timeFncDatFiles_) {
    return ValTimeFuncDatFile(time,fncname);
  }
  else {
    if ((intervalTF_a <= time) && (intervalTF_b >= time) )
          return ptTimeFnc(2*3.1416*argTimeFnc*time);
    else return 0;
  }
}

Double TimeFunc::ValTimeFuncDatFile(const Double time, const std::string fncname)
{
 Double help,help1,help2;
 Integer i, sizefnc, numfnc;

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

 sizefnc=maxvalTimeFunc[numfnc];

 if ( time < timeTimeFunc[numfnc][0] )
    Error("Wrong time in TimeFuncAtTime",__FILE__,__LINE__);
 
 //if time larger as defined in time function, return the last value
 if (time > timeTimeFunc[numfnc][sizefnc-1]) 
   return valTimeFunc[numfnc][sizefnc-1];
 
 for (i=0; i<sizefnc; i++)
 {
    help=timeTimeFunc[numfnc][i];
    if (time < help) break;
    if (help==time) return valTimeFunc[numfnc][i];
    if (time > help) continue;
 }
 
  help1=help-timeTimeFunc[numfnc][i-1];
  help2=((help-time)/help1)*valTimeFunc[numfnc][i-1]+
        ((time-timeTimeFunc[numfnc][i-1])/help1)*valTimeFunc[numfnc][i];
  return help2;
}
 
// ----------------- Deconstructor of TimeFunc -------------------------------
 
TimeFunc :: ~TimeFunc()
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::~TimeFunc" << std::endl;
#endif
 
  if (timeFncDatFile_) 
    {
      if (timeTimeFunc) 
	{ 
	  for (Integer i=0; i < maxnumTimeFunc_; i++ )
	    delete [] timeTimeFunc[i];
	  delete [] timeTimeFunc;
	}
      if (valTimeFunc) 
	{ 
	  for (Integer i=0; i < maxnumTimeFunc_; i++ )
	    delete [] valTimeFunc[i];
	  delete [] valTimeFunc;
	}
      if (maxvalTimeFunc) 
	delete [] maxvalTimeFunc;
    }
}

// -------------- Print TimeFunc -------------------------------------------
 
 void TimeFunc::Print(std::ostream * outfileDat) const
{
  (*outfileDat) << "------------- Print Time function ----------------" 
		<< std::endl;  for (Integer i=0; i < maxnumTimeFunc_; i++)
  {
  (*outfileDat) << " Number of time function is " << i+1 << std::endl;
   (*outfileDat) << "\t" << "time" <<
                       "\t" << "value" << std::endl;
  for (Integer ii=0; ii < maxvalTimeFunc[i]; ii++ )
    {
      (*outfileDat) << "\t" << timeTimeFunc[i][ii] <<
                       "\t" << valTimeFunc[i][ii] << std::endl;
    }
  } // end of 1st for
 
} // end of TimeFunc::Print

} // end of namespace
