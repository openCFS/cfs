#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

#include "timefunc.hh"
#include "conffile.hh"

namespace CoupledField
{

TimeFunc :: TimeFunc(FileType * aptFileType)
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::TimeFunc" << std::endl;
#endif

  ptFileType=aptFileType;

  timeFncDatFile_=FALSE;
  std::string nametf;
//   if (conf->is_there("time_data_file")) {
//      conf->get("time_data_file",nametf);
//    timeFncDatFile_=TRUE;    
//    }
//   else {
//     if (conf->is_there("time_fnc")) {
      conf->get("time_fnc",nametf);
      ptTimeFnc=FncReader(nametf);
      conf->get("time_arg",argTimeFnc);
      conf->get ("time_interval_fnc_a",intervalTF_a);
      conf->get ("time_interval_fnc_b",intervalTF_b);
//     }             
//    else Error("There is no information about time function in conf-file");
//    }
    
  if (timeFncDatFile_) ReadTimeFunc(nametf); 
  std::cout << " boolean: " << timeFncDatFile_ << std::endl;
}

void TimeFunc :: ReadTimeFunc(const std::string nametf)
{
#ifdef TRACE
  (*trace) << " entering TimeFunc :: ReadTimeFunc " << std::endl;
#endif

    std::ifstream timefile;
    timefile.open(nametf.c_str());
    if (!timefile)  Error("Can't open file with data of time function");

    std::string buffer;
    std::getline(timefile,buffer,'\n');

    maxnumTimeFunc=1; // at present we are working only with 1 function

    maxvalTimeFunc = new Integer [maxnumTimeFunc];

    Integer ibuf;
    timefile >> ibuf >> maxvalTimeFunc[0];
 
    timeTimeFunc = new Double * [maxnumTimeFunc];  // timeTF and valTF

    Integer i,j;
    for (i=0; i < maxnumTimeFunc; i++)
      timeTimeFunc[i] = new Double[maxvalTimeFunc[i]];
    valTimeFunc = new Double * [maxnumTimeFunc];
    for (i=0; i < maxnumTimeFunc; i++)
      valTimeFunc[i] = new Double[maxvalTimeFunc[i]];
    for (i=0; i < maxnumTimeFunc; i++)
      for (j=0; j < maxvalTimeFunc[i]; j++)
	timefile >> timeTimeFunc[i][j] >> valTimeFunc[i][j];

    timefile.close();
}

// Double TimeFunc::TimeFuncWaveSt(const Double time)
// {
//   Double valtf=sin(2*3.1416*50000*time);
//   testtf << time << "    " << valtf << std::endl;
//   return valtf;
// }

Double TimeFunc::TimeFuncAtTime(const Double time, const Integer num)
{
#ifdef TRACE
  (*trace) << " entering TimeFunc::TimeFuncAtTime " << std::endl;
#endif 
 
  if (timeFncDatFile_) {
    return ValTimeFuncDatFile(time,num);
  }
  else {
    if ((intervalTF_a <= time) && (intervalTF_b >= time) )
          return ptTimeFnc(2*3.1416*argTimeFnc*time);
    else return 0;
  }
}

Double TimeFunc::ValTimeFuncDatFile(const Double time, const Integer num)
{
 Double help,help1,help2;
 Integer i,n;
 
 n=maxvalTimeFunc[num];

 if ( time < timeTimeFunc[num][0] )
    Error("Wrong time in TimeFuncAtTime",__FILE__,__LINE__);
 
 if (time > timeTimeFunc[num][n-1]) return 0;
 
 for (i=0; i<n; i++)
 {
    help=timeTimeFunc[num][i];
    if (time < help) break;
    if (help==time) return valTimeFunc[num][i];
    if (time > help) continue;
 }
 
  help1=help-timeTimeFunc[num][i-1];
  help2=((help-time)/help1)*valTimeFunc[num][i-1]+
        ((time-timeTimeFunc[num][i-1])/help1)*valTimeFunc[num][i];
  return help2;
}
 
// ----------------- Deconstructor of TimeFunc -------------------------------
 
TimeFunc :: ~TimeFunc()
{
#ifdef TRACE
  (*trace) << "entering TimeFunc::~TimeFunc" << std::endl;
#endif
 
  if (timeFncDatFile_) {
    if (timeTimeFunc) { for (Integer i=0; i < maxnumTimeFunc; i++ )
      delete [] timeTimeFunc[i];
    delete [] timeTimeFunc;
    }
    if (valTimeFunc) { for (Integer i=0; i < maxnumTimeFunc; i++ )
      delete [] valTimeFunc[i];
    delete [] valTimeFunc;
    }
    if (maxvalTimeFunc) delete [] maxvalTimeFunc;
  }
}

// -------------- Print TimeFunc -------------------------------------------
 
 void TimeFunc::Print(std::ostream * outfileDat) const
{
  (*outfileDat) << "------------- Print Time function ----------------" << std::endl;  for (Integer i=0; i < maxnumTimeFunc; i++)
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
