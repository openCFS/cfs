#include <string>
#include <fstream>
#include <iostream>
 
#include "myclock.hh"

namespace CoupledField
{

  MyClock::MyClock (Char * title)
  {
 
    if (title)
      {
        Char * namefile=new Char[20];
        strcpy(namefile,title);
 
        filetime.open(strcat(namefile,".time"));
 
        if (!filetime) std::cerr<< "ERROR ("<<__FILE__<<" "<<__LINE__<<"): Cann't create  file :" 
                                << namefile << std::endl;
 
        delete namefile;
        InFile=TRUE;
 
      }
    else InFile=FALSE;
  };
 
  MyClock::~MyClock()
  {

    if (InFile) filetime.close();
 
  };

  void MyClock::ClockCount(enum status n, const std::string title)
  {
    switch(n)
      { 
      case beg:
        tm = time(0);
        ck = clock();
 
        break;
 
      case end:
        tm_tmp=time(0);
        ck_tmp=clock();
 
        if (InFile) 
          {
            if (title.empty()) filetime << "You don't specified title !" << std::endl;
 
            filetime << title << "::\t" << "Wall time: " <<
              difftime(tm_tmp, tm) << " c." << "\t" << "CPU time: " <<
              Double(ck_tmp - ck)/CLOCKS_PER_SEC << std::endl;
 
            tm=tm_tmp;
            ck=ck_tmp;
          }
        else 
          {  
            std::cout << " " << title << "::\t" << "Wall time: " <<
              difftime(tm_tmp, tm) << " c." << "\t" << "CPU time: " <<
              Double(ck_tmp - ck)/CLOCKS_PER_SEC << std::endl << std::endl;
 
            tm=tm_tmp;
            ck=ck_tmp; 
          }
 
        break;
      };
  }
} // end of namespace  
