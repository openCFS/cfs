#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>

#include "writeresults.hh"
#include "conffile.hh"

namespace CoupledField
{

WriteResults::WriteResults(const Char * const filename)
{
#ifdef TRACE
 if (trace) (*trace)<< "entering WriteResults::WriteResults()" << std::endl;
#endif
  namefile_=new Char[20];
  strcpy(namefile_,filename);

  NeedHistory_=TRUE;
  conf->getlist(nodeshist_,"history_node");

  if (nodeshist_.empty()) NeedHistory_=FALSE;

  if (NeedHistory_)
   {
     std::string S="mkdir -p history";
     system(S.c_str());

     Char * name=new Char[30];     

     std::string namedir="history/";
   
     Integer nnodhist=nodeshist_.size(); 
     historyfile=new std::ofstream[nnodhist];
     Integer i;
     for (i=0; i<nodeshist_.size(); i++) 
    {
     sprintf(name,"%s%s.%i.hist",namedir.c_str(),namefile_,nodeshist_[i]);

     historyfile[i].open(name);

     if (!historyfile[i]) 
          Error("Can't open history file",__FILE__,__LINE__);
    }
      
    lastsavetime.Resize(nnodhist);
    lastsavetime[0]=-1;

    delete [] name;
   }
}

void WriteResults::AddInHistory(const Double time, const Double val,const Integer ifile)
{ 
 lastsavetime[ifile]=time;
 historyfile[ifile] << time << "  " << val << std::endl;
}

WriteResults::~WriteResults()
{
#ifdef TRACE
  (*trace) << "entering WriteResults::~WriteResults" << std::endl;
#endif
  Integer i;
  for (i=0; i<nodeshist_.size(); i++)
     historyfile[i].close();

  if (historyfile) delete [] historyfile;
  delete [] namefile_; 

  delete output;
}

} // end of namespace
