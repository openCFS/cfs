#include <string>
#include <fstream>
#include <iostream>

#include "writeresults.hh"
#include "conffile.hh"

namespace CoupledField
{

template<class Dim>
WriteResults<Dim>::WriteResults()
{
#ifdef TRACE
 if (trace) (*trace)<< "entering WriteResults<Dim>::WriteResults()" << std::endl;
#endif

  NeedHistory_=TRUE;
  history_node_=-1;
  conf->get("history_node",history_node_);
  std::cout << history_node_ << std::endl;

  if (history_node_==-1) NeedHistory_=FALSE; 

  std::cout << " ok " << std::endl;

  if (NeedHistory_)
   { historyfile.open("history.hist");
      if (!historyfile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open history file" << std::endl;
                exit(1);}
   }
}

template<class Dim>
void WriteResults<Dim>::AddInHistory(const Double time, const Double val)
{
 historyfile << time << "  " << val << std::endl;
}

} // end of namespace
