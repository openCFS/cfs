#include <string>
#include <fstream>
#include <iostream>

#include "writeresults.hh"
#include "conffile.hh"

namespace CoupledField
{

WriteResults::WriteResults()
{
#ifdef TRACE
 if (trace) (*trace)<< "entering WriteResults::WriteResults()" << std::endl;
#endif

  NeedHistory_=TRUE;
  history_node_=-1;
  conf->get("history_node",history_node_);

  if (history_node_==-1) NeedHistory_=FALSE; 

  if (NeedHistory_)
   { historyfile.open("history.hist");
      if (!historyfile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open history file" << std::endl;
                exit(1);}
   }
}

void WriteResults::AddInHistory(const Double time, const Double val)
{
 historyfile << time << "  " << val << std::endl;
}

WriteResults::~WriteResults()
{
#ifdef TRACE
  (*trace) << "entering WriteResults::~WriteResults" << std::endl;
#endif

 if (!output) delete output;
 if (NeedHistory_) historyfile.close();

}

} // end of namespace
