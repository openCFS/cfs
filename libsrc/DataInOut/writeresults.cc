#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>

#include "writeresults.hh"
#include "conffile.hh"
#include <DataInOut/AnsysFile/ansysfile.hh>
#include <DataInOut/WriteInfo.hh>

namespace CoupledField
{

WriteResults::WriteResults(const Char * const filename, Boolean withHistory, FileType * const aInFile)
  :NeedHistory_(withHistory)
{
#ifdef TRACE
 if (trace) (*trace)<< "entering WriteResults::WriteResults()" << std::endl;
#endif
  namefile_=new Char[20];
  strcpy(namefile_,filename);

  ascii_=TRUE;

  historyfile=NULL;

  pt2Inputfile_ = aInFile;

  if (withHistory) 
    InitHistoryFiles();
}





void WriteResults::AddInHistory(const Double time, const Double val,const Integer ifile)
{ 
 lastsavetime[ifile]=time;
 historyfile[ifile] << time << "  " << val << std::endl;
}


void WriteResults::AddVecInHistory(const Double time, const std::vector<Double> val,const Integer ifile)
{ 
 lastsavetime[ifile]=time;
 historyfile[ifile] << time;
 
 for (Integer i=0; i<val.size(); i++)
   historyfile[ifile] <<  "  " << val[i];
 
 historyfile[ifile] <<  std::endl;
}


WriteResults::~WriteResults()
{
#ifdef TRACE
  (*trace) << "entering WriteResults::~WriteResults" << std::endl;
#endif
  
  if (historyfile) {
  Integer i;
  for (i=0; i<nodeshist_.size(); i++)
     historyfile[i].close();
  }
  if (historyfile) delete [] historyfile;
  delete [] namefile_; 
}




void WriteResults::ReadSaveNodes()
{
#ifdef TRACE
  (*trace)<< "entering WriteResults::ReadSaveNodes" << std::endl;
#endif
 std::vector<std::string> historyList; 
 std::list<Integer> * histNodes;
 
 conf->ifgetliststr("save_nodes", historyList);  

 if (historyList.size())
   histNodes = new std::list<Integer>[historyList.size()];

 pt2Inputfile_->ReadSaveNodes(histNodes, historyList);

 if (historyList.size())
   {
     Info->PrintVec("Area names, in which save nodes are stored:", historyList);
     
     for(int i=0; i<historyList.size(); i++)
       for (std::list<Integer>::const_iterator p=histNodes[i].begin(); p!=histNodes[i].end(); p++)
	 nodeshist_.push_back(*p);
   }
}






void WriteResults::InitHistoryFiles()
{
#ifdef TRACE
 if (trace) (*trace)<< "entering WriteResults::InitHistoryFiles()" << std::endl;
#endif

 // read nodes "by name" from the config-file command "save_nodes"
 ReadSaveNodes();
 
 std::vector<Integer> nodesTmp;

 conf->getlist(nodesTmp,"history_node");


 for (int i=0; i < nodesTmp.size(); i++)
   // there are allready elements in nodeshist_
   nodeshist_.push_back(nodesTmp[i]);


 if (nodeshist_.empty()) 
   NeedHistory_=FALSE;
 else
   {
     std::string S="mkdir -p history";
     system(S.c_str());

     Char * name=new Char[30];     

     std::string namedir="history/";
   
     Integer nnodhist=nodeshist_.size(); 
     historyfile=new std::ofstream[nnodhist];


     // write save nodes to info file
     Info->PrintVec("List of node numbers, in which results have to be saved:", nodeshist_);     

     for (Integer i=0; i<nodeshist_.size(); i++) 
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

void WriteResults::WriteSolMatrix(Grid * ptgrid, const Integer level, const Vector<Double> sol, 
				  const std::string matFileName, const Integer nrDofs)
{
  //get and write number of nodes on the level
  Integer numnodes=ptgrid->GetMaxnumnodes(level);
  Integer dim=ptgrid->GetDim();
  Integer i;

  std::ofstream * matrixOut = new std::ofstream(matFileName.c_str());
  // use scientific output, formatting is much better
  matrixOut->setf(std::ios::scientific, std::ios::floatfield);
  
  if (dim==2) 
    {
      Point<2> point;
      
      // write x,y,z-coordinate
      for (i=0; i<numnodes; i++)
	{
	  ptgrid->GetCoordinateNode(i,level,point);
	  (*matrixOut) << " \t" << point[0] << " \t" << point[1] << " \t" << 0 << " \t";
	    
	  for (Integer actDof =0; actDof < nrDofs; actDof++)
	    (*matrixOut) << sol[i*nrDofs + actDof] << " \t";
	  (*matrixOut) << std::endl;
	}
    }
  else 
    {
      Point<3> point;
      
      // write x,y,z-coordinate
      for (i=0; i<numnodes; i++)
	{
	  ptgrid->GetCoordinateNode(i,level,point);
	  (*matrixOut) << " \t" << point[0] << " \t" << point[1] << " \t" << point[2] << " \t";

	  for (Integer actDof =0; actDof < nrDofs; actDof++)
	    (*matrixOut) << sol[i*nrDofs + actDof] << " \t";
	  (*matrixOut) << std::endl;
	} 
    }
}

} // end of namespace
