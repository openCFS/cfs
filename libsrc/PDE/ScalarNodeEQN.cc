#include <fstream>
#include <iostream>
#include <string>

#include "ScalarNodeEQN.hh"

namespace CoupledField
{

ScalarNodeEQN :: ScalarNodeEQN(Grid * aptgrid, BCs *aptBCs, std::vector<std::string>& asubdoms, 
			       std::vector<std::string>& bcs_hd, Integer actlevel)
  : BaseEQN(aptgrid,aptBCs,asubdoms,actlevel)
{
#ifdef TRACE
  (*trace) << "entering ScalarNodeEQN::ScalarNodeEQN" << std::endl;
#endif

  bcs_hd_ = bcs_hd;
  InitEQN();
  SetHomoDBCs();
  ApplyEqnNrs();
}

ScalarNodeEQN :: ~ScalarNodeEQN()
{
#ifdef TRACE
  (*trace) << "entering ScalarNodeEQN::~ScalarNodeEQN" << std::endl;
#endif

}


void ScalarNodeEQN :: InitEQN()
{
#ifdef TRACE
  (*trace) << "entering ScalarNodeEQN::InitEQN" << std::endl;
#endif

  EQN_.Resize(ptgrid_->GetMaxnumnodes(actlevel_));
  EQN_.Init(0);

  std::vector<Elem*> SD;

  // Iterate over Subdomains
  for (Integer i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(SD,subdoms_[i],actlevel_);
      // Iterate over all elements in subdomain
      for (Integer j=0; j<SD.size(); j++)
	{
	  // Iterate over all element nodes
	  for (Integer NumNodes=0; NumNodes<SD[j]->connect.size(); NumNodes++)
	    {
	      // Check if node was already assigned
	      if (EQN_[SD[j]->connect[NumNodes] - 1] == 0)
		{
		  EQN_[SD[j]->connect[NumNodes] - 1] = 1;
		  EQN2MeshNode_.push_back(SD[j]->connect[NumNodes]);
		}
	    }
	}
    }

  numPDENodes_ = EQN2MeshNode_.size();

}

void ScalarNodeEQN :: SetHomoDBCs()
{
#ifdef TRACE
  (*trace) << "entering ScalarNodeEQN::SetHomoDBCs" << std::endl;
#endif

  std::list<Integer> nodes;
  Integer node;

  for (Integer i=0; i<bcs_hd_.size(); i++)
    {  
      nodes=ptbcs_->GetNodesLevel(bcs_hd_[i]);
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	EQN_[*p-1] = 0;
    }
}

void ScalarNodeEQN :: ApplyEqnNrs()
{
#ifdef TRACE
  (*trace) << "entering ScalarNodeEQN::ApplyEqnNrs" << std::endl;
#endif

  numEqns_ = 0;
  for (Integer i=0; i<EQN_.size(); i++)
    if (EQN_[i] > 0)
      {
	numEqns_++;
	EQN_[i] = numEqns_;
      }

}


void ScalarNodeEQN :: Mesh2Eqn(Vector<Integer>& Eqns, Vector<Integer>& nodes)
{
#ifdef TRACE
  (*trace) << "entering ScalarNodeEQN::Mesh2Eqn" << std::endl;
#endif

  Eqns.Resize(nodes.size());
  
  for (Integer i=0; i<Eqns.size(); i++) 
     Eqns[i] = EQN_[nodes[i]-1];

}


void ScalarNodeEQN :: Print()
{
#ifdef TRACE
  (*trace) << "entering ScalarNodeEQN::Print" << std::endl;
#endif

#ifdef DEBUG
  (*debug) << std::endl << "PDE-nodes: " << numPDENodes_ << std::endl;
  (*debug) << "EQN-Array:" << std::endl;
  for (Integer i=0; i<EQN_.size(); i++)
    (*debug) << i+1 << " " << EQN_[i] << std::endl;

#endif
}

} // end of namespace
