#include "itercoupledpde.hh"

namespace CoupledField
{

IterCoupledPDE::IterCoupledPDE(std::vector<BasePDE*> & PDEs,
			      std::vector<PDECoupling*> & Couplings,
			      Grid *aptgrid, 
			      BCs *aptBCs, 
			      FileType *aInFile, 
			      WriteResults *aOutFile) 
  : BaseCoupledPDE(PDEs, Couplings, aptgrid, aptBCs, aInFile, aOutFile)
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::InitCoupling" << std::endl;
#endif

} 
  

IterCoupledPDE::~IterCoupledPDE()
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::InitCoupling" << std::endl;
#endif
}

void IterCoupledPDE::InitCoupling(Integer level)
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::InitCoupling" << std::endl;
#endif
  
  std::vector<std::string> CouplingTerms;
  std::vector<std::string> NodeCouplings;
  std::vector<std::string> SubdomainCouplings;
  std::vector<std::string> Elem1DCouplings;
  std::vector<std::string> Elem2DCouplings;

  //std::cerr << "IterCoupledPDE::InitCoupling..." << std::endl;

  // Iterate over all PDEs
  for (Integer i=0; i<PDEs_.size(); i++)
    {
      CouplingTerms.clear();

      //std::cerr << "processing PDE" << PDEs_[i]->GetName() << std::endl;

      // read in "input_coupling_terms"
      if (conf->ifgetliststr("input_coupling_terms", CouplingTerms, PDEs_[i]->GetName()))
	{
	  //std::cerr << "found " << CouplingTerms.size() << " Coupling Terms" << std::endl;
	  for (Integer j=0; j<CouplingTerms.size(); j++)
	    {

	      NodeCouplings.clear();
	      SubdomainCouplings.clear();
	      Elem1DCouplings.clear();
	      Elem2DCouplings.clear();
	      // Read in node coupling terms
	      if (conf->ifgetliststr(CouplingTerms[j], NodeCouplings, PDEs_[i]->GetName(), "node_coupling"))
		{
		  //std::cerr << "found " << NodeCouplings.size() << " Node Couplings" << std::endl;
		  for (Integer k=0; k<NodeCouplings.size(); k++)
		  {
		    //std::cerr << "Adding Node CouplingTerms " << NodeCouplings[k] << std::endl;
		    Couplings_[i]->AddInput(CouplingTerms[j], NodeCouplings[k], NODES, actlevel_, Couplings_);
		  }
		}

	      // Read in subdomain coupling terms
	      if (conf->ifgetliststr(CouplingTerms[j], SubdomainCouplings, PDEs_[i]->GetName(), "subdomain_coupling"))
		for (Integer k=0; k<SubdomainCouplings.size(); k++)
		  {
		    //std::cerr << "SubDomain Coupling_size: " << SubdomainCouplings.size() << std::endl;
		    //std::cerr << "Adding SubDomain CouplingTerms " << SubdomainCouplings[k] << std::endl;
		    Couplings_[i]->AddInput(CouplingTerms[j], SubdomainCouplings[k], SUBDOMAIN, actlevel_, Couplings_);
		  }
	      
	      // Read in elem1D coupling terms
	      if (conf->ifgetliststr(CouplingTerms[j], Elem1DCouplings, PDEs_[i]->GetName(), "elem1d_coupling"))
		for (Integer k=0; k<Elem1DCouplings.size(); k++)
		  {
		    //std::cerr << "Adding Elem1D CouplingTerms " << Elem1DCouplings[k] << std::endl;
		    Couplings_[i]->AddInput(CouplingTerms[j], Elem1DCouplings[k], ELEMS1D, actlevel_, Couplings_);
		  }
	      
	      // Read in elem2D coupling terms
	      if (conf->ifgetliststr(CouplingTerms[j], Elem2DCouplings, PDEs_[i]->GetName(), "elem2d_coupling"))
		for (Integer k=0; k<Elem2DCouplings.size(); k++)
		  {
		    //std::cerr << "Adding Elem2D CouplingTerms " << Elem2DCouplings[k] << std::endl;
		    Couplings_[i]->AddInput(CouplingTerms[j], Elem2DCouplings[k], ELEMS2D, actlevel_, Couplings_);
		  }
	    }
	}
    }
  
  // Initialize each PDEs coupling terms
  for (Integer i=0; i<PDEs_.size(); i++)
    PDEs_[i]->InitCoupling(Couplings_[i]); 

  //for (Integer i=0; i<PDEs_.size(); i++)
  //  PDEs_[i]->BasePDE::InitCoupling(Couplings_[i]);
  
}


void IterCoupledPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::SolveStepStatic" << std::endl;
#endif
  
  for (Integer iter=0; iter<1; iter++)
    {

      for (Integer i=0; i<PDEs_.size(); i++)
	{
	  std::cerr << "*************************************" << std::endl;
	  std::cerr << "Processing PDE number " << i << std::endl;
	  PDEs_[i]->CalcInputCoupling();
	  PDEs_[i]->SolveStepStatic(actlevel_);
	  PDEs_[i]->CalcOutputCoupling();
	}
    }
}

void IterCoupledPDE::SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat)
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::SolveStepTrans" << std::endl;
#endif
}
  

void IterCoupledPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::WriteResultsInFile" << std::endl;
#endif

  for (Integer i=0; i<PDEs_.size(); i++)
    PDEs_[i]->WriteResultsInFile();

}



} // end of namespace
