#include "itercoupledpde.hh"
#include "DataInOut/WriteInfo.hh"

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
    (*trace) << "entering IterCoupledPDE::IterCoupledPDE" << std::endl;
#endif
 
    maxiter_ = 100;
   
    //if values are defined in conf-file, take these
    conf->ifget("maxiter", maxiter_, couplingSectionName_); // maximal number of iterations
  } 




  
  IterCoupledPDE::~IterCoupledPDE()
  {
#ifdef TRACE
    (*trace) << "entering IterCoupledPDE::~IterCoupledPDE" << std::endl;
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

    // Iterate over all PDEs
    for (Integer i=0; i<PDEs_.size(); i++)
      {
	CouplingTerms.clear();

	// read in "input_coupling_terms"
	if (conf->ifgetliststr("input_coupling_terms", CouplingTerms, PDEs_[i]->GetName()))
	  {
	    for (Integer j=0; j<CouplingTerms.size(); j++)
	      {

		NodeCouplings.clear();
		SubdomainCouplings.clear();
		Elem1DCouplings.clear();
		Elem2DCouplings.clear();

		// Read in node coupling terms
		if (conf->ifgetliststr(CouplingTerms[j], NodeCouplings, PDEs_[i]->GetName(), "node_coupling"))
		  {
		  for (Integer k=0; k<NodeCouplings.size(); k++)
		    {
		      Couplings_[i]->AddInput(CouplingTerms[j], NodeCouplings[k], NODES, actlevel_, Couplings_);
		      norms_.push_back(1.0);
		    }
		  }
		
		else if (conf->ifgetliststr(CouplingTerms[j], SubdomainCouplings, PDEs_[i]->GetName(), 
					    "subdomain_coupling"))
		  {
		    // Read in subdomain coupling terms
		    if (conf->ifgetliststr(CouplingTerms[j], SubdomainCouplings, PDEs_[i]->GetName(), 
					 "subdomain_coupling"))
		      for (Integer k=0; k<SubdomainCouplings.size(); k++)
			{
			  Couplings_[i]->AddInput(CouplingTerms[j], SubdomainCouplings[k], SUBDOMAIN, actlevel_, 
						  Couplings_);
			  norms_.push_back(1.0);
			}
		  }
		
		else if (conf->ifgetliststr(CouplingTerms[j], Elem1DCouplings, PDEs_[i]->GetName(), 
					    "elem1d_coupling"))
		  {
		    // Read in elem1D coupling terms
		    for (Integer k=0; k<Elem1DCouplings.size(); k++)
		      {
			Couplings_[i]->AddInput(CouplingTerms[j], Elem1DCouplings[k], ELEMS1D, actlevel_, 
						Couplings_);
			norms_.push_back(1.0);
		      }
		  }
		
		else if (conf->ifgetliststr(CouplingTerms[j], Elem2DCouplings, PDEs_[i]->GetName(), 
					     "elem2d_coupling"))
		  {
		    // Read in elem2D coupling terms
		    for (Integer k=0; k<Elem2DCouplings.size(); k++)
		      {
			Couplings_[i]->AddInput(CouplingTerms[j], Elem2DCouplings[k], ELEMS2D, actlevel_, 
						Couplings_);
			norms_.push_back(1.0);
		      }
		  }
		
		else
		  {
		    std::string message = "No specification for input_coupling_term defined for PDE: " 
		                          + PDEs_[i]->GetName();
		    Error(message.c_str());
		  }    

	      }
	  }
	else
	  {
	    std::string message = "No input_coupling_terms defined for PDE: " + PDEs_[i]->GetName();
	    Error(message.c_str());
	  }
	
      }
  
    // Initialize each PDEs coupling terms
    for (Integer i=0; i<PDEs_.size(); i++)
      PDEs_[i]->InitCoupling(Couplings_[i]); 
  
    // write coupling data in .info-file
    WriteCouplingInfo();
  }


  void IterCoupledPDE::SolveStepStatic(const Integer level, const Double aTime)
  {
#ifdef TRACE
    (*trace) << "entering  IterCoupledPDE::SolveStepStatic" << std::endl;
#endif
  
    Array<Double> *val, *oldVal;
    Integer iter = 0;
    Integer counter = 0;
    Boolean normsReached = FALSE;


    while (iter < maxiter_ &&  (! normsReached))
      {
	Info->PrintF(coupledpdename_,""); 
	Info->PrintF(coupledpdename_, " COUPLED ITERATION %i =================================", 
		     iter+1);
	

	counter = 0;
	normsReached = TRUE;
      
	for (Integer i=0; i<PDEs_.size(); i++)
	  {
	    Info->PrintF(coupledpdename_, " Processing PDE %s", 
			 (PDEs_[i]->GetName()).c_str());

	    PDEs_[i]->PreStepStatic(actlevel_);
	    PDEs_[i]->CalcInputCoupling();
	    PDEs_[i]->SolveStepStatic(actlevel_, aTime);
	    PDEs_[i]->PostStepStatic(actlevel_);
	    PDEs_[i]->CalcOutputCoupling();

	    // Calculate Norms
	    for (Integer k=0; k<Couplings_[i]->GetNumOutputCouplings(); k++)
	      {
		Couplings_[i]->GetOutputValues(k, val);
		Couplings_[i]->GetOutputOldValues(k, oldVal);
		norms_[counter] = CalcNorm(Couplings_[i]->GetOutputNormType(k), *val, *oldVal);

		Info->PrintF(coupledpdename_, " %s : Norm of %s = %g", 
			     (Couplings_[i]->GetPDEName()).c_str(),
			     (Couplings_[i]->GetOutputQuantity(k)).c_str(), norms_[counter]);
	      
		if (norms_[counter] > Couplings_[i]->GetOutputEpsilon(k))
		  normsReached = FALSE;
		
		*oldVal = *val;
		counter++;	      
	      }
	  }

	iter++;
	Info->PrintF(coupledpdename_, "\n"); 
      }

    // now we are converged and can compute any postprocessing-quantities
    for (Integer i=0; i<PDEs_.size(); i++)
      PDEs_[i]->PostProcess(actlevel_);

  }



  


  void IterCoupledPDE::SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
				      const Boolean updatesysmat)
  {
#ifdef TRACE
    (*trace) << "entering  IterCoupledPDE::SolveStepTrans" << std::endl;
#endif

    Integer pdenumber = 0;
    Integer nsys      = 0;
    Double  steptime  = asteptime;
    Integer stepsave  = isavebegin_-1;
    Double  dt        = firstdt_;  


    Integer iter = 0;
    Boolean normsReached = FALSE;

    while (iter < maxiter_ &&  (! normsReached))
      {
	Info->PrintF(coupledpdename_,""); 
	Info->PrintF(coupledpdename_, " COUPLED ITERATION %i =================================", 
		     iter+1);

	Integer counter = 0;
	normsReached = TRUE;
      
	for (Integer i=0; i<PDEs_.size(); i++)
	  {
	    Info->PrintF(coupledpdename_, " Processing PDE %s", 
			 (PDEs_[i]->GetName()).c_str());

	    PDEs_[i]->PreStepTrans(kstep, steptime, level, updatesysmat);
	    PDEs_[i]->CalcInputCoupling();
	    PDEs_[i]->SolveStepTrans(kstep, steptime, level, updatesysmat);
	    PDEs_[i]->PostStepTrans(level);
	    PDEs_[i]->CalcOutputCoupling();

	    // Calculate Norms
	    for (Integer k=0; k<Couplings_[i]->GetNumOutputCouplings(); k++)
	      {
		Array<Double> *val, *oldVal;
		Couplings_[i]->GetOutputValues(k, val);
		Couplings_[i]->GetOutputOldValues(k, oldVal);
		norms_[counter] = CalcNorm(Couplings_[i]->GetOutputNormType(k), *val, *oldVal);

		Info->PrintF(coupledpdename_, " %s : Norm of %s = %g", 
			     (Couplings_[i]->GetPDEName()).c_str(),
			     (Couplings_[i]->GetOutputQuantity(k)).c_str(), norms_[counter]);
	      
		if (norms_[counter] > Couplings_[i]->GetOutputEpsilon(k))
		  normsReached = FALSE;

		*oldVal = *val;
		
		counter++;	      
	      }
	  }

	iter++;
	Info->PrintF(coupledpdename_, "\n"); 
      }
  }

  

void IterCoupledPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::WriteResultsInFile" << std::endl;
#endif

  for (Integer i=0; i<PDEs_.size(); i++)
    PDEs_[i]->WriteResultsInFile();
}


void IterCoupledPDE::WriteCouplingInfo()
{

#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::WriteCouplingInfo" << std::endl;
#endif 

  Array<Double> *val, *oldval;
  std::vector<Integer> * nodes;

  if (!debug)
    return;

  // write information in .debug-file
  (*debug) << "=======================" << std::endl;
  (*debug) << " COUPLING INFORMATION  " << std::endl;
  (*debug) << "=======================" << std::endl;
  (*debug) << std::endl;

  for (Integer ipde=0; ipde<PDEs_.size(); ipde++)
    {
      
      (*debug) << "Entering " << Couplings_[ipde]->GetPDEName() << ".InitCoupling" << std::endl;
      (*debug) << "=====================================" << std::endl;
      
      // Show InputCouplings
      for (Integer i=0; i<Couplings_[ipde]->GetNumInputCouplings(); i++)
	{
	  Couplings_[ipde]->GetInputNodes(i, nodes);
	  Couplings_[ipde]->GetInputValues(i,val);
	  (*debug) << std::endl;
	  (*debug) << "Input Coupling " << i+1 << ":" << std::endl;
	  (*debug) << "---------------------" << std::endl;
	  (*debug) << "Coupling Type: " << Couplings_[ipde]->GetInputType(i) << std::endl;
	  (*debug) << "InputQuantity: " << Couplings_[ipde]->GetInputQuantity(i) << std::endl;
	  (*debug) << "Region: " << Couplings_[ipde]->GetInputRegion(i) << std::endl;
	  (*debug) << "RegionType: " << Couplings_[ipde]->GetInputRegionType(i) << std::endl;
	  (*debug) << "Size of Input Values: " << val->size() << std::endl;
	  (*debug) << "Dof of Input Values: " << Couplings_[ipde]->GetInputDof(i) << std::endl;
	  (*debug) << "Number of Input Nodes: " << Couplings_[ipde]->GetInputNumNodes(i) << std::endl;
	  (*debug) << "Number of Input Elems: " << Couplings_[ipde]->GetInputNumElems(i) << std::endl;
	  (*debug) << "NormType: " << Couplings_[ipde]->GetInputNormType(i) << std::endl;
	  (*debug) << "Tolerance: " << Couplings_[ipde]->GetInputEpsilon(i) << std::endl;
	  
	
	}
      (*debug) << std::endl;
      
      // Show OutputCouplings
      nodes = 0;
      for (Integer i=0; i<Couplings_[ipde]->GetNumOutputCouplings(); i++)
	{
	  Couplings_[ipde]->GetOutputNodes(i, nodes);
	  Couplings_[ipde]->GetOutputValues(i,val);
	  (*debug) << std::endl;
	  (*debug) << "Output Coupling " << i+1 << ":" << std::endl;
	  (*debug) << "---------------------" << std::endl;
	  (*debug) << "Coupling Type: " << Couplings_[ipde]->GetOutputType(i) << std::endl;
	  (*debug) << "OutputQuantity: " << Couplings_[ipde]->GetOutputQuantity(i) << std::endl;
	  (*debug) << "Region: " << Couplings_[ipde]->GetOutputRegion(i) << std::endl;
	  (*debug) << "RegionType: " << Couplings_[ipde]->GetOutputRegionType(i) << std::endl;
	  (*debug) << "Size of Output Values: " << val->size() << std::endl;
	  (*debug) << "Dof of Output Values: " << Couplings_[ipde]->GetOutputDof(i) << std::endl;
 	  (*debug) << "Number of Output Nodes: " << Couplings_[ipde]->GetOutputNumNodes(i) << std::endl;
	  (*debug) << "Number of Output elems: " << Couplings_[ipde]->GetOutputNumElems(i) << std::endl;
	  (*debug) << "NormType: " << Couplings_[ipde]->GetOutputNormType(i) << std::endl;
	  (*debug) << "Tolerance: " << Couplings_[ipde]->GetOutputEpsilon(i) << std::endl;
	  
	}
      (*debug) << std::endl;
      
    }

}


Double IterCoupledPDE::CalcNorm(NormType normtype, Array<Double> & val, Array<Double> & oldval)
{
#ifdef TRACE
  (*trace) << "entering  IterCoupledPDE::CalcNorm" << std::endl;
#endif

  Array<Double> delta;
  Double norm, valNorm2;

  delta = val - oldval;

  switch (normtype)
    {
    case L2ABS:
      norm = delta.normL2();
      break;

    case L2REL:
      valNorm2 =  val.normL2();
      if (valNorm2 > 0)
	norm = delta.normL2() / valNorm2;
      else
	norm = delta.normL2();

      break;
    }

  return norm;
}


} // end of namespace
