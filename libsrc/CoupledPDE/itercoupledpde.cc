#include "itercoupledpde.hh"
#include "DataInOut/WriteInfo.hh"
#include "Utils/vector.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

  IterCoupledPDE::IterCoupledPDE(StdVector<BasePDE*> & PDEs,
				 StdVector<PDECoupling*> & Couplings,
				 Grid *aptgrid, 
				 BCs *aptBCs, 
				 FileType *aInFile, 
				 WriteResults *aOutFile) 
    : BaseCoupledPDE(PDEs, Couplings, 
		     aptgrid, aptBCs, aInFile, aOutFile)
  {
    ENTER_FCN( "IterCoupledPDE::IterCoupledPDE" );

    maxiter_ = 100;
    nonLinLogging_ = TRUE;
    
    // Per default all PDEs are solved
    solvePDE_.Resize(PDEs_.GetSize());
    solvePDE_.Init(TRUE);
   
    //if values are defined in conf-file, take these
    params->Get("maxNumIters", maxiter_, "iterative", "nonLinear");
    nonLinLogging_ = params->IsSet("logging", "iterative", "nonLinear");
  } 




  
  IterCoupledPDE::~IterCoupledPDE()
  {
    ENTER_FCN( "IterCoupledPDE::~IterCoupledPDE" );
  }




  void IterCoupledPDE::InitCoupling(Integer level) {

    ENTER_FCN ("IterCoupledPDE::InitCoupling" );
  
    StdVector<std::string> quantities;
    StdVector<std::string> interfaceTypes;
    StdVector<std::string> interfaceNames;
    StdVector<std::string> stopCritQuantities;
    StdVector<std::string> neighbourRegions;
    std::string normtype;
    std::string errMsg;
    Double epsilon;
    
    // vectors containing each quantity only one
    // time
    StdVector<std::string> quantitiesSorted;
    StdVector<std::string> interfaceTypesSorted;
    StdVector<StdVector<std::string> > interfaceNamesSorted;
    StdVector<std::string> nameVectorAux;

    
    // Get list for all quantities, which are listed
    // in section "nonLinear" and have a stopping
    // criterion
    params->GetList("quantity", stopCritQuantities, "nonLinear", "stopCrit");
    
    // Iterate over all PDEs
    for ( Integer iPDE = 0; iPDE < PDEs_.GetSize(); iPDE++ ) {

      quantities.Clear();
      interfaceTypes.Clear();
      interfaceNames.Clear();
      quantitiesSorted.Clear();
      interfaceTypesSorted.Clear();
      interfaceNamesSorted.Clear();
	
      // Check for coupling quantities for current PDE
      params->GetList("quantity", quantities,
		      PDEs_[iPDE]->GetName(), "coupling");

      params->GetList("type", interfaceTypes, 
		      PDEs_[iPDE]->GetName(),	"coupling");

      params->GetList("name", interfaceNames,
		      PDEs_[iPDE]->GetName(), "coupling");

      // Check if definition is consistent
      if (quantities.GetSize() != interfaceTypes.GetSize() ||
	  quantities.GetSize() != interfaceNames.GetSize() ||
	  interfaceNames.GetSize() != interfaceTypes.GetSize()) {

	errMsg  = "IterCoupledPDE::InitCoupling: Inconsistent definition ";
	errMsg += "of Coupling interfaces for PDE '";
	errMsg += PDEs_[iPDE]->GetName();
	errMsg += "'. Check your parameter file!";
	Error(errMsg.c_str(), __FILE__, __LINE__);
      }
	
      std::string typeAux;
      std::string quantityTemp;
      std::string nameAux;

      // Check if for one quantity different kinds of 
      // interfaces were specified
      for (Integer i=0; i<quantities.GetSize(); i++) {

	typeAux = interfaceTypes[i];
	quantityTemp = quantities[i];
	
	for (Integer j=0; j<quantities.GetSize(); j++) {
	  if (quantityTemp == quantities[j]) {
	    if (typeAux != interfaceTypes[j]) {
	      errMsg  = "IterCoupledPDE::InitCoupling:";
	      errMsg += "Per coupling quantity only one type (node, region";
	      errMsg += ", surface) of interface is allowed. \n You specified ";
	      errMsg += "for quantity '";
	      errMsg += quantityTemp;
	      errMsg += "' as interface '";
	      errMsg += typeAux;
	      errMsg +="' and '";
	      errMsg += interfaceTypes[j];
	      errMsg += "'. Please correct the parameter file!";
	      Error(errMsg.c_str(), __FILE__, __LINE__);
	    }
	  }
	}
      }
	
      // Now merge for each coupling quantity all interface names
      // of the same type (node, interface, region)
      Boolean found = FALSE;
      for (Integer iQuant=0; iQuant<quantities.GetSize(); iQuant++) {
	typeAux = interfaceTypes[iQuant];
	quantityTemp = quantities[iQuant];
	nameAux = interfaceNames[iQuant];
	found = FALSE;
	for (Integer iSorted=0; iSorted<quantitiesSorted.GetSize(); iSorted++){
	  if(quantityTemp == quantitiesSorted[iSorted]) {
	    interfaceNamesSorted[iSorted].Push_back(nameAux);
	    found = TRUE;
	  }
	}
	if (!found) {
	  quantitiesSorted.Push_back(quantityTemp);
	  interfaceTypesSorted.Push_back(typeAux);
	  nameVectorAux.Clear();
	  nameVectorAux.Push_back(nameAux);
	  interfaceNamesSorted.Push_back(nameVectorAux);
	}	
      }

      NormType normTypeAux;
      CouplingRegionType regionTypeAux;
      SolutionType quantityAux;

      // Construct vectors for restricted search parameter
      StdVector<std::string> keyVec;
      StdVector<std::string> attrVec;
      StdVector<std::string> valVec;
      attrVec = "", "", "quantity";

      // Get for each quantity the according stopping Criteria and normtype
      for( Integer iQuant = 0; iQuant < quantitiesSorted.GetSize(); iQuant++ ){

	epsilon = 0.0;
	normtype.clear();

	// Quantity for which we are interested in value and l2norm
	 attrVec = "", "", "quantity";
	 valVec = "", "", quantitiesSorted[iQuant];

	String2Enum(quantitiesSorted[iQuant], quantityAux);

	if ( stopCritQuantities.Find(quantitiesSorted[iQuant]) != -1 ) {

	  // First get the value
	  keyVec  = "couplingList", "nonLinear", "stopCrit", "value";
	  params->Get( keyVec, attrVec, valVec, epsilon );

	  // Now get the l2norm
	  keyVec  = "couplingList", "nonLinear", "stopCrit", "l2Norm";
	  params->Get( keyVec, attrVec, valVec, normtype );

	  // if quantity is elecForceVWP or magForceVWP, get neighbouring region
	  neighbourRegions.Clear();
	  if (quantityAux == MAG_FORCE_VWP ||
	      quantityAux == ELEC_FORCE_VWP) {
	    keyVec = "couplingList", "coupling", "neighbourRegion";
	    attrVec = "", "quantity";
	    valVec = "",  quantitiesSorted[iQuant];
	    params->GetList(keyVec, attrVec, valVec, neighbourRegions);

	   //  std::cerr << "IterCoupledPDE::InitCoupling: neighbourRegions = " << std::endl;
// 	    std::cerr << neighbourRegions << std::endl;
	  }

	}
	else {
	  epsilon = -1.0;
	  normtype = "no";
	}
	
	  String2Enum(normtype, normTypeAux);
	  String2Enum(interfaceTypesSorted[iQuant], regionTypeAux);
	  
	

	// register the interface at the according coupling-object
	Couplings_[iPDE]->AddInput(quantityAux, 
				   interfaceNamesSorted[iQuant],
				   regionTypeAux, neighbourRegions, level, 
				   epsilon, normTypeAux, Couplings_);
	norms_.Push_back(1.0);
      }
    }
    
    // Initialize each PDEs coupling terms
    for ( Integer i = 0; i < PDEs_.GetSize(); i++ ) {
      PDEs_[i]->InitCoupling(Couplings_[i]); 
    }

    // write coupling data in .info-file
    WriteCouplingInfo(*debug);
  }


  void IterCoupledPDE::DefineSolvingPDEs(StdVector<BasePDE*> & pdes)
  {
    ENTER_FCN( "IterCoupledPDE::DefineSolvingPDEs" );

    solvePDE_.Init(FALSE);
    Boolean nameFound;
    std::string errMsg;

    for (Integer i=0; i<pdes.GetSize(); i++){
      nameFound = FALSE;
      for (Integer j=0; j<PDEs_.GetSize(); j++)
	if (pdes[i]->GetName() == PDEs_[j]->GetName()) {
	  nameFound = TRUE;
	  solvePDE_[j] = TRUE;
	  break;}

      if (! nameFound){
	errMsg  = "IterCoupledPDE:DefineSolvingPDEs: The PDE with name'";
	errMsg += pdes[i]->GetName();
	errMsg += "' is not defined in current set of coupling!";
	Error(errMsg.c_str(), __FILE__, __LINE__);
      }
    }
  }



  void IterCoupledPDE::SolveStepStatic(const Integer kstep, const Double aTime,
				       const Integer level,
				       const Boolean updatesysmat ) {

    ENTER_FCN ( "entering  IterCoupledPDE::SolveStepStatic" );
  
    CFSVector *val, *oldVal;
    Integer iter = 0;
    Integer counter = 0;
    Boolean normsReached = FALSE;
    std::string quantityConv;


    while (iter < maxiter_ &&  (! normsReached))
      {
	if (nonLinLogging_)
	  {
	Info->PrintF(coupledpdename_,""); 
	Info->PrintF(coupledpdename_, " COUPLED ITERATION %i =================================", 
		     iter+1);
	  }
	
	counter = 0;
	normsReached = TRUE;
      
	for (Integer i=0; i<PDEs_.GetSize(); i++)
	  {
	    if (nonLinLogging_)
	    Info->PrintF(coupledpdename_, " Processing PDE %s", 
			 (PDEs_[i]->GetName()).c_str());

	    // Only solve current PDE, if the corresponding
	    // flag in 'solvePDE_' is set to TRUE
	    if (solvePDE_[i] == TRUE) {
	      PDEs_[i]->PreStepStatic(kstep,aTime,actlevel_,updatesysmat);
	      PDEs_[i]->CalcInputCoupling();
	      PDEs_[i]->SolveStepStatic(kstep,aTime,actlevel_,updatesysmat);
	      PDEs_[i]->PostStepStatic(kstep,aTime,actlevel_);
	      PDEs_[i]->CalcOutputCoupling();
	      
	      // Calculate Norms
	      for (Integer k=0; k<Couplings_[i]->GetNumOutputCouplings(); k++)
		{
		  Couplings_[i]->GetOutputValues(k, val);
		  Couplings_[i]->GetOutputOldValues(k, oldVal);
		  norms_[counter] = CalcNorm(Couplings_[i]->GetOutputNormType(k), *val, *oldVal);

		  if (nonLinLogging_)
		    {
		      Enum2String(Couplings_[i]->GetOutputQuantity(k), quantityConv);
		      Info->PrintF(coupledpdename_, " %s : Norm of %s = %g", 
				   (Couplings_[i]->GetPDEName()).c_str(),
				   quantityConv.c_str(), norms_[counter]);
		    }
		  
		  if (norms_[counter] > Couplings_[i]->GetOutputEpsilon(k) && 
		      Couplings_[i]->GetOutputNormType(k) != NO_NORM)
		    normsReached = FALSE;
		
		  //copy values of new solution to old one
		  *oldVal = *val;

		} // end of if 
		counter++;	      
	      }
	  }

	iter++;
	
	if(nonLinLogging_)
	Info->PrintF(coupledpdename_, "\n"); 
      }

    // now we are converged and can compute any postprocessing-quantities
    for (Integer i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->PostProcess(actlevel_);

  }



  


  void IterCoupledPDE::SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
				      const Boolean updatesysmat)
  {
    ENTER_FCN( "IterCoupledPDE::SolveStepTrans" );

    Double  steptime  = asteptime;

    Integer iter = 0;
    Boolean normsReached = FALSE;
    std::string quantityConv;

    while (iter < maxiter_ &&  (! normsReached))
      {
	if (nonLinLogging_)
	  {
	    Info->PrintF(coupledpdename_,""); 
	    Info->PrintF(coupledpdename_, " COUPLED ITERATION %i =================================", 
			 iter+1);
	  }

	Integer counter = 0;
	normsReached = TRUE;
      
	for (Integer i=0; i<PDEs_.GetSize(); i++)
	  {
	    if (nonLinLogging_)
	      Info->PrintF(coupledpdename_, " Processing PDE %s", 
			   (PDEs_[i]->GetName()).c_str());

	    // Only solve current PDE, if the corresponding
	    // flag in 'solvePDE_' is set to TRUE
	    if (solvePDE_[i] == TRUE) {
	      
	      PDEs_[i]->PreStepTrans(kstep, steptime, level, updatesysmat);
	      PDEs_[i]->CalcInputCoupling();
	      PDEs_[i]->SolveStepTrans(kstep, steptime, level, updatesysmat);
	      PDEs_[i]->PostStepTrans(kstep, steptime, level);
	      PDEs_[i]->CalcOutputCoupling();
	      
	      // Calculate Norms
	      for (Integer k=0; k<Couplings_[i]->GetNumOutputCouplings(); k++)
		{
		  CFSVector *val, *oldVal;
		  Couplings_[i]->GetOutputValues(k, val);
		  Couplings_[i]->GetOutputOldValues(k, oldVal);
		  norms_[counter] = CalcNorm(Couplings_[i]->GetOutputNormType(k), *val, *oldVal);
		  
		  if (nonLinLogging_)
		    {
		      Enum2String(Couplings_[i]->GetOutputQuantity(k), quantityConv);
		      Info->PrintF(coupledpdename_, " %s : Norm of %s = %g", 
				   (Couplings_[i]->GetPDEName()).c_str(),
				   quantityConv.c_str(), norms_[counter]);
		    }
		  if (norms_[counter] > Couplings_[i]->GetOutputEpsilon(k))
		    normsReached = FALSE;
		  
		  *oldVal = *val;
		} // end if
	      counter++;	      
	      }
	  }

	iter++;
	if (nonLinLogging_)
	  Info->PrintF(coupledpdename_, "\n"); 
      }
  }

  

void IterCoupledPDE::WriteResultsInFile(Integer stepOffset,
					Double timeOffset)
{
  ENTER_FCN( "IterCoupledPDE::WriteResultsInFile" );

  for (Integer i=0; i<PDEs_.GetSize(); i++)
    PDEs_[i]->WriteResultsInFile(stepOffset, timeOffset);
}


void IterCoupledPDE::WriteCouplingInfo(std::ostream &out)
{
  ENTER_FCN( "IterCoupledPDE::WriteCouplingInfo" );

  CFSVector *val;
  StdVector<Integer> * nodes;
  StdVector<std::string> couplingRegions;
  std::string couplingTypeAux, quantityAux, regionTypeAux, normTypeAux;

  if (!debug)
    return;

  // write information in .debug-file
  out << "=======================" << std::endl;
  out << " COUPLING INFORMATION  " << std::endl;
  out << "=======================" << std::endl;
  out << std::endl;

  for (Integer ipde=0; ipde<PDEs_.GetSize(); ipde++)
    {
      
      out << "Entering " << Couplings_[ipde]->GetPDEName() << ".InitCoupling" << std::endl;
      out << "=====================================" << std::endl;
      
      
      // Show InputCouplings
      for (Integer i=0; i<Couplings_[ipde]->GetNumInputCouplings(); i++)
	{
	  Couplings_[ipde]->GetInputNodes(i, nodes);
	  Couplings_[ipde]->GetInputValues(i,val);

	  // Convert enum-types in strings
	  Enum2String(Couplings_[ipde]->GetInputType(i), couplingTypeAux);
	  Enum2String(Couplings_[ipde]->GetInputQuantity(i), quantityAux);
	  Enum2String(Couplings_[ipde]->GetInputRegionType(i), regionTypeAux);
	  Enum2String(Couplings_[ipde]->GetInputNormType(i), normTypeAux);

	  out << std::endl;
	  out << "Input Coupling " << i+1 << ":" << std::endl;
	  out << "---------------------" << std::endl;
	  out << "Coupling Type: " << couplingTypeAux << std::endl;
	  out << "InputQuantity: " << quantityAux<< std::endl;
	  Couplings_[ipde]->GetInputRegions(i, couplingRegions);
	  out << "RegionNames: ";
	  for (Integer j=0; j<couplingRegions.GetSize()-1; j++)
	    out << couplingRegions[j] << ", ";
	  out << couplingRegions[couplingRegions.GetSize() -1];
	  out << std::endl;
	  out << "RegionType: " << regionTypeAux << std::endl;
	  out << "Number of Input coupling values: " << val->GetSize() << std::endl;
	  out << "Dof of Input Values: " << Couplings_[ipde]->GetInputDof(i) << std::endl;
	  out << "Number of Input Nodes: " << Couplings_[ipde]->GetInputNumNodes(i) << std::endl;
	  out << "Number of Input Elems: " << Couplings_[ipde]->GetInputNumElems(i) << std::endl;
	  out << "NormType: " << normTypeAux << std::endl;
	  out << "Tolerance: " << Couplings_[ipde]->GetInputEpsilon(i) << std::endl;
	  
	
	}
      out << std::endl;
      
      couplingRegions.Clear();
      // Show OutputCouplings
      nodes = 0;
      for (Integer i=0; i<Couplings_[ipde]->GetNumOutputCouplings(); i++)
	{
	  Couplings_[ipde]->GetOutputNodes(i, nodes);
	  Couplings_[ipde]->GetOutputValues(i,val);

	  // Convert enum-types in strings
	  Enum2String(Couplings_[ipde]->GetOutputType(i), couplingTypeAux);
	  Enum2String(Couplings_[ipde]->GetOutputQuantity(i), quantityAux);
	  Enum2String(Couplings_[ipde]->GetOutputRegionType(i), regionTypeAux);
	  Enum2String(Couplings_[ipde]->GetOutputNormType(i), normTypeAux);

	  out << std::endl;
	  out << "Output Coupling " << i+1 << ":" << std::endl;
	  out << "---------------------" << std::endl;
	  out << "Coupling Type: " << couplingTypeAux << std::endl;
	  out << "OutputQuantity: " << quantityAux << std::endl;
	  out << "RegionNames: ";
	  Couplings_[ipde]->GetOutputRegions(i, couplingRegions);
	  if (couplingRegions.GetSize() > 0)
	    {
	      for (Integer j=0; j<couplingRegions.GetSize()-1; j++)
		out << couplingRegions[j] << ", ";
	      out << couplingRegions[couplingRegions.GetSize() -1];
	      out << std::endl;   
	    }
	  out << "RegionType: " << regionTypeAux << std::endl;
	  out << "Number of Output coupling Values: " << val->GetSize() << std::endl;
	  out << "Dof of Output Values: " << Couplings_[ipde]->GetOutputDof(i) << std::endl;
 	  out << "Number of Output Nodes: " << Couplings_[ipde]->GetOutputNumNodes(i) << std::endl;
	  out << "Number of Output elems: " << Couplings_[ipde]->GetOutputNumElems(i) << std::endl;
	  out << "NormType: " << normTypeAux << std::endl;
	  out << "Tolerance: " << Couplings_[ipde]->GetOutputEpsilon(i) << std::endl;
	  
	}
      out << std::endl;
      
    }

}


Double IterCoupledPDE::CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval)
{
  ENTER_FCN( "IterCoupledPDE::CalcNorm" );

  // ATTENTION: Currently only working with Double-values
  // will be changed as soon as dynamic type information
  // is available

  Vector<Double> delta;
 
  Double norm, valNorm2;
  

  Vector<Double> & val_vec =\
    dynamic_cast<Vector<Double>& >(val);

  Vector<Double> & oldval_vec =\
    dynamic_cast<Vector<Double>& >(oldval);
  
  delta = val_vec - oldval_vec;

  switch (normtype)
    {
    case NO_NORM:
      return 0;
      break;
      
    case L2ABS:
      norm = delta.NormL2();
      break;

    case L2REL:
      valNorm2 =  val_vec.NormL2();
      if (valNorm2 > 0)
	norm = delta.NormL2() / valNorm2;
      else
	norm = delta.NormL2();

      break;
    }

  return norm;
}


} // end of namespace
