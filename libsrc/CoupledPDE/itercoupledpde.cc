#include "itercoupledpde.hh"

#include "pdecoupling.hh"

#include "DataInOut/WriteInfo.hh"
#include "PDE/SinglePDE.hh"
#include "Driver/iterSolveStep.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

  IterCoupledPDE::IterCoupledPDE(StdVector<StdPDE*> & PDEs,
                                 StdVector<SinglePDE*> & singlePDEs,
				 StdVector<PDECoupling*> & Couplings,
				 std::string sequenceTag) 
    : BasePDE()
  {
    ENTER_FCN( "IterCoupledPDE::IterCoupledPDE" );
    
    PDEs_       = PDEs;
    singlePDEs_ = singlePDEs;
    Couplings_  = Couplings;
    
    NumPDEs_ = PDEs.GetSize();
    sequenceTag_ = sequenceTag;
    solveStep_ = NULL;
    
    
    // get analysis type
    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    
    //  if (analysis=="static") 
    // analysype_ = STATIC;
    // else if (analysis=="transient")
    //   analysistype_ = TRANSIENT;
    //else
    // Error("Analysis Type not supported",__FILE__,__LINE__);
    
    pdename_ = "CoupledPDE: ";


    for (UInt actPDE=0; actPDE < PDEs.GetSize()-1; actPDE++)
      pdename_ += PDEs[actPDE] -> GetName() + "+";
    
    pdename_ += PDEs[PDEs.GetSize()-1] -> GetName();
      
    maxiter_ = 100;
    nonLinLogging_ = TRUE;
    std::string loggingString;
    
    // Per default all PDEs are solved
    solvePDE_.Resize(PDEs_.GetSize());
    solvePDE_.Init(TRUE);
    
    //if values are defined in conf-file, take these
    StdVector<std::string> keyVec, attrVec, valVec;
    keyVec = "couplingList", "iterative", "nonLinear", "maxNumIters";
    attrVec = "tag", "", "";
    valVec = sequenceTag_, "", "";
    
    params->Get( keyVec, attrVec, valVec, maxiter_ );
    keyVec = "couplingList", "iterative", "nonLinear", "logging";
    params->Get( keyVec, attrVec, valVec, loggingString  );
    
    if (loggingString == "yes")
      nonLinLogging_ = TRUE;
    else
      nonLinLogging_ = FALSE;
    
  } 




  
  IterCoupledPDE::~IterCoupledPDE()
  {
    ENTER_FCN( "IterCoupledPDE::~IterCoupledPDE" );

    // delete solveStep-object
    delete solveStep_;

    // delete coupling objects
    for (UInt i=0; i<Couplings_.GetSize(); i++)
      if (Couplings_[i] != NULL)
	delete Couplings_[i];

    Couplings_.Clear();
  }




  void IterCoupledPDE::InitCoupling() {

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

    // vectors for restricted parameter search
    StdVector<std::string> keyVec, valVec, attrVec;
    
    // Get list for all quantities, which are listed
    // in section "nonLinear" and have a stopping
    // criterion
    keyVec = "couplingList", "iterative", "nonLinear", "stopCrit", "quantity";
    attrVec = "tag", "", "", "";
    valVec = sequenceTag_, "", "", "";
    params->GetList( keyVec, attrVec, valVec, stopCritQuantities);
    
    // Iterate over all PDEs
    for ( UInt iPDE = 0; iPDE < singlePDEs_.GetSize(); iPDE++ ) {

      quantities.Clear();
      interfaceTypes.Clear();
      interfaceNames.Clear();
      quantitiesSorted.Clear();
      interfaceTypesSorted.Clear();
      interfaceNamesSorted.Clear();
	
      // Check for coupling quantities for current PDE
      attrVec = "tag", "", "", "";
      valVec = sequenceTag_, "", "", "";
      
      keyVec = "couplingList", "iterative", singlePDEs_[iPDE]->GetName(), 
	       "coupling", "quantity";
      params->GetList( keyVec, attrVec, valVec, quantities);

      keyVec = "couplingList", "iterative", singlePDEs_[iPDE]->GetName(),
               "coupling", "type";
      params->GetList( keyVec, attrVec, valVec, interfaceTypes);

      keyVec = "couplingList", "iterative", singlePDEs_[iPDE]->GetName(),
               "coupling", "name";
      params->GetList( keyVec, attrVec, valVec, interfaceNames);

      // Check if definition is consistent
      if (quantities.GetSize() != interfaceTypes.GetSize() ||
	  quantities.GetSize() != interfaceNames.GetSize() ||
	  interfaceNames.GetSize() != interfaceTypes.GetSize()) {

	errMsg  = "IterCoupledPDE::InitCoupling: Inconsistent definition ";
	errMsg += "of Coupling interfaces for PDE '";
	errMsg += singlePDEs_[iPDE]->GetName();
	errMsg += "'. Check your parameter file!";
	Error(errMsg.c_str(), __FILE__, __LINE__);
      }
	
      std::string typeAux;
      std::string quantityTemp;
      std::string nameAux;

      // Check if for one quantity different kinds of 
      // interfaces were specified
      for (UInt i=0; i<quantities.GetSize(); i++) {

	typeAux = interfaceTypes[i];
	quantityTemp = quantities[i];
	
	for (UInt j=0; j<quantities.GetSize(); j++) {
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
      for (UInt iQuant=0; iQuant<quantities.GetSize(); iQuant++) {
	typeAux = interfaceTypes[iQuant];
	quantityTemp = quantities[iQuant];
	nameAux = interfaceNames[iQuant];
	found = FALSE;
	for (UInt iSorted=0; iSorted<quantitiesSorted.GetSize(); iSorted++){
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
      for( UInt iQuant = 0; iQuant < quantitiesSorted.GetSize(); iQuant++ ){

	epsilon = 0.0;
	normtype.clear();

	// Quantity for which we are interested in value and l2norm
	 attrVec = "tag", "", "","quantity";
	 valVec = sequenceTag_, "" ,"", quantitiesSorted[iQuant];

	String2Enum(quantitiesSorted[iQuant], quantityAux);

	if ( stopCritQuantities.Find(quantitiesSorted[iQuant]) != -1 ) {

	  // First get the value
	  keyVec  = "couplingList", "iterative", "nonLinear", "stopCrit", "value";
	  params->Get( keyVec, attrVec, valVec, epsilon );

	  // Now get the l2norm
	  keyVec  = "couplingList", "iterative", "nonLinear", "stopCrit", "l2Norm";
	  params->Get( keyVec, attrVec, valVec, normtype );

	  // if quantity is elecForceVWP or magForceVWP, get neighbouring region
	  neighbourRegions.Clear();
	  if (quantityAux == MAG_FORCE_VWP ||
	      quantityAux == ELEC_FORCE_VWP) {
	    keyVec = "couplingList", "iterative", "coupling", "neighbourRegion";
	    attrVec = "tag", "", "quantity";
	    valVec = sequenceTag_, "",  quantitiesSorted[iQuant];
	    params->GetList(keyVec, attrVec, valVec, neighbourRegions);

	  }

	}
	else {
	  epsilon = 0.0;
	  normtype = "no";
	}

	  String2Enum(normtype, normTypeAux);
	  String2Enum(interfaceTypesSorted[iQuant], regionTypeAux);
	  
	

	// register the interface at the according coupling-object
	Couplings_[iPDE]->AddInput(quantityAux, 
				   interfaceNamesSorted[iQuant],
				   regionTypeAux, neighbourRegions, 
				   epsilon, normTypeAux, Couplings_);
	norms_.Push_back(1.0);
      }
    }
    
    // Initialize each PDEs coupling terms
    for ( UInt i = 0; i < Couplings_.GetSize(); i++ ) {
      Couplings_[i]->GetPDE()->InitCoupling(Couplings_[i]); 
    }

#ifdef DEBUG    
    // write coupling data in .debug-file
    WriteCouplingInfo(*debug);
#endif

    // create solve step object
    solveStep_ = new IterSolveStep( *this );
  }


  void IterCoupledPDE::DefineSolvingPDEs(StdVector<StdPDE*> & pdes)
  {
    ENTER_FCN( "IterCoupledPDE::DefineSolvingPDEs" );

    solvePDE_.Init(FALSE);
    Boolean nameFound;
    std::string errMsg;

    for (UInt i=0; i<pdes.GetSize(); i++){
      nameFound = FALSE;
      for (UInt j=0; j<PDEs_.GetSize(); j++)
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



  // void IterCoupledPDE::SolveStepStatic(const UInt kstep, const Double aTime,
// 				       const Boolean updatesysmat ) {

//     ENTER_FCN ( "entering  IterCoupledPDE::SolveStepStatic" );
  
//     CFSVector *val, *oldVal;
//     UInt iter = 0;
//     UInt counter = 0;
//     Boolean normsReached = FALSE;
//     std::string quantityConv;


//     while (iter < maxiter_ &&  (! normsReached))
//       {
// 	if (nonLinLogging_)
// 	  {
// 	Info->PrintF(pdename_,"\n"); 
// 	Info->PrintF(pdename_, " COUPLED ITERATION %i =================================\n", 
// 		     iter+1);
// 	  }
	
// 	counter = 0;
// 	normsReached = TRUE;
      
// 	for (UInt i=0; i<PDEs_.GetSize(); i++)
// 	  {
// 	    if (nonLinLogging_)
// 	    Info->PrintF(pdename_, " Processing PDE %s\n", 
// 			 (PDEs_[i]->GetName()).c_str());

// 	    // Only solve current PDE, if the corresponding
// 	    // flag in 'solvePDE_' is set to TRUE
// 	    if (solvePDE_[i] == TRUE) {
// 	      PDEs_[i]->GetSolveStep()->PreStepStatic(kstep,aTime,updatesysmat);
// 	      PDEs_[i]->CalcInputCoupling();
// 	      PDEs_[i]->GetSolveStep()->SolveStepStatic(kstep,aTime,updatesysmat);
// 	      PDEs_[i]->GetSolveStep()->PostStepStatic(kstep,aTime);
// 	      PDEs_[i]->CalcOutputCoupling();
	      
// 	      // Calculate Norms
// 	      for (UInt k=0; k<Couplings_[i]->GetNumOutputCouplings(); k++)
// 		{
// 		  Couplings_[i]->GetOutputValues(k, val);
// 		  Couplings_[i]->GetOutputOldValues(k, oldVal);
// 		  norms_[counter] = CalcNorm(Couplings_[i]->GetOutputNormType(k), *val, *oldVal);

// 		  if (nonLinLogging_)
// 		    {
// 		      Enum2String(Couplings_[i]->GetOutputQuantity(k), quantityConv);
// 		      Info->PrintF(pdename_, " %s : Norm of %s = %g\n", 
// 				   (Couplings_[i]->GetPDEName()).c_str(),
// 				   quantityConv.c_str(), norms_[counter]);
// 		    }
		  
// 		  if (norms_[counter] > Couplings_[i]->GetOutputEpsilon(k) && 
// 		      Couplings_[i]->GetOutputNormType(k) != NO_NORM)
// 		    normsReached = FALSE;
		
// 		  //copy values of new solution to old one
// 		  *oldVal = *val;

// 		} // end of if 
// 		counter++;	      
// 	      }
// 	  }

// 	iter++;
	
// 	if(nonLinLogging_)
// 	Info->PrintF(pdename_, "\n"); 
//       }

//     // now we are converged and can compute any postprocessing-quantities
//     for (UInt i=0; i<PDEs_.GetSize(); i++)
//       PDEs_[i]->PostProcess(actlevel_);

//   }



  


//   void IterCoupledPDE::SolveStepTrans(const UInt kstep, const Double asteptime, const UInt level, 
// 				      const Boolean updatesysmat)
//   {
//     ENTER_FCN( "IterCoupledPDE::SolveStepTrans" );

//     Double  steptime  = asteptime;

//     UInt iter = 0;
//     Boolean normsReached = FALSE;
//     std::string quantityConv;

//     // In the beginning of each time step
//     // the coupling data has to be reseted
//     for (UInt i=0; i<PDEs_.GetSize(); i++)
//       PDEs_[i]->ResetCoupling();

//     while (iter < maxiter_ &&  (! normsReached))
//       {
// 	if (nonLinLogging_)
// 	  {
// 	    Info->PrintF(pdename_,"\n"); 
// 	    Info->PrintF(pdename_, " COUPLED ITERATION %i =================================\n", 
// 			 iter+1);
// 	  }

// 	UInt counter = 0;
// 	normsReached = TRUE;
      
// 	for (UInt i=0; i<PDEs_.GetSize(); i++)
// 	  {
// 	    if (nonLinLogging_)
// 	      Info->PrintF(pdename_, " Processing PDE %s\n", 
// 			   (PDEs_[i]->GetName()).c_str());

// 	    // Only solve current PDE, if the corresponding
// 	    // flag in 'solvePDE_' is set to TRUE
// 	    if (solvePDE_[i] == TRUE) {
	      
// 	      PDEs_[i]->GetSolveStep()->PreStepTrans(kstep, steptime, level, updatesysmat);
// 	      PDEs_[i]->CalcInputCoupling();
// 	      PDEs_[i]->GetSolveStep()->SolveStepTrans(kstep, steptime, level, updatesysmat);
// 	      PDEs_[i]->GetSolveStep()->PostStepTrans(kstep, steptime, level);
// 	      PDEs_[i]->CalcOutputCoupling();
	      
// 	      // Calculate Norms
// 	      for (UInt k=0; k<Couplings_[i]->GetNumOutputCouplings(); k++)
// 		{
// 		  CFSVector *val, *oldVal;
// 		  Couplings_[i]->GetOutputValues(k, val);
// 		  Couplings_[i]->GetOutputOldValues(k, oldVal);
// 		  norms_[counter] = CalcNorm(Couplings_[i]->GetOutputNormType(k), *val, *oldVal);

// 		  if (nonLinLogging_)
// 		    {
// 		      Enum2String(Couplings_[i]->GetOutputQuantity(k), quantityConv);
// 		      Info->PrintF(pdename_, " %s : Norm of %s = %g\n", 
// 				   (Couplings_[i]->GetPDEName()).c_str(),
// 				   quantityConv.c_str(), norms_[counter]);
// 		    }
// 		  if (norms_[counter] > Couplings_[i]->GetOutputEpsilon(k)) 
// 		    normsReached = FALSE;
		  
// 		  *oldVal = *val;
// 		} // end if
// 	      counter++;	      
// 	      }
// 	  }

// 	iter++;
// 	if (nonLinLogging_)
// 	  Info->PrintF(pdename_, "\n"); 
//       }
//   }

  

  void IterCoupledPDE::WriteResultsInFile(const UInt kstep,
					  const Double asteptime,
					  UInt stepOffset,
					  Double timeOffset)
{
  ENTER_FCN( "IterCoupledPDE::WriteResultsInFile" );

  for (UInt i=0; i<PDEs_.GetSize(); i++)
    PDEs_[i]->WriteResultsInFile(kstep, asteptime, stepOffset, timeOffset);
}


void IterCoupledPDE::WriteCouplingInfo(std::ostream &out)
{
  ENTER_FCN( "IterCoupledPDE::WriteCouplingInfo" );

  CFSVector *val;
  StdVector<UInt> * nodes;
  StdVector<std::string> couplingRegions;
  std::string couplingTypeAux, quantityAux, regionTypeAux, normTypeAux;

  // write information in .debug-file
  out << "=======================" << std::endl;
  out << " COUPLING INFORMATION  " << std::endl;
  out << "=======================" << std::endl;
  out << std::endl;

  for (UInt ipde=0; ipde<PDEs_.GetSize(); ipde++)
    {
      
      out << "Entering " << Couplings_[ipde]->GetPDE()->GetName() 
          << ".InitCoupling" << std::endl;
      out << "=====================================" << std::endl;
      
      
      // Show InputCouplings
      for (UInt i=0; i<Couplings_[ipde]->GetNumInputCouplings(); i++)
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
	  for (UInt j=0; j<couplingRegions.GetSize()-1; j++)
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
      for (UInt i=0; i<Couplings_[ipde]->GetNumOutputCouplings(); i++)
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
	      for (UInt j=0; j<couplingRegions.GetSize()-1; j++)
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


// Double IterCoupledPDE::CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval)
// {
//   ENTER_FCN( "IterCoupledPDE::CalcNorm" );

//   // ATTENTION: Currently only working with Double-values
//   // will be changed as soon as dynamic type information
//   // is available

//   Vector<Double> delta;
 
//   Double norm, valNorm2;
  

//   Vector<Double> & val_vec =\
//     dynamic_cast<Vector<Double>& >(val);

//   Vector<Double> & oldval_vec =\
//     dynamic_cast<Vector<Double>& >(oldval);
  
//   delta = val_vec - oldval_vec;

//   switch (normtype)
//     {
//     case NO_NORM:
//       return 0;
//       break;
      
//     case L2ABS:
//       norm = delta.NormL2();
//       break;

//     case L2REL:
//       valNorm2 =  val_vec.NormL2();
//       if (valNorm2 > 0)
// 	norm = delta.NormL2() / valNorm2;
//       else
// 	norm = delta.NormL2();

//       break;
//     }

//   return norm;
// }


void IterCoupledPDE::WriteGeneralPDEdefines()
{
  ENTER_FCN( "IterCoupledPDE::WriteGeneralPDEdefines" );
  
  for (UInt i=0; i<PDEs_.GetSize(); i++)
    PDEs_[i]->WriteGeneralPDEdefines();
}

BaseSolveStep * IterCoupledPDE::GetSolveStep() {
  ENTER_FCN( "IterCoupledPDE::GetSolveStep" );
  return solveStep_;
}

// ======================================================
// POSTPROC SECTION
// ======================================================

// Do Postprocessing as descriped in conf file
void IterCoupledPDE::PostProcess() 
{
  ENTER_FCN( "IterCoupledPDE::PostProcess" );

    for (UInt i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->PostProcess();
}



} // end of namespace
