#include "pdecoupling.hh"

#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Domain/bcs.hh"
#include "list"
#include "PDE/basePDE.hh"

namespace CoupledField
{


PDECoupling::CouplingInterface::CouplingInterface()
{
  
  level = 0;
  dof = 0;
  numNodes = 0;    
  numElems = 0;
  epsilon = 0;
  values = NULL;
  oldValues = NULL;
  }

PDECoupling::CouplingInterface::~CouplingInterface()
{
  if(values)
    delete values;
  if(oldValues)
    delete oldValues;
}


PDECoupling::PDECoupling(Grid * aptgrid, BCs * aptBCs)
{
  ENTER_FCN("PDECoupling::PDECoupling")
  
  ptGrid_ = aptgrid;
  ptBCs_ = aptBCs;
  
  defaultEpsilon = 1e-5;
  defaultNormType = L2REL;
  
}

PDECoupling::~PDECoupling()
{
  ENTER_FCN("PDECoupling::~PDECoupling")

  for (Integer i=0; i<outputInterfaces_.GetSize(); i++)
    if (outputInterfaces_[i]) delete outputInterfaces_[i];
  
  for (Integer i=0; i<inputInterfaces_.GetSize(); i++)
    if (inputInterfaces_[i]) delete inputInterfaces_[i];
  
}

void PDECoupling::RegisterInput(CouplingInputType InType, SolutionType Quantity)

{
  ENTER_FCN("PDECoupling::RegisterInput")
  
  inputTypes_.Push_back(InType);
  inputQuantities_.Push_back(Quantity);  
  inputInterfaces_.Push_back(0);
}


void PDECoupling::AddInput(SolutionType quantity, 
			   StdVector<std::string> &region, 
			   CouplingRegionType regionType,
			   Integer level,
			   Double epsilon,
			   NormType normtype,
			   StdVector<PDECoupling*> & couplings)
{
  ENTER_FCN("PDECoupling::AddInput");
  
  
  CouplingInterface *myInterface = 0; 
  Integer myNum = -1;
  std::string quantityConv;
  
  // search matching  Quantity
  for (Integer i=0; i<inputQuantities_.GetSize(); i++)
    {
      if (inputQuantities_[i] == quantity)
	myNum = i; 
    }
  
  if (myNum == -1)
    {
      Enum2String(quantity, quantityConv);
      std::string ErrMsg = "Quantity \'" + quantityConv + "\' not registered for PDE \'" + myPDE_->GetName()+ "\'";
      Error(ErrMsg.c_str(),__FILE__,__LINE__);
    }
  
  
  // Add coupling input as output to pde, which can calculate specified quantity
  CouplingOutputType myOutputType = Input2OutputType(inputTypes_[myNum]);
  
  Integer i=0;
  while (myInterface == 0 && i<couplings.GetSize())
    {
      myInterface = couplings[i++]->AddOutput(myOutputType, quantity, region, regionType, level);
    }
  
  // If no pde has the specified quantity as output
  if (myInterface == NULL)
    {
      Enum2String(quantity, quantityConv);
      std::string ErrMsg = "Qantity \'" + quantityConv 
	+ "\' can not be calculated with current set of PDEs";
      Error(ErrMsg.c_str(),__FILE__,__LINE__);
    }
  
  // i now contains the number of
  // the coupling object, which shares the data
  // of this coupling object
  i--;  // because i has been incremented before;
  BasePDE * oppositePDE = couplings[i]->myPDE_;
  
  
  // Set dof according to myPDE's dof
  // if the inputType is COORD (= moving geometry)
  // then the number of dofs is equal to the
  // dimension of the grid, otherwise the number dofs
  // must be equal to that of my PDE
  if (inputTypes_[myNum] == COORD)
    myInterface->dof = myPDE_->dim_;
  else
    myInterface->dof = myPDE_->dofspernode_;
  
  // Initialize the values and oldValues arrays
  //myInterface->values->SetDof(myPDE_->dofspernode_);
  //myInterface->values->Init(0.0);
  //myInterface->oldValues->SetDof(myPDE_->dofspernode_);
  //myInterface->oldValues->Init(0.0);
  
  // set normtype and epsilon
  //myInterface->epsilon = defaultEpsilon;
  //conf->ifget(inputQuantities_[myNum], myInterface->epsilon, "coupling", "tolerance"); 

  myInterface->epsilon = epsilon;
  myInterface->normtype = normtype;

  // std::string normtype;
//   myInterface->normtype = defaultNormType;
  
  
  // if (conf->ifget(inputQuantities_[myNum],normtype,"coupling","normtype") == TRUE)
//     {
      
//       if (normtype == "L2rel")
// 	myInterface->normtype = L2REL;
//       else if (normtype == "L2abs")
// 	myInterface->normtype = L2ABS;
//       else 
// 	{
// 	  std::string errMsg = "Normtype \'" + normtype + "\' not known!";
// 	  Error(errMsg.c_str(),__FILE__,__LINE__);
// 	} 
//     }
  
  if (myInterface->elements.GetSize() != 0)
    {
      // Set the material of the interface.
      // Since the Inputcoupling values are computed by another PDE
      // as Output values, it might need the material paramter
      // for some integrators.
      
      // 1. Step: get the neighbouring elements 
      
      const StdVector<Elem*> * interfaceElems = &(myInterface->elements);
      //StdVector<Elem*>  actSubdomain;
      //StdVector<Elem*> possibleNeighbours;
      
      
      
      // for (Integer iSd=0; iSd < myPDE_->subdoms_.GetSize(); iSd++)
// 	{GetVolNeighboursForSurf
// 	  ptGrid_->GetElemSD(actSubdomain, myPDE_->subdoms_[iSd], level);
// 	  for (Integer j=0; j<actSubdomain.GetSize(); j++)
// 	    possibleNeighbours.Push_back(actSubdomain[j]);
// 	}
      
      ptGrid_->GetVolNeighboursForSurf(*interfaceElems,  myPDE_->subdoms_, 
				       myInterface->neighbours, level);
      if (!myInterface->neighbours.GetSize())
	Error("No neighbours for element coupling found!",  __FILE__,__LINE__);
      
      
      
      
      // 2. Step: For each interface element, set the
      //          material parameter according to its neighbour
      
      for (Integer i=0; i<myInterface->neighbours.GetSize(); i++)
	{
	  Boolean subdomFound = FALSE;
	  Integer subDomNr = 0;
	  
	  for (subDomNr=0; subDomNr<myPDE_->subdoms_.GetSize(); subDomNr++)
	    if (myPDE_->subdoms_[subDomNr] == myInterface->neighbours[i]->namesd)
	      {
		subdomFound = TRUE;
		break;
	      }
	  
	  
	  if (subdomFound == FALSE)
	    Error("Subdomain name of neighbouring elements was not found",__FILE__,__LINE__);
	  
	  myInterface->materials[i] = &(myPDE_->materialData_[subDomNr]);	    
	}
      
      
      
      // Set the material of the neighbours elements of the "opposite" PDE
      // 1. Step: get the neighbouring elements
//       possibleNeighbours.Clear();
      
      
//       for (Integer iSd=0; iSd < oppositePDE->subdoms_.GetSize(); iSd++)
// 	{
// 	  ptGrid_->GetElemSD(actSubdomain, oppositePDE->subdoms_[iSd], level);
// 	  for (Integer j=0; j<actSubdomain.GetSize(); j++)
// 	    possibleNeighbours.Push_back(actSubdomain[j]);
// 	}
      
      ptGrid_->GetVolNeighboursForSurf(*interfaceElems, oppositePDE->subdoms_, 
				       myInterface->oppositePdeNeighbours, level);
      if (!myInterface->oppositePdeNeighbours.GetSize())
	Error("No opposite neighbours for element coupling found!",  __FILE__,__LINE__);
      
      for (Integer i=0; i<myInterface->oppositePdeNeighbours.GetSize(); i++)
	{
	  Boolean subdomFound = FALSE;
	  Integer subDomNr = 0;
	  
	  for (subDomNr=0; subDomNr<oppositePDE->subdoms_.GetSize(); subDomNr++)
	    if (oppositePDE->subdoms_[subDomNr] == myInterface->oppositePdeNeighbours[i]->namesd)
	      {
		subdomFound = TRUE;
		break;
	      }
	  
	  
	  
	  if (subdomFound == FALSE)
	    Error("Subdomain name of neighbouring elements was not found",__FILE__,__LINE__);
	  
	  myInterface->oppositePdeMaterials[i] = &(oppositePDE->materialData_[subDomNr]);
	}
      
      } // end if
  
  inputInterfaces_[myNum] = myInterface;
  
  
}

PDECoupling::CouplingInterface* PDECoupling::AddOutput(CouplingOutputType outputType, 
						       SolutionType quantity, 
						       StdVector<std::string> & regions,
						       CouplingRegionType regionType,
						       Integer level)
{
  ENTER_FCN("PDECoupling::AddOutput");
  
  // search if output exists already
  for (Integer i=0; i<outputTypes_.GetSize(); i++)
    if (outputTypes_[i] == outputType 
	&& outputQuantities_[i] == quantity
	&& outputInterfaces_[i]->regions == regions)
      return outputInterfaces_[i];
  
  
  // if not find out if PDE can calculate desired output
  if (myPDE_->HasOutput(quantity) == FALSE)
    return 0;
  
  //  std::cerr << "\n found output " << quantity << " in PDE " << myPDE_->GetName() << std::endl;
  //  std::cerr << "Type: " << regionType << std::endl;

  // create new Coupling Output
  outputTypes_.Push_back(outputType);
  outputQuantities_.Push_back(quantity);

  CouplingInterface *myInterface = new CouplingInterface;  
  myInterface->regions = regions;
  myInterface->regionType = regionType;
  myInterface->level = level;

  std::list<Integer> nodesConverted;
  StdVector<Elem*> SD;
  std::list<Integer>::iterator it;
  Integer inode = 0;
  Integer numNodes = 0;
  StdVector<Elem*> auxElems;

  // Get Elements/nodes of coupling region
  SD.Clear();
  switch (regionType)
    {
    case REGION:

      for (Integer iSD=0; iSD<regions.GetSize(); iSD++)
	{
	  ptGrid_->GetElemSD(auxElems, regions[iSD], level);
	  for (Integer iElem=0; iElem<auxElems.GetSize(); iElem++)
	    SD.Push_back(auxElems[iElem]);
	}


      ptGrid_->CalcNumberOfNodesInPatch(SD, myInterface->nodes);
      myInterface->numNodes = myInterface->nodes.GetSize();
      //myInterface->values.resize(myInterface->nodes.GetSize());
      //myInterface->oldValues.resize(myInterface->nodes.GetSize());
      //  myInterface->values->SetNumNodes(myInterface->nodes.GetSize());
      //       myInterface->oldValues->SetNumNodes(myInterface->nodes.GetSize());
	  
      
      break;

    case NODES:
      

      // count complete number of nodes
      for (Integer i=0; i<regions.GetSize(); i++)
	numNodes+= ptBCs_->GetNumNodesLevel(regions[i], level);

      myInterface->nodes.Resize(numNodes);
      
      inode =0;

      // get for each nodeslist all nodes
      for (Integer i=0; i<regions.GetSize(); i++)
	{
	  nodesConverted = ptBCs_->GetNodesLevel(regions[i], level);
	  
	  it = nodesConverted.begin();
	  
	  for (it=nodesConverted.begin(); it != nodesConverted.end(); it++, inode++)
	    myInterface->nodes[inode] = *it;
	}

      myInterface->numNodes = myInterface->nodes.GetSize();
      //myInterface->values.resize(myInterface->nodes.GetSize());
      //myInterface->oldValues.resize(myInterface->nodes.GetSize());
      // myInterface->values->SetNumNodes(myInterface->nodes.GetSize());
//       myInterface->oldValues->SetNumNodes(myInterface->nodes.GetSize());
      break;

    case SURFACE:

      for (Integer iSD=0; iSD<regions.GetSize(); iSD++)
	{
	  if (ptGrid_->GetDim() == 2)
	    auxElems = ptBCs_->getEdgesBC(regions[iSD], level);
	  else
	    auxElems = ptBCs_->getFacesBC(regions[iSD], level);
	  for (Integer iElem=0; iElem<auxElems.GetSize(); iElem++)
	    SD.Push_back(auxElems[iElem]);
	}

      
      // Get the nodes from BCs
      //myInterface->elements = ptBCs_->getEdgesBC(region, level);
      myInterface->elements = SD;

      // Get the nodes contained in the list of elements
      ptGrid_ -> CalcNumberOfNodesInPatch(myInterface->elements, myInterface->nodes);

      myInterface->numNodes = myInterface->nodes.GetSize();
      myInterface->numElems = myInterface->elements.GetSize();
      //myInterface->oldValues.resize(myInterface->nodes.GetSize());   
      //myInterface->values.resize(myInterface->nodes.GetSize());      
      // myInterface->values->SetNumNodes(myInterface->nodes.GetSize());
      // myInterface->oldValues->SetNumNodes(myInterface->nodes.GetSize());
      myInterface->materials.Resize(myInterface->elements.GetSize());
      myInterface->oppositePdeMaterials.Resize(myInterface->elements.GetSize());
	    

      break;

    // case ELEMS2D:
//       for (Integer iSD=0; iSD<regions.GetSize(); iSD++)
// 	{
// 	  auxElems = ptBCs_->getFacesBC(regions[iSD], level);
// 	  for (Integer iElem=0; iElem<auxElems.GetSize(); iElem++)
// 	    SD.Push_back(auxElems[iElem]);
// 	}
      
//       myInterface->elements = SD;
//       ptGrid_ -> CalcNumberOfNodesInPatch(myInterface->elements, myInterface->nodes);
      
//       myInterface->numNodes = myInterface->nodes.GetSize();
//       myInterface->numElems = myInterface->elements.GetSize();
//       //myInterface->values.resize(myInterface->nodes.GetSize());
//       //myInterface->oldValues.resize(myInterface->nodes.GetSize());
//       // myInterface->values->SetNumNodes(myInterface->nodes.GetSize());
// //       myInterface->oldValues->SetNumNodes(myInterface->nodes.GetSize());
//       myInterface->materials.Resize(myInterface->elements.GetSize());
//       break;
    }

  outputInterfaces_.Push_back(myInterface);
  
  return myInterface;
}




void PDECoupling::SetPDE(BasePDE * aPDE)
{
  ENTER_FCN("PDECoupling::SetPDE");
  myPDE_ = aPDE;
}

void PDECoupling::SetOutputNumNodes(Integer i, Integer nnodes)
{
  ENTER_FCN("PDECoupling::SetOutputNumNodes");
  
  outputInterfaces_[i]->numNodes = nnodes;
    
}  

void PDECoupling::SetOutputNumElems(Integer i, Integer nelems)
{
  ENTER_FCN("PDECoupling::SetOutputNumElems");
  outputInterfaces_[i]->numElems = nelems;
  
}  


void PDECoupling::SetOutputDof(Integer i, ShortInt dof)
{
  ENTER_FCN("PDECoupling::SetOutputDof");
  outputInterfaces_[i]->dof = dof;
}

std::string PDECoupling::GetPDEName()
{
  ENTER_FCN("PDECoupling::GetPDEName");
  return myPDE_->GetName();
}


Integer PDECoupling::GetNumInputCouplings()
{
  ENTER_FCN("PDECoupling::GetnumInputCouplings");
  return inputInterfaces_.GetSize();
}


Integer PDECoupling::GetNumOutputCouplings()
{
  ENTER_FCN("PDECoupling::GetnumOutputCouplings");
  return outputQuantities_.GetSize();
}


void PDECoupling::CreateCouplingVector(Integer i,
				       Boolean isComplex)
{
  ENTER_FCN("PDECoupling::CreateCouplingVector");

  Integer numNodes = outputInterfaces_[i]->numNodes;
  Integer dof =  outputInterfaces_[i]->dof;

  if (outputInterfaces_[i]->values != NULL || \
      outputInterfaces_[i]->oldValues != NULL)
    Error("PDECoupling::CreateCouplingVector: This function may only\
           be called for Output-Coupling interfaces!",__FILE__,__LINE__);

  if (isComplex)
    {
      outputInterfaces_[i]->values = new Vector<Complex>;
      outputInterfaces_[i]->oldValues = new Vector<Complex>;
    } 
  else
    {
      outputInterfaces_[i]->values = new Vector<Double>;
      outputInterfaces_[i]->oldValues = new Vector<Double>;
    }
  
  // resize vector with coupling values
  outputInterfaces_[i]->values->Resize(dof*numNodes);
  
  // resize vector with old coupling values
  outputInterfaces_[i]->oldValues->Resize(dof*numNodes);

 }


CouplingOutputType PDECoupling::Input2OutputType(CouplingInputType inputType)
{
  ENTER_FCN("PDECoupling::Input2OutputType");

  switch (inputType)
    {
    case COORD:
      // Coordinate Coupling -> node values needed
      return NODE;
      break;
      
    case RHS:
      // RHS coupling -> node values needed
      return NODE;
      break;
      
    case ID_BC:
      // Inhomogeneous Dirichlet BC -> node values needed
      return NODE;
      break;
      
    case MAT:
      // Material parameter coupling -> element values needed
      return ELEM;
      break;
    }
  
}



//   std::ostream& operator << ( std::ostream & out , const CouplingInterface & inter)
//   {
//     out << myendl 
// 	<< "Interface: coupling region: " << inter.region << myendl
// 	<< "           nr. coupling nodes: " <<  inter.nodes.GetSize() << myendl
// 	<< "           nr. coupling elements: " <<  inter.elements.GetSize() << myendl
// 	<< "           name of 1. coupling element: " <<  inter.elements[0].namesd << myendl;
// 	<< "           nr. neighbours elements: " <<  inter.neighbours.GetSize() << myendl
// 	<< "           name of 1. neighbours element: " <<  inter.neighbours[0].namesd << myendl;

//   } 

} // end of namespace
