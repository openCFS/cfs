#include "pdecoupling.hh"

#include <Domain/elem.hh>
#include <Domain/grid.hh>
#include <Domain/bcs.hh>
#include <list>
#include <Utils/storesol.hh>

#include <PDE/basePDE.hh>

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

  for (Integer i=0; i<outputInterfaces_.size(); i++)
    if (outputInterfaces_[i]) delete outputInterfaces_[i];
  
  for (Integer i=0; i<inputInterfaces_.size(); i++)
    if (inputInterfaces_[i]) delete inputInterfaces_[i];
  
}

void PDECoupling::RegisterInput(CouplingInputType InType, std::string Quantity)
{
  ENTER_FCN("PDECoupling::RegisterInput")
  
  inputTypes_.push_back(InType);
  inputQuantities_.push_back(Quantity);  
  inputInterfaces_.push_back(0);
}


void PDECoupling::AddInput(std::string quantity, 
			   std::string region, 
			   CouplingRegionType regionType,
			   Integer level,
			   std::vector<PDECoupling*> & couplings)
{
  ENTER_FCN("PDECoupling::AddInput");
  
  
  CouplingInterface *myInterface = 0; 
  Integer myNum = -1;
  
  // search matching  Quantity
  for (Integer i=0; i<inputQuantities_.size(); i++)
    {
      if (inputQuantities_[i] == quantity)
	myNum = i; 
    }
  
  if (myNum == -1)
    {
      std::string ErrMsg = "Quantity \'" + quantity + "\' not registered for PDE \'" + myPDE_->GetName()+ "\'";
      Error(ErrMsg.c_str(),__FILE__,__LINE__);
    }
  
  
  // Add coupling input as output to pde, which can calculate specified quantity
  CouplingOutputType myOutputType = Input2OutputType(inputTypes_[myNum]);
  
  Integer i=0;
  while (myInterface == 0 && i<couplings.size())
    {
      myInterface = couplings[i++]->AddOutput(myOutputType, quantity, region, regionType, level);
    }
  
  // If no pde has the specified quantity as output
  if (myInterface == NULL)
    {
      std::string ErrMsg = "Qantity \'" + quantity 
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
  myInterface->epsilon = defaultEpsilon;
  conf->ifget(inputQuantities_[myNum], myInterface->epsilon, "coupling", "tolerance"); 
  
  std::string normtype;
  myInterface->normtype = defaultNormType;
  
  
  if (conf->ifget(inputQuantities_[myNum],normtype,"coupling","normtype") == TRUE)
    {
      
      if (normtype == "L2rel")
	myInterface->normtype = L2REL;
      else if (normtype == "L2abs")
	myInterface->normtype = L2ABS;
      else 
	{
	  std::string errMsg = "Normtype \'" + normtype + "\' not known!";
	  Error(errMsg.c_str(),__FILE__,__LINE__);
	} 
    }
  
  if (myInterface->elements.size() != 0)
    {
      // Set the material of the interface.
      // Since the Inputcoupling values are computed by another PDE
      // as Output values, it might need the material paramter
      // for some integrators.
      
      // 1. Step: get the neighbouring elements 
      
      const std::vector<Elem*> * interfaceElems = &(myInterface->elements);
      std::vector<Elem*>  actSubdomain, possibleNeighbours;
      
      
      
      for (Integer iSd=0; iSd < myPDE_->subdoms_.size(); iSd++)
	{
	  ptGrid_->GetElemSD(actSubdomain, myPDE_->subdoms_[iSd], level);
	  for (Integer j=0; j<actSubdomain.size(); j++)
	    possibleNeighbours.push_back(actSubdomain[j]);
	}
      
      ptGrid_->DefineBelonging4Elems(*interfaceElems, possibleNeighbours, myInterface->neighbours);
      if (!myInterface->neighbours.size())
	Error("No neighbours for element coupling found!",  __FILE__,__LINE__);
      
      
      
      
      // 2. Step: For each interface element, set the
      //          material parameter according to its neighbour
      
      for (Integer i=0; i<myInterface->neighbours.size(); i++)
	{
	  Boolean subdomFound = FALSE;
	  Integer subDomNr = 0;
	  
	  for (subDomNr=0; subDomNr<myPDE_->subdoms_.size(); subDomNr++)
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
      possibleNeighbours.clear();
      
      
      for (Integer iSd=0; iSd < oppositePDE->subdoms_.size(); iSd++)
	{
	  ptGrid_->GetElemSD(actSubdomain, oppositePDE->subdoms_[iSd], level);
	  for (Integer j=0; j<actSubdomain.size(); j++)
	    possibleNeighbours.push_back(actSubdomain[j]);
	}
      
      ptGrid_->DefineBelonging4Elems(*interfaceElems, possibleNeighbours, myInterface->oppositePdeNeighbours);
      if (!myInterface->oppositePdeNeighbours.size())
	Error("No opposite neighbours for element coupling found!",  __FILE__,__LINE__);
      
      for (Integer i=0; i<myInterface->oppositePdeNeighbours.size(); i++)
	{
	  Boolean subdomFound = FALSE;
	  Integer subDomNr = 0;
	  
	  for (subDomNr=0; subDomNr<oppositePDE->subdoms_.size(); subDomNr++)
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
						       std::string quantity, 
						       std::string region,
						       CouplingRegionType regionType,
						       Integer level)
{
  ENTER_FCN("PDECoupling::AddOutput");
  
  // search if output exists already
  for (Integer i=0; i<outputTypes_.size(); i++)
    if (outputTypes_[i] == outputType 
	&& outputQuantities_[i] == quantity
	&& outputInterfaces_[i]->region == region)
      return outputInterfaces_[i];
  
  
  // if not find out if PDE can calculate desired output
  if (myPDE_->HasOutput(quantity) == FALSE)
    return 0;
  
  //std::cerr << "found output " << Quantity << " in PDE " << myPDE_->GetName() << std::endl;

  // create new Coupling Output
  outputTypes_.push_back(outputType);
  outputQuantities_.push_back(quantity);

  CouplingInterface *myInterface = new CouplingInterface;  
  myInterface->region = region;
  myInterface->regionType = regionType;
  myInterface->level = level;

  std::list<Integer> nodesConverted;
  std::vector<Elem*> SD;
  std::list<Integer>::iterator it;
  Integer inode = 0;


  // Get Elements/nodes of coupling region
  switch (regionType)
    {
    case SUBDOMAIN:
      ptGrid_->GetElemSD(SD, region, level);
      ptGrid_->CalcNumberOfNodesInPatch(SD, myInterface->nodes);
      myInterface->numNodes = myInterface->nodes.size();
      //myInterface->values.resize(myInterface->nodes.size());
      //myInterface->oldValues.resize(myInterface->nodes.size());
     //  myInterface->values->SetNumNodes(myInterface->nodes.size());
//       myInterface->oldValues->SetNumNodes(myInterface->nodes.size());
	    
      
      break;

    case NODES:

      nodesConverted = ptBCs_->GetNodesLevel(region, level);
      myInterface->nodes.resize(nodesConverted.size());
      it = nodesConverted.begin();
      
      inode = 0;
      for (it=nodesConverted.begin(); it != nodesConverted.end(); it++, inode++)
	myInterface->nodes[inode] = *it;

      myInterface->numNodes = myInterface->nodes.size();
      //myInterface->values.resize(myInterface->nodes.size());
      //myInterface->oldValues.resize(myInterface->nodes.size());
      // myInterface->values->SetNumNodes(myInterface->nodes.size());
//       myInterface->oldValues->SetNumNodes(myInterface->nodes.size());
      break;

    case ELEMS1D:

      // Get the nodes from BCs
      myInterface->elements = ptBCs_->getEdgesBC(region, level);

      // Get the nodes contained in the list of elements
      ptGrid_ -> CalcNumberOfNodesInPatch(myInterface->elements, myInterface->nodes);

            
      myInterface->numNodes = myInterface->nodes.size();
      myInterface->numElems = myInterface->elements.size();
      //myInterface->oldValues.resize(myInterface->nodes.size());   
      //myInterface->values.resize(myInterface->nodes.size());      
      // myInterface->values->SetNumNodes(myInterface->nodes.size());
      // myInterface->oldValues->SetNumNodes(myInterface->nodes.size());
      myInterface->materials.resize(myInterface->elements.size());
      myInterface->oppositePdeMaterials.resize(myInterface->elements.size());
	    

      break;

    case ELEMS2D:
      myInterface->elements = ptBCs_->getFacesBC(region, level);
      ptGrid_ -> CalcNumberOfNodesInPatch(myInterface->elements, myInterface->nodes);
      
      myInterface->numNodes = myInterface->nodes.size();
      myInterface->numElems = myInterface->elements.size();
      //myInterface->values.resize(myInterface->nodes.size());
      //myInterface->oldValues.resize(myInterface->nodes.size());
      // myInterface->values->SetNumNodes(myInterface->nodes.size());
//       myInterface->oldValues->SetNumNodes(myInterface->nodes.size());
      myInterface->materials.resize(myInterface->elements.size());
      break;
    }

  outputInterfaces_.push_back(myInterface);
  
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
  return inputInterfaces_.size();
}


Integer PDECoupling::GetNumOutputCouplings()
{
  ENTER_FCN("PDECoupling::GetnumOutputCouplings");
  return outputQuantities_.size();
}


void PDECoupling::CreateStoreSol(Integer i,
				 SolutionType solType, 
				 Boolean isComplex)
{
  ENTER_FCN("PDECoupling::CreateStoreSol");

  Integer numNodes = outputInterfaces_[i]->numNodes;
  Integer dof =  outputInterfaces_[i]->dof;

  if (outputInterfaces_[i]->values != NULL || \
      outputInterfaces_[i]->oldValues != NULL)
    Error("PDECoupling::CreateStoreSol: This function may only\
           be called for Output-Coupling interfaces!",__FILE__,__LINE__);

  if (isComplex)
    {
      outputInterfaces_[i]->values = new StoreSol<Complex>;
      outputInterfaces_[i]->oldValues = new StoreSol<Complex>;
    } 
  else
    {
      outputInterfaces_[i]->values = new StoreSol<Double>;
      outputInterfaces_[i]->oldValues = new StoreSol<Double>;
    }
  
  // initialize storesol-object values
  outputInterfaces_[i]->values->SetNumSolutions(1);
  outputInterfaces_[i]->values->SetNumNodes(numNodes);
  outputInterfaces_[i]->values->SetSolutionType(solType);
  outputInterfaces_[i]->values->SetNumDofs(dof);
  outputInterfaces_[i]->values->Init(0.0);
  
  // initialize storesol-object oldValues
  outputInterfaces_[i]->oldValues->SetNumSolutions(1);
  outputInterfaces_[i]->oldValues->SetNumNodes(numNodes);
  outputInterfaces_[i]->oldValues->SetSolutionType(solType);
  outputInterfaces_[i]->oldValues->SetNumDofs(dof);
  outputInterfaces_[i]->oldValues->Init(0.0);

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
// 	<< "           nr. coupling nodes: " <<  inter.nodes.size() << myendl
// 	<< "           nr. coupling elements: " <<  inter.elements.size() << myendl
// 	<< "           name of 1. coupling element: " <<  inter.elements[0].namesd << myendl;
// 	<< "           nr. neighbours elements: " <<  inter.neighbours.size() << myendl
// 	<< "           name of 1. neighbours element: " <<  inter.neighbours[0].namesd << myendl;

//   } 

} // end of namespace
