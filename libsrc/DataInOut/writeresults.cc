#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <list>

#include "writeresults.hh"
#include "ParamHandling/ConfFile.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"


namespace CoupledField {

  WriteResults::WriteResults(const Char * const filename,
			     FileType * const aInFile)
  {

    ENTER_FCN( "WriteResults::WriteResults" );

    namefile_ = filename;
   //  namefile_ = new Char[strlen(filename)+1+MAXPOSTFIX];
//     strcpy(namefile_,filename);

  ascii_=TRUE;

  pt2Inputfile_ = aInFile;

  InitHistoryFiles();
}


WriteResults::~WriteResults()
{
  ENTER_FCN( "WriteResults::~WriteResults" );
  
  for (Integer i=0; i<historyFiles_.GetSize(); i++)
    for (Integer j=0; j<historyFiles_[i].GetSize(); j++)
      if (historyFiles_[i][j])
	delete historyFiles_[i][j];
  
}



void WriteResults::InitHistoryFiles()
{
  ENTER_FCN( "WriteResults::InitHistoryFiles" );

 
 StdVector<Integer> nodesTmp;

#ifndef XMLPARAMS
 conf->getlist(nodesTmp,"history_node");
 Error("History files are not supported anymore for .conf-file format!");
#else
 
 StdVector<std::string> keyVec, attrVec, valVec;
 StdVector<std::string> nodeVec, quantVec;
 SolutionType solType;
 
 attrVec = "", "";
 valVec = "", "";

 keyVec = "storeResults", "nodeHistory", "saveNodes";
 params->GetList(keyVec, attrVec, valVec, nodeVec);

 keyVec = "storeResults", "nodeHistory", "type";
 params->GetList(keyVec, attrVec, valVec, quantVec);

//  std::cerr << "nodeVec = " << std::endl << nodeVec << std::endl;
//  std::cerr << "quantVec = " << std::endl << quantVec << std::endl;


 NeedHistory_ = FALSE;
 if (nodeVec.GetSize() > 0)
   NeedHistory_ = TRUE;

 if (NeedHistory_ == FALSE)
   return;

 // Loop over all read in quantities and nodes
 for (Integer iQuant=0; iQuant<quantVec.GetSize(); iQuant++) {
   
   Integer quantityFound = -1;

   // Look in all histQuantites, if quantity was already read in
   for (Integer j=0; j<histQuantities_.GetSize(); j++) {
     String2Enum(quantVec[iQuant], solType);
     if (histQuantities_[j] == solType) {
       quantityFound = j;
       break;
     }
   }
  
   
   if (quantityFound == -1) {
     String2Enum(quantVec[iQuant], solType);
     histQuantities_.Push_back(solType);
     histNodesPerPDE_.Push_back(StdVector<Integer>());
     quantityFound = histNodesPerPDE_.GetSize() -1;
   }
   
   // Add new nodes
   StdVector<Integer> tempNodes;
   pt2Inputfile_->ReadSaveNodes(tempNodes, nodeVec[iQuant]);
   for (Integer i=0; i<tempNodes.GetSize(); i++)
     histNodesPerPDE_[quantityFound].Push_back(tempNodes[i]);

//    std::list<Integer> tempNodes;
//    tempNodes= ptBCs_->GetNodesLevel(nodeVec[iQuant]);
//    for (std::list<Integer>::const_iterator p=tempNodes.begin();
// 	p!=tempNodes.end(); p++)
//      histNodesPerPDE_[quantityFound].Push_back(*p);

    } // iQuant

//  std::cerr << "histQuantities_ " <<histQuantities_ << std::endl; 
//  std::cerr << "length of nodes" <<  histNodesPerPDE_.GetSize() << std::endl;
//  std::cerr << "nodes[0] = " << histNodesPerPDE_[0] << std::endl;
 // Create history directory
 std::string S="mkdir -p history";
 system(S.c_str());
 
 std::string namePrefix="history/" + namefile_ + "-";
 std::string namePostfix = ".hist";
 std::string totalName;
 
 // Resize number of history files
 historyFiles_.Resize(histQuantities_.GetSize());
 
 // Iterate over all quantities
 for (Integer iQuant=0; iQuant<histQuantities_.GetSize(); iQuant++) {
   
   std::string comment = "Nodehistory for quantity '";
   comment += histQuantities_[iQuant] + "'";
   Info->PrintVec(comment.c_str(), histNodesPerPDE_[iQuant]);

   historyFiles_[iQuant].Resize(histNodesPerPDE_[iQuant].GetSize());
   for (Integer iNode=0; iNode<histNodesPerPDE_[iQuant].GetSize(); iNode++) {

     // Create total filename
     std::string quantString;
     Enum2String(histQuantities_[iQuant], quantString);
     totalName = namePrefix + quantString;
     totalName += "-node-";
     totalName += Info->GenStr(histNodesPerPDE_[iQuant][iNode]);
     totalName += namePostfix;

     historyFiles_[iQuant][iNode] = NULL;
//      std::cerr << "creating hist-file" << totalName << std::endl;
     historyFiles_[iQuant][iNode] = new std::ofstream(totalName.c_str());
     
     if (!historyFiles_[iQuant][iNode]) {
       std::string errMsg = "InitHistory: Can not open history file '";
       errMsg += totalName + "'";
       Error (errMsg.c_str(), __FILE__, __LINE__);
     }
     
   }
 }
 
#endif
}

void WriteResults::WriteSolMatrix(Grid * ptgrid, const Integer level, const Vector<Double> sol, 
				  const std::string matFileName, const Integer nrDofs)
{
  ENTER_FCN( "WriteResults::WriteSolMatrix" );

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

void WriteResults::WriteNodeHistoryTransient(const NodeStoreSol<Double>& data, 
					     const Integer step, 
					     const Double time)
{
  ENTER_FCN( "WriteResults::WriteNodeHistoryTransien" );

  std::ofstream * myHist;
  SolutionType actSolType;
  StdVector<SolutionType> solTypes;
  Integer iQuant, actDof;
  std::string quantity;
  Double val;


  data.GetSolutionTypes(solTypes);
  
  // Iterate over all solutiontypes
  for (Integer iSol=0; iSol<solTypes.GetSize(); iSol++) {
   
    actDof = data.GetDof(solTypes[iSol]);
    
    // Find the according quantitiy
    iQuant = -1;
    for (Integer i=0; i<histQuantities_.GetSize(); i++) {
     
 std::cerr << "Quantity = " << histQuantities_[i] << std::endl;
      if (histQuantities_[i] == solTypes[iSol]){
	iQuant = i;
	break;
      }
    }
    

    if (iQuant == -1){
      Enum2String(solTypes[iSol], quantity);
      std::string errMsg = "Quantity '";
      errMsg += quantity;
      errMsg += "' for history not found!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
    
    // Iterate over all history nodes
    for (Integer iNode=0; 
	 iNode<histNodesPerPDE_[iQuant].GetSize(); 
	 iNode++) {
      myHist = historyFiles_[iQuant][iNode];
      (*myHist) << time;
      
      // Iterate over all dofs
      for (Integer iDof=0; iDof<actDof; iDof++) {
// 	std::cerr << "Writing out Sol " << iSol << ", Node " << histNodesPerPDE_[iQuant][iNode];
// 	std::cerr << ", and dof " << iDof << std::endl;
	data.Get(solTypes[iSol],histNodesPerPDE_[iQuant][iNode]-1, iDof, val);
	(*myHist) << "  " << val;
      } // iDof
      (*myHist) << std::endl;
    }  // iNode
  } // iSol
}

void WriteResults::WriteNodeHistoryHarmonic(const NodeStoreSol<Complex>& data, 
					    const Integer step,
					    const Double frequency,
					    const ComplexFormat format)
{
  ENTER_FCN( "WriteResults::WriteNodeHistoryHarmonic" );
}

} // end of namespace
