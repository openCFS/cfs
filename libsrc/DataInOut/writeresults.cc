#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <list>

#include "writeresults.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"


namespace CoupledField {

  // ===============
  //   Constructor
  // ===============
  WriteResults::WriteResults( const Char * const filename,
			      FileType * const aInFile ) {

    ENTER_FCN( "WriteResults::WriteResults" );

    namefile_ = filename;
    NeedHistory_ = FALSE;
    ascii_=TRUE;
    
    pt2Inputfile_ = aInFile;
  }


  // ======================
  //   Default Destructor
  // ======================
  WriteResults::~WriteResults() {
    ENTER_FCN( "WriteResults::~WriteResults" );
  
    for (Integer i=0; i<historyFiles_.GetSize(); i++) {
      for (Integer j=0; j<historyFiles_[i].GetSize(); j++) {
	if (historyFiles_[i][j]) {
	  delete historyFiles_[i][j];
	}
      }
    }
  }


  // ==========================
  //   Initialisation Routine
  // ==========================
  void WriteResults::InitHistoryFiles() {

    ENTER_FCN( "WriteResults::InitHistoryFiles" );
 
    StdVector<Integer> nodesTmp;

    // *****************************************************************
    //   Determine for which output quantities there are history nodes
    //   and what their names are
    // *****************************************************************
    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<std::string> nodeVec, quantVec;
    SolutionType solType;
 
    attrVec = "", "";
    valVec = "", "";

    keyVec = "storeResults", "nodeHistory", "saveNodes";
    params->GetList(keyVec, attrVec, valVec, nodeVec);

    keyVec = "storeResults", "nodeHistory", "type";
    params->GetList(keyVec, attrVec, valVec, quantVec);

    // Check if there are any history nodes at all
    // and set NeedHistory_ flag accordingly
    NeedHistory_ = FALSE;
    if (nodeVec.GetSize() > 0) {
      NeedHistory_ = TRUE;
    }

    // If there is aren't any history nodes,
    // then there is nothing to do
    if ( NeedHistory_ == FALSE ) {
      return;
    }

    // *****************************************************************
    //   Generate vectors holding information on the output quantities
    //   and nodes for which the user wants to see the history
    //
    //   histQuantities_ ....... holds identifiers for all quantities
    //   histNodesPerQuant_ .... node numbers of the history nodes for
    //                           each of the output quantities
    //
    // *****************************************************************

    // Make a sanity check
    if ( histQuantities_.GetSize() != 0 || histNodesPerQuant_.GetSize() != 0 ){
      Error( "Repeated call to WriteResults::InitHistoryFiles ?!?!",
	     __FILE__, __LINE__ );
    }

    // Define some loop counters and auxilliary variables
    Integer iQuant, i, j;
    Integer quantityFound;

    // Local data structure for combining identifiers of history nodes
    // with the respective node numbers.
    // For each output quantity we store a StdVector that stores triples
    // containing
    //  a) index of identifier in nodeVec
    //  b) start index of node numbers for identifier in histNodesPerQuant_
    //  c) stop index of node numbers for identifier in histNodesPerQuant_
    StdVector< StdVector<Integer> > histNodeNumIdentCoup;

    // In the vector of quantities there may be duplicates
    // these must be eliminated!
    for ( iQuant = 0; iQuant < quantVec.GetSize(); iQuant++ ) {

      // (re-)initialise
      quantityFound = -1;

      // Determine if current output quantity was already assigned
      for ( j = 0; j < histQuantities_.GetSize(); j++ ) {
	String2Enum( quantVec[iQuant], solType );
	if ( histQuantities_[j] == solType ) {
	  quantityFound = j;
	  break;
	}
      }

      // The quantity was not yet assigned, so do this now
      // and also make room for storing the curresponding
      // number of history nodes
      if ( quantityFound == -1 ) {
	String2Enum( quantVec[iQuant], solType );
	histQuantities_.Push_back(solType);
	histNodesPerQuant_.Push_back(StdVector<Integer>());
	histNodeNumIdentCoup.Push_back(StdVector<Integer>());
	quantityFound = histNodesPerQuant_.GetSize() - 1;
      }

      // Whenever an output quantity appears in quantVec it is
      // related to a saveNodes identifier. Add the corresponding
      // nodes to the histNodesPerQuant_ vector and the indices of
      // the identifier and node numbers to histNodeNumIdentCoup
      StdVector<Integer> tempNodes;
      pt2Inputfile_->ReadSaveNodes( tempNodes, nodeVec[iQuant] );
      histNodeNumIdentCoup[quantityFound].Push_back(iQuant);
      histNodeNumIdentCoup[quantityFound].
	Push_back(histNodesPerQuant_[quantityFound].GetSize());
      for ( i = 0; i < tempNodes.GetSize(); i++ ) {
	histNodesPerQuant_[quantityFound].Push_back(tempNodes[i]);
      }
      histNodeNumIdentCoup[quantityFound].
	Push_back(histNodesPerQuant_[quantityFound].GetSize()-1);
    }

    // ************************
    //   Prepare output files
    // ************************

    // Create history directory
    std::string S="mkdir -p history";
    system(S.c_str());
 
    std::string namePrefix="history/" + namefile_ + "-";
    std::string namePostfix = ".hist";
    std::string totalName;

    // Resize number of history files
    historyFiles_.Resize(histQuantities_.GetSize());

    // Iterate over all output quantities
    std::string quantityName;
    for ( iQuant = 0; iQuant < histQuantities_.GetSize(); iQuant++ ) {

      // Log information on this output quantity to info file.
      // This will also tell the user the relation between the
      // node numbers and the identifiers
      // std::string comment = "Identifiers for history nodes of quantity '";
      Enum2String( histQuantities_[iQuant], quantityName );

      Info->PrintF( "", "\nHistory nodes for quantity '%s':\n",
		    quantityName.c_str() );
      for ( i = 0; i < histNodeNumIdentCoup[iQuant].GetSize() / 3; i++ ) {

	Info->PrintF( "", " identifier = %s, node(s) =", 
		      nodeVec[histNodeNumIdentCoup[iQuant][i*3]].c_str() );

	for ( j  = histNodeNumIdentCoup[iQuant][i*3+1];
	      j <= histNodeNumIdentCoup[iQuant][i*3+2]; j++ ) {
	  Info->PrintF( "", " %d", histNodesPerQuant_[iQuant][j] );
	}
	Info->PrintF( "", "\n" );
      }

      // For each node there will be one history file
      historyFiles_[iQuant].Resize(histNodesPerQuant_[iQuant].GetSize());

      for ( Integer iNode = 0; iNode < histNodesPerQuant_[iQuant].GetSize();
	    iNode++ ) {

	// Create complete filename
	std::string quantString;
	Enum2String(histQuantities_[iQuant], quantString);
	totalName = namePrefix + quantString;
	totalName += "-node-";
	totalName += Info->GenStr(histNodesPerQuant_[iQuant][iNode]);
	totalName += namePostfix;

	// Open output file
	historyFiles_[iQuant][iNode] = NULL;
	historyFiles_[iQuant][iNode] = new std::ofstream(totalName.c_str());
	if ( !historyFiles_[iQuant][iNode] ) {
	  std::string errMsg = "InitHistory: Cannot open history file '";
	  errMsg += totalName + "'";
	  Error (errMsg.c_str(), __FILE__, __LINE__);
	}
      }
    }
 
  }


  void WriteResults::WriteSolMatrix(Grid * ptgrid, const Integer level,
				    const Vector<Double> sol, 
				    const std::string matFileName,
				    const Integer nrDofs) {

    ENTER_FCN( "WriteResults::WriteSolMatrix" );

    //get and write number of nodes on the level
    Integer numnodes=ptgrid->GetMaxnumnodes(level);
    Integer dim=ptgrid->GetDim();
    Integer i;

    std::ofstream * matrixOut = new std::ofstream(matFileName.c_str());
    // use scientific output, formatting is much better
    matrixOut->setf(std::ios::scientific, std::ios::floatfield);
  
    if (dim==2) {
      Point<2> point;
      
      // write x,y,z-coordinate
      for (i=0; i<numnodes; i++) {
	ptgrid->GetCoordinateNode(i,level,point);
	(*matrixOut) << " \t" << point[0] << " \t" << point[1] << " \t"
		     << 0 << " \t";

	for (Integer actDof =0; actDof < nrDofs; actDof++) {
	  (*matrixOut) << sol[i*nrDofs + actDof] << " \t";
	}
	(*matrixOut) << std::endl;
      }
    }
    else {
      Point<3> point;
      
      // write x,y,z-coordinate
      for (i=0; i<numnodes; i++) {
	ptgrid->GetCoordinateNode(i,level,point);
	(*matrixOut) << " \t" << point[0] << " \t" << point[1] << " \t"
		     << point[2] << " \t";

	for (Integer actDof =0; actDof < nrDofs; actDof++) {
	  (*matrixOut) << sol[i*nrDofs + actDof] << " \t";
	}
	(*matrixOut) << std::endl;
      } 
    }
  }


  void WriteResults::WriteNodeHistoryTransient(const NodeStoreSol<Double>&data,
					       const Integer step,
					       const Double time) {

    ENTER_FCN( "WriteResults::WriteNodeHistoryTransient" );

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
    
      // Find the related quantity
      iQuant = -1;
      for (Integer i=0; i<histQuantities_.GetSize(); i++) {

	if (histQuantities_[i] == solTypes[iSol]){
	  iQuant = i;
	  break;
	}
      }
      
      if ( iQuant == -1 ){

	// Report what quantities we have
	std::cerr << "Found the following quantities:\n";
	for ( Integer i = 0; i < histQuantities_.GetSize(); i++ ) {
	  std::cerr << " " << histQuantities_[i] << "\n";
	}

	Enum2String( solTypes[iSol], quantity );
	(*error) << "Quantity '" << quantity << "' for history not found!";
	Error( __FILE__, __LINE__ );
      }
    
      // Iterate over all history nodes
      for ( Integer iNode = 0; iNode < histNodesPerQuant_[iQuant].GetSize();
	    iNode++ ) {
	myHist = historyFiles_[iQuant][iNode];
	(*myHist) << time;
      
	// Iterate over all dofs
	for ( Integer iDof = 0; iDof < actDof; iDof++ ) {
	  data.Get( solTypes[iSol], histNodesPerQuant_[iQuant][iNode]-1,
		    iDof, val );

	  (*myHist) << "  " << val;
	}
	(*myHist) << std::endl;
      }
    }
  }


  void WriteResults::WriteNodeHistoryHarmonic(const NodeStoreSol<Complex>&data,
					      const Integer step,
					      const Double frequency,
					      const ComplexFormat format) {
    ENTER_FCN( "WriteResults::WriteNodeHistoryHarmonic" );

    std::ofstream * myHist;
    SolutionType actSolType;
    StdVector<SolutionType> solTypes;
    Integer iQuant, actDof;
    std::string quantity;
    Complex val;
    Double val1, val2;

    data.GetSolutionTypes(solTypes);

    // Iterate over all solutiontypes
    for (Integer iSol=0; iSol<solTypes.GetSize(); iSol++) {
   
      actDof = data.GetDof(solTypes[iSol]);
    
      // Find the related quantity
      iQuant = -1;
      for (Integer i=0; i<histQuantities_.GetSize(); i++) {

	if (histQuantities_[i] == solTypes[iSol]){
	  iQuant = i;
	  break;
	}
      }
      
      if ( iQuant != -1 ) {
    
	// Iterate over all history nodes
	for ( Integer iNode = 0; iNode < histNodesPerQuant_[iQuant].GetSize();
	      iNode++ ) {
	  myHist = historyFiles_[iQuant][iNode];
	  (*myHist) << std::scientific << frequency;
	  
	  // Iterate over all dofs
	  for ( Integer iDof = 0; iDof < actDof; iDof++ ) {
	    data.Get( solTypes[iSol], histNodesPerQuant_[iQuant][iNode]-1,
		      iDof, val );
	    
	    if (format == REAL_IMAG)
	      {
		val1 = val.real();
		val2 = val.imag();
		
		(*myHist) << "  " << val1 << "  " << val2;
	      }
	    else if (format == AMPLITUDE_PHASE) {
	      val1 = std::abs(val); 
	      if (abs(val.imag()) > 1e-16) {
		val2 = std::arg(val)*180/PI;
	      }
	      else {
		val2 = 0;
	      }
	      
	      (*myHist) << "  " << val1 << "  " << val2;
	    }
	  }
	  (*myHist) << std::endl;
	}
      }
    }

  }

} // end of namespace
