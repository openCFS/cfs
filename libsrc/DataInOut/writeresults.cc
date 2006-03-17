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
  WriteResults::WriteResults( const Char * const filename ){

    ENTER_FCN( "WriteResults::WriteResults" );

    namefile_ = filename;
    NeedHistory_ = FALSE;
    ptGrid_ = NULL;
    totalNumHistFiles_ = 0;
    
  }


  // ======================
  //   Default Destructor
  // ======================
  WriteResults::~WriteResults() {
    ENTER_FCN( "WriteResults::~WriteResults" );
  
    // delete element history files
    for (UInt i=0; i<histElemFiles_.GetSize(); i++) {
      for (UInt j=0; j<histElemFiles_[i].GetSize(); j++) {
        if (histElemFiles_[i][j]) {
          delete histElemFiles_[i][j];
        }
      }
    }

    // delete nodal history files
    for (UInt i=0; i<histNodeFiles_.GetSize(); i++) {
      for (UInt j=0; j<histNodeFiles_[i].GetSize(); j++) {
        if (histNodeFiles_[i][j]) {
          delete histNodeFiles_[i][j];
        }
      }
    }

  }


  // ==========================
  //   Initialisation Routine
  // ==========================
  void WriteResults::InitHistoryFiles() {

    ENTER_FCN( "WriteResults::InitHistoryFiles" );
    
    // Check for element history
    InitHistoryData(histElemQuantities_,
                    histElemsPerQuant_,
                    histElemFiles_,
                    ELEMENT);
    
    // Check for nodal history
    InitHistoryData(histNodeQuantities_,
                    histNodesPerQuant_,
                    histNodeFiles_,
                    NODE);
  }




  // ================================
  //  Private Initialisation Routine
  // ================================
  void WriteResults::
  InitHistoryData(StdVector<SolutionType> & quantities,
                  StdVector<StdVector<UInt> > & entitiesPerQuant,
                  StdVector<StdVector<std::ofstream*> > & histEntFiles,
                  entityType entType ) {
  
    ENTER_FCN( "WriteResults::InitHistoryData" );
 
    // *****************************************************************
    //   Determine for which output quantities there are history nodes
    //   and what their names are
    // *****************************************************************
    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<std::string> entityVec, regionVec, quantVec;
    SolutionType solType;
    RegionIdType regionId;
    
    attrVec = "", "";
    valVec = "", "";
    
    // First, get all quantities, which have to be saved

    if ( entType == ELEMENT ) {
      keyVec = "storeResults", "elemHistory", "type";
      params->GetList(keyVec, attrVec, valVec, quantVec);

      keyVec = "storeResults", "elemHistory", "saveElems";
      params->GetList(keyVec, attrVec, valVec, entityVec);

      keyVec = "storeResults", "elemHistory", "saveRegion";
      params->GetList(keyVec, attrVec, valVec, regionVec);
	}
    else {
      keyVec = "storeResults", "nodeHistory", "type";
      params->GetList(keyVec, attrVec, valVec, quantVec);

      keyVec = "storeResults", "nodeHistory", "saveNodes";
      params->GetList(keyVec, attrVec, valVec, entityVec);
      
      keyVec = "storeResults", "nodeHistory", "saveRegion";
      params->GetList(keyVec, attrVec, valVec, regionVec);
    }
    
    // Check if there are any history entites at all
    // and set NeedHistory_ flag accordingly
    if ( quantVec.GetSize() > 0 ) {
      NeedHistory_ = TRUE;
    } else {
      return;
    }

    // *****************************************************************
    //   Generate vectors holding information on the output quantities
    //   and entities for which the user wants to see the history
    //
    //   quantities .......... holds identifiers for all quantities
    //   entitiesPerQuant .... numbers of the history entities (elements, 
    //                         nodes) for each of the output quantities
    //
    // *****************************************************************

    // Make a sanity check
    if ( quantities.GetSize() != 0 || entitiesPerQuant.GetSize() != 0 ){
      Error( "Repeated call to WriteResults::InitHistoryData ?!?!",
             __FILE__, __LINE__ );
    }

    // Define some loop counters and auxilliary variables
    UInt iQuant, i, j;
    Integer quantityFound;

    // Local data structure for combining identifiers of history entites
    // with the respective entity numbers.
    // For each output quantity we store a StdVector that stores triples
    // containing
    //  a) index of identifier in entiyVec
    //  b) start index of entity numbers for identifier in entitiesPerQuant
    //  c) stop index of entity numbers for identifier in entitiesPerQuant
    StdVector< StdVector<UInt> > histEntityNumIdentCoup;

    // In the vector of quantities there may be duplicates
    // these must be eliminated!
    for ( iQuant = 0; iQuant < quantVec.GetSize(); iQuant++ ) {

      // (re-)initialise
      quantityFound = -1;

      // Determine if current output quantity was already assigned
      for ( j = 0; j < quantities.GetSize(); j++ ) {
        String2Enum( quantVec[iQuant], solType );
        if ( quantities[j] == solType ) {
          quantityFound = j;
          break;
        }
      }

      // The quantity was not yet assigned, so do this now
      // and also make room for storing the curresponding
      // number of history entities
      if ( quantityFound == -1 ) {
        String2Enum( quantVec[iQuant], solType );
        quantities.Push_back(solType);
        entitiesPerQuant.Push_back(StdVector<UInt>());
        histEntityNumIdentCoup.Push_back(StdVector<UInt>());
        quantityFound = entitiesPerQuant.GetSize() - 1;
      }

      // Whenever an output quantity appears in quantVec it is
      // related to a saveNodes/saveElements identifier. Add the 
      // corresponding nodes / elements  to the entitesPerQuant 
      //! vector and the indices of the identifier and entity numbers 
      // to histEntityNumIdentCoup
      StdVector<UInt> tempEntities;
      
      // Check, if any region name or element/node name was given
      if ( regionVec[iQuant] == ""
           && entityVec[iQuant] == "" ) {
        (*error) << "There was no region-/element-/nodeName given "
                 << "for Quantity '" << quantVec[iQuant] <<"'!\n";
        Error( __FILE__, __LINE__ );
      }

      // Get elements/regions
      if ( entType == ELEMENT ) {

        // === GET ELEMENTS ===
        StdVector<Elem*> tempElems;

        // get elements of region          
        if ( regionVec[iQuant] != "") {
          regionId = ptGrid_->RegionNameToId( regionVec[iQuant] );
		  std::cout << "Getting element for region " << regionId << std::endl;
          ptGrid_->GetElems( tempElems, regionId);
          for (UInt iElem = 0; iElem < tempElems.GetSize(); iElem++) {
            tempEntities.Push_back(tempElems[iElem]->elemNum);
          }
		  std::cout << "-> found " << tempEntities.GetSize() << std::endl;
        }
        // get elements by name        
        if ( entityVec[iQuant] != "" ) {
          ptGrid_->GetElemsByName( tempElems, entityVec[iQuant] );
          for (UInt iElem = 0; iElem < tempElems.GetSize(); iElem++) {
            tempEntities.Push_back(tempElems[iElem]->elemNum);
          }
        }
      } else {
        // == GET NODES ==
        StdVector<UInt> tempNodes;
        
        // get nodes of region
        if ( regionVec[iQuant] != "") {
          regionId = ptGrid_->RegionNameToId( regionVec[iQuant] );
          ptGrid_->GetNodesByRegion( tempNodes, regionId);
          for (UInt iNode = 0; iNode < tempNodes.GetSize(); iNode++) {
            tempEntities.Push_back( tempNodes[iNode] );
          }
        }

        // get nodes by name
        if ( entityVec[iQuant] != "" ) {
          ptGrid_->GetNodesByName( tempEntities, entityVec[iQuant] );
          for (UInt iNode = 0; iNode < tempNodes.GetSize(); iNode++) {
            tempEntities.Push_back( tempNodes[iNode] );
          }
        }
      }
      
      // sum up number of eneities ( which corresponds to the number
      // of file streams to open)
      totalNumHistFiles_ += tempEntities.GetSize();

      histEntityNumIdentCoup[quantityFound].Push_back(iQuant);
      histEntityNumIdentCoup[quantityFound].
        Push_back(entitiesPerQuant[quantityFound].GetSize());
      for ( i = 0; i < tempEntities.GetSize(); i++ ) {
        entitiesPerQuant[quantityFound].Push_back(tempEntities[i]);
      }
      histEntityNumIdentCoup[quantityFound].
        Push_back(entitiesPerQuant[quantityFound].GetSize()-1);
    }

    // Check if the maximum number of open file streams would be 
    // violated
    if ( totalNumHistFiles_ > MAX_NUM_HIST_FILES ) {
      (*error) << "Due to file system limitations only "
               << MAX_NUM_HIST_FILES << " history nodes/elements "
               << "can be written.\nHowever, you want to write "
               << totalNumHistFiles_ << " results.";
      Error( __FILE__, __LINE__ );
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
    histEntFiles.Resize(quantities.GetSize());

    // Iterate over all output quantities
    std::string quantityName;
    for ( iQuant = 0; iQuant < quantities.GetSize(); iQuant++ ) {

      // Log information on this output quantity to info file.
      // This will also tell the user the relation between the
      // entity numbers and the identifiers
      // std::string comment = "Identifiers for history nodes of quantity '";
      Enum2String( quantities[iQuant], quantityName );

      if ( entType == ELEMENT ) {
        Info->PrintF( "", "\nHistory elements for quantity '%s':\n",
                      quantityName.c_str() );
      } else {
        Info->PrintF( "", "\nHistory nodes for quantity '%s':\n",
                      quantityName.c_str() );
      }

      for ( i = 0; i < histEntityNumIdentCoup[iQuant].GetSize() / 3; i++ ) {

        // *** ELEMENT HISTORY ***
        if ( entType == ELEMENT ) {

          if (entityVec[histEntityNumIdentCoup[iQuant][i*3]] != "" ) {
            Info->PrintF( "", " identifier = %s, element(s) =", 
                          entityVec[histEntityNumIdentCoup[iQuant][i*3]].c_str() );
          } else {
            Info->PrintF( "", " identifier = %s, element(s) =", 
                          regionVec[histEntityNumIdentCoup[iQuant][i*3]].c_str() );
          }

        } else {
          // *** NODE HISTORY ***
          if (entityVec[histEntityNumIdentCoup[iQuant][i*3]] != "" ) {
            Info->PrintF( "", " identifier = %s, node(s) =", 
                          entityVec[histEntityNumIdentCoup[iQuant][i*3]].c_str() );
          } else {
            Info->PrintF( "", " identifier = %s, node(s) =", 
                          regionVec[histEntityNumIdentCoup[iQuant][i*3]].c_str() );
          }
        }

        for ( j  = histEntityNumIdentCoup[iQuant][i*3+1];
              j <= histEntityNumIdentCoup[iQuant][i*3+2]; j++ ) {
          Info->PrintF( "", " %d", entitiesPerQuant[iQuant][j] );
        }
        Info->PrintF( "", "\n" );
      }

      // For each node there will be one history file
      histEntFiles[iQuant].Resize(entitiesPerQuant[iQuant].GetSize());

      for ( UInt iEntity = 0; iEntity <entitiesPerQuant[iQuant].GetSize();
            iEntity++ ) {

        // Create complete filename
        std::string quantString;
        Enum2String(quantities[iQuant], quantString);
        totalName = namePrefix + quantString;
        if ( entType == ELEMENT ) {
          totalName += "-elem-";
        } else {
          totalName += "-node-";
        }
        totalName += GenStr(entitiesPerQuant[iQuant][iEntity]);
        totalName += namePostfix;

        // Open output file
        histEntFiles[iQuant][iEntity] = NULL;
        histEntFiles[iQuant][iEntity] = new std::ofstream(totalName.c_str());
        if ( !histEntFiles[iQuant][iEntity] ) {
          std::string errMsg = "InitHistoryData: Cannot open history file '";
          errMsg += totalName + "'";
          Error (errMsg.c_str(), __FILE__, __LINE__);
        }
      }
    }
  }



  void WriteResults::WriteSolMatrix(Grid * ptgrid, 
                                    const Vector<Double> sol, 
                                    const std::string matFileName,
                                    const UInt nrDofs) {

    ENTER_FCN( "WriteResults::WriteSolMatrix" );

    //get and write number of nodes on the level
    UInt numnodes=ptgrid->GetNumNodes();
    UInt dim=ptgrid->GetDim();
    UInt i;

    std::ofstream * matrixOut = new std::ofstream(matFileName.c_str());
    // use scientific output, formatting is much better
    matrixOut->setf(std::ios::scientific, std::ios::floatfield);
  
    if (dim==2) {
      Point<2> point;
      
      // write x,y,z-coordinate
      for (i=0; i<numnodes; i++) {
        ptgrid->GetNodeCoordinate(point,i);
        (*matrixOut) << " \t" << point[0] << " \t" << point[1] << " \t"
                     << 0 << " \t";

        for (UInt actDof =0; actDof < nrDofs; actDof++) {
          (*matrixOut) << sol[i*nrDofs + actDof] << " \t";
        }
        (*matrixOut) << std::endl;
      }
    }
    else {
      Point<3> point;
      
      // write x,y,z-coordinate
      for (i=0; i<numnodes; i++) {
        ptgrid->GetNodeCoordinate(point,i);
        (*matrixOut) << " \t" << point[0] << " \t" << point[1] << " \t"
                     << point[2] << " \t";

        for (UInt actDof =0; actDof < nrDofs; actDof++) {
          (*matrixOut) << sol[i*nrDofs + actDof] << " \t";
        }
        (*matrixOut) << std::endl;
      } 
    }
  }


  void WriteResults::WriteNodeHistoryTransient(const NodeStoreSol<Double>&data,
                                               const UInt step,
                                               const Double time) {

    ENTER_FCN( "WriteResults::WriteNodeHistoryTransient" );

    std::ofstream * myHist;
    StdVector<SolutionType> solTypes;
    Integer iQuant; 
    UInt actDof;
    std::string quantity;
    Double val;

    data.GetSolutionTypes(solTypes);
  
    // Iterate over all solutiontypes
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
   
      actDof = data.GetDof(solTypes[iSol]);
    
      // Find the related quantity
      iQuant = histNodeQuantities_.Find(solTypes[iSol]);
      if ( iQuant != -1 ) {
        
        // Iterate over all history nodes
        for ( UInt iNode = 0; iNode < histNodesPerQuant_[iQuant].GetSize();
              iNode++ ) {
          myHist = histNodeFiles_[iQuant][iNode];
          (*myHist) << time;
          
          // Iterate over all dofs
          for ( UInt iDof = 0; iDof < actDof; iDof++ ) {
            data.Get( solTypes[iSol], histNodesPerQuant_[iQuant][iNode]-1,
                      iDof, val );
            
            (*myHist) << "  " << val;
          }
          (*myHist) << std::endl;
        }
      }
    }
  }

  void WriteResults::WriteElemHistoryTransient(const ElemStoreSol<Double>&data,
                                               const UInt step,
                                               const Double time) {
    
    ENTER_FCN( "WriteResults::WriteElemHistoryTransient" );

    std::ofstream * myHist;
    StdVector<SolutionType> solTypes;
    Integer iQuant; 
    UInt actDof;
    std::string quantity;
    Double val;
    
    data.GetSolutionTypes(solTypes);

    // Iterate over all solutiontypes
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
   
      actDof = data.GetDof(solTypes[iSol]);
      
      // Find the related quantity
      iQuant = histElemQuantities_.Find(solTypes[iSol]);
      if ( iQuant != -1 ) {
        
        // Iterate over all history elements
        for ( UInt iElem = 0; iElem < histElemsPerQuant_[iQuant].GetSize();
              iElem++ ) {
          myHist = histElemFiles_[iQuant][iElem];
          (*myHist) << time;
          
          // Iterate over all dofs
          for ( UInt iDof = 0; iDof < actDof; iDof++ ) {
            data.Get( solTypes[iSol], histElemsPerQuant_[iQuant][iElem]-1,
                      iDof, val );
            
            (*myHist) << "  " << val;
          }
          (*myHist) << std::endl;
        }
      }
    }
  }



  void WriteResults::WriteNodeHistoryHarmonic(const NodeStoreSol<Complex>&data,
                                              const UInt step,
                                              const Double frequency,
                                              const ComplexFormat format) {

    ENTER_FCN( "WriteResults::WriteNodeHistoryHarmonic" );

    std::ofstream * myHist;
    StdVector<SolutionType> solTypes;
    Integer iQuant;
    UInt actDof;
    std::string quantity;
    Complex val;
    Double val1, val2;

    data.GetSolutionTypes(solTypes);

    // Iterate over all solutiontypes
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
   
      actDof = data.GetDof(solTypes[iSol]);
    
      // Find the related quantity
      iQuant = histNodeQuantities_.Find(solTypes[iSol]);
      if ( iQuant != -1 ) {    
        // Iterate over all history nodes
        for ( UInt iNode = 0; iNode < histNodesPerQuant_[iQuant].GetSize();
              iNode++ ) {
          myHist = histNodeFiles_[iQuant][iNode];
          (*myHist) << std::scientific << frequency;
          
          // Iterate over all dofs
          for ( UInt iDof = 0; iDof < actDof; iDof++ ) {
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
  

  void WriteResults::WriteElemHistoryHarmonic(const ElemStoreSol<Complex>&data,
                                              const UInt step,
                                              const Double frequency,
                                              const ComplexFormat format) {
    ENTER_FCN( "WriteResults::WriteElemHistoryHarmonic" );
    std::ofstream * myHist;
    StdVector<SolutionType> solTypes;
    Integer iQuant;
    UInt actDof;
    std::string quantity;
    Complex val;
    Double val1, val2;

    data.GetSolutionTypes(solTypes);

    // Iterate over all solutiontypes
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
   
      actDof = data.GetDof(solTypes[iSol]);
    
      // Find the related quantity
      iQuant = histElemQuantities_.Find(solTypes[iSol]);
      if ( iQuant != -1 ) {    
        // Iterate over all history elements
        for ( UInt iElem = 0; iElem < histElemsPerQuant_[iQuant].GetSize();
              iElem++ ) {
          myHist = histElemFiles_[iQuant][iElem];
          (*myHist) << std::scientific << frequency;
          
          // Iterate over all dofs
          for ( UInt iDof = 0; iDof < actDof; iDof++ ) {
            data.Get( solTypes[iSol], histElemsPerQuant_[iQuant][iElem]-1,
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
