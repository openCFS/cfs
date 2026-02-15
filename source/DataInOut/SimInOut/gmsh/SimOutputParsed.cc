// =====================================================================================
// 
//       Filename:  SimOutputParsed.cc
// 
//    Description:  Implementation file for the POS writer
//                  We proceed like the following:
//
//                  The first approach to a transient simulation will be,
//                  that each timestep will get its own file storing also the mesh.
//
//                  for each result in each timestep we create a new file naming
//                  <SimName>_<ResultName>.pos
//
//                  Each region is represented with a view. The used may specify
//                  if the interpolation matrix is to be stored for each element which
//                  which will be very complex or constant within each region (most cases)
//
//                  To create a view we do the following:
//                   1. if we store region constant storage
//                      a. loop over each element
//                         get its type and its interpolation matrices
//                         if we do not have this elementtype we store it in a map
//                         and save its elemIdx
//                   2. if we have storage for each element
//                      proceed like above but associate for each element this datastructure
//
//                  When we write out the results, we create the elements and fill the timestep
//                  values with undefinied variables of the form:
//                  x_<eqnNr>_<stepnum>
//                  in each timestep we fill the variables to the top of the file
//                  This means, that gmsh will open the file, displays some errors about undefinied
//                  variables and assumes 0 as a value.
//                  The lower part of the file will not be altered lateron.
//
// 
//        Version:  1.0
//        Created:  11/19/2011 04:56:24 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#include "SimOutputParsed.hh"

#include "Domain/Mesh/Grid.hh"

#include <sstream>
#include <cmath>

namespace CoupledField{

  SimOutputParsed::SimOutputParsed(std::string fileName,
                                 PtrParamNode outputNode,
                                 PtrParamNode infoNode, bool isRestart ) :
                SimOutput( fileName, outputNode, infoNode, isRestart )
  {

    // The restart case is currently not implemented, i.e. resuls from a 
    // partial simulation get overwritten.
    if( isRestart_ ) {
      WARN( "The GMSH-Parsed-Writer is currently not adapted to write restarted "
            "results correctly, thus the results of the previous run get"
            " overwritten." );
    }
    
    formatName_ = "gmshParsed";

    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    dirName_ = "results_" + formatName_;

    try{
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }
    firstStep_ = true;

  }

  SimOutputParsed::~SimOutputParsed(){
    //delete input file streams
    std::map<SolutionType, std::fstream*>::iterator inIt = infiles_.begin();
    while(inIt != infiles_.end()){
      inIt->second->close();
      delete inIt->second;
      inIt->second = NULL;
      inIt++;
    }
    infiles_.clear();

    //delete output file streams
    std::map<SolutionType, std::fstream*>::iterator outIt = outfiles_.begin();
    while(outIt != outfiles_.end()){
      outIt->second->close();
      delete outIt->second;
      outIt->second = NULL;
      outIt++;
    }
    outfiles_.clear();
  }

  void SimOutputParsed::Init(Grid* ptGrid, bool printGridOnly ){
    myGrid_ = ptGrid;
    FillElemMapping();
  }

  //! Begin multisequence step
  void SimOutputParsed::BeginMultiSequenceStep( UInt step,
                                       BasePDE::AnalysisType type,
                                       UInt numSteps){
    numSteps_ = numSteps;
    SolutionRegionMap::iterator solIt = mySolutions_.begin();
    while(solIt!=mySolutions_.end()){
      RegionInterpMap::iterator regIter = solIt->second.begin();
      while(regIter!=solIt->second.end()){
        std::list<ElemInterpolation>::iterator interIter = regIter->second.begin();
        while(interIter!= regIter->second.end()){
          interIter->saveEnd_ = (numSteps_<interIter->saveEnd_)?  numSteps_ : interIter->saveEnd_;
          interIter++;
        }
        regIter++;
      }
      solIt++;
    }
  }

  //! Register result (within one multisequence step)
  void SimOutputParsed::RegisterResult( shared_ptr<BaseResult> sol,
                               UInt saveBegin, UInt saveInc,
                               UInt saveEnd,
                               bool isHistory ){
    //obtain the entity list of the result
    shared_ptr<EntityList> curList = sol->GetEntityList();
    std::map<Elem::ShapeType,std::list<UInt> > sortedElems;

    if(curList->GetDefineType() == EntityList::REGION){

      shared_ptr<BaseFeFunction> curFnc = sol->GetResultInfo()->GetFeFunction().lock();

      if(!curFnc){
        WARN("Going to ignore result as feFunction is not set");
        return;
      }
      shared_ptr<EntityList> actList = myGrid_->GetEntityList( EntityList::ELEM_LIST,
                                                             curList->GetName() );
      EntityIterator entIt = actList->GetIterator();
      for(;!entIt.IsEnd();entIt++){
        Elem::ShapeType curType = entIt.GetElem()->GetShapeType(entIt.GetElem()->type);
        sortedElems[curType].push_back(entIt.GetElem()->elemNum);
      }
      //That was easy, now we add our lists
      std::map<Elem::ShapeType,std::list<UInt> >::iterator eIt = sortedElems.begin();
      while(eIt != sortedElems.end()){
        ElemInterpolation curInterp;
        curInterp.eIndices_ = eIt->second;
        curInterp.curType = eTypeMap_[eIt->first];
        //rude but works
        if(sol->GetResultInfo()->entryType == ResultInfo::VECTOR){
          UInt enumNr = (UInt)eTypeMap_[eIt->first];
          enumNr++;
          curInterp.curType = (GmeshParsedElemTypes)enumNr;
        }else if(sol->GetResultInfo()->entryType == ResultInfo::TENSOR){
          UInt enumNr = (UInt)eTypeMap_[eIt->first];
          enumNr +=2;
          curInterp.curType = (GmeshParsedElemTypes)enumNr;
        }
        curInterp.saveBegin_ = saveBegin;
        curInterp.saveInc_ = saveInc;
        curInterp.saveEnd_ = saveEnd;
        mySolutions_[sol->GetResultInfo()->resultType][actList->GetRegion()].push_back(curInterp);
        eIt++;
      }
    }else{
      Exception("Define Type of result is not supported by Writer");
    }

  }

  //! Begin single analysis step
  void SimOutputParsed::BeginStep( UInt stepNum, Double stepVal ){
    curStep_ = stepNum;
    curTime_ = stepVal;
  }

  //! Add result to current step
  void SimOutputParsed::AddResult( shared_ptr<BaseResult> sol ){
    shared_ptr<BaseFeFunction> myFunct = sol->GetResultInfo()->GetFeFunction().lock();


    if(!myFunct){
      return;
    }
    //cast down to double FeFunction
    FeFunction<Double>* feFnc = dynamic_cast<FeFunction<Double>* >(myFunct.get());
    if(firstStep_){
      //now we already prepare everything
      PrepareResultFile(sol);
      firstStep_ = false;
    }
    std::fstream* outFile =  outfiles_[sol->GetResultInfo()->resultType];
    std::fstream* inFile =  infiles_[sol->GetResultInfo()->resultType];
    if(!outFile){
      WARN("cannot obtain file pointer");
      return;
    }
    //now loop over every element in Region Interp map
    std::stringstream solVars;

    //write down current time variable
    solVars << "t_" << curStep_ << "=" << curTime_ << ";" << std::endl;

    RegionInterpMap::iterator regInterIter = mySolutions_[sol->GetResultInfo()->resultType].begin();
    while(regInterIter!= mySolutions_[sol->GetResultInfo()->resultType].end()){
      std::list<ElemInterpolation>::iterator eInter = regInterIter->second.begin();
      while(eInter!=regInterIter->second.end()){
        std::list<UInt>::iterator elemIter = (*eInter).eIndices_.begin();
        while(elemIter!=(*eInter).eIndices_.end()){
          const Elem* curE = myGrid_->GetElem(*elemIter);
          Vector<Double> eSol;
          StdVector<Integer> eqns;
          feFnc->GetFeSpace()->GetElemEqns(eqns,curE);
          feFnc->GetElemSolution(eSol,curE);
          for(UInt i = 0;i<eSol.GetSize();i++){
            solVars << "x_" << eqns[i] << "_" << curStep_ << "=" << eSol[i] << ";" << std::endl;
          }
          elemIter++;
        }
        eInter++;
      }
      regInterIter++;
    }
    //outFile->seekp(10,std::ios::end);
    PutVarsToResultFile(outFile,inFile, solVars.str(),resFilePosition_[sol->GetResultInfo()->resultType]);
  }

  //! End single analysis step
  void SimOutputParsed::FinishStep( ){

  }

  //! End multisequence step
  void SimOutputParsed::FinishMultiSequenceStep( ){

  }

  //! Finalize the output
  void SimOutputParsed::Finalize(){

  }

  void SimOutputParsed::FillElemMapping(){
    GmeshParsedElemTypesEnum.Add(SP,"SP");
    GmeshParsedElemTypesEnum.Add(VP,"VP");
    GmeshParsedElemTypesEnum.Add(TP,"TP");
    GmeshParsedElemTypesEnum.Add(SL,"SL");
    GmeshParsedElemTypesEnum.Add(VL,"VL");
    GmeshParsedElemTypesEnum.Add(TL,"TL");
    GmeshParsedElemTypesEnum.Add(ST,"ST");
    GmeshParsedElemTypesEnum.Add(VT,"VT");
    GmeshParsedElemTypesEnum.Add(TT,"TT");
    GmeshParsedElemTypesEnum.Add(SQ,"SQ");
    GmeshParsedElemTypesEnum.Add(VQ,"VQ");
    GmeshParsedElemTypesEnum.Add(TQ,"TQ");
    GmeshParsedElemTypesEnum.Add(SS,"SS");
    GmeshParsedElemTypesEnum.Add(VS,"VS");
    GmeshParsedElemTypesEnum.Add(TS,"TS");
    GmeshParsedElemTypesEnum.Add(SH,"SH");
    GmeshParsedElemTypesEnum.Add(VH,"VH");
    GmeshParsedElemTypesEnum.Add(TH,"TH");
    GmeshParsedElemTypesEnum.Add(SI,"SI");
    GmeshParsedElemTypesEnum.Add(VI,"VI");
    GmeshParsedElemTypesEnum.Add(TI,"TI");
    GmeshParsedElemTypesEnum.Add(SY,"SY");
    GmeshParsedElemTypesEnum.Add(VY,"VY");
    GmeshParsedElemTypesEnum.Add(TY,"TY");
    GmeshParsedElemTypesEnum.Add(T2,"T2");
    GmeshParsedElemTypesEnum.Add(T3,"T3");


    //we store only the enum of the scalar results
    //to the a valid enum index later when we notice
    //a vector or tensor result we add 1 or 2 respecively
    eTypeMap_[Elem::ST_UNDEF] = UNDEF;
    eTypeMap_[Elem::ST_POINT] = SP;
    eTypeMap_[Elem::ST_LINE] = SL;
    eTypeMap_[Elem::ST_TRIA] = ST;
    eTypeMap_[Elem::ST_QUAD] = SQ;
    eTypeMap_[Elem::ST_TET] = SS;
    eTypeMap_[Elem::ST_HEXA] = SH;
    eTypeMap_[Elem::ST_PYRA] = SY;
    eTypeMap_[Elem::ST_WEDGE] = SI;
  }

  void SimOutputParsed::PrepareResultFile(shared_ptr<BaseResult> sol){
    //TODO: CLEAN UP THIS MESS
    //  ok this function is very big which does not has to be the case.
    //  So along with the time marching schemes we will refactor this and perhaps think about
    //  new data-structures

    //now lets go over each solution and create a file
    SolutionRegionMap::iterator solIt = mySolutions_.begin();
    StdVector<std::string> regionNames;
    myGrid_->GetRegionNames(regionNames);

    while(solIt != mySolutions_.end()){
      std::string solName = SolutionTypeEnum.ToString(solIt->first);

      std::string name = fileName_ + "_" + solName + ".pos";

      // Generate basename for output file
      fs::path filePath = dirName_ / name;
      
      outfiles_[solIt->first]  = new std::fstream();
      infiles_[solIt->first]  = new std::fstream();

      std::fstream* outFile =  outfiles_[solIt->first];
      std::fstream* inFile =  infiles_[solIt->first];

      outFile->open(filePath.c_str(),std::ios::trunc | std::ios::out);
      inFile->open(filePath.c_str(),std::ios::in);

      if ( !outFile && outFile->is_open() ) {
        EXCEPTION("Could not open file " << filePath
                  << " for writing parsed Gmsh output");
      }

      if ( !inFile && inFile->is_open() ) {
        EXCEPTION("Could not open file " << filePath
                  << " for reading parsed Gmsh output");
      }
      //write a newline at the beginning of the file, just in case...
      (*outFile) << std::endl;

      //file is opened now loop over the regions and create a view for each region
      RegionInterpMap curregInter = solIt->second;
      RegionInterpMap::iterator curIt = curregInter.begin();
      shared_ptr<BaseFeFunction> myFunct = sol->GetResultInfo()->GetFeFunction().lock();

      while(curIt != curregInter.end()){
        (*outFile) << "View \"" << regionNames[curIt->first] << "\" {" << std::endl;
        //now we loop over the InterpolationList
        std::list<ElemInterpolation>::iterator elemIt = curIt->second.begin();
        while(elemIt != curIt->second.end()){
          ElemInterpolation eInterpol = (*elemIt);
          std::string eType = GmeshParsedElemTypesEnum.ToString(eInterpol.curType);
          std::list<UInt>::iterator eIdx = eInterpol.eIndices_.begin();
          while(eIdx != eInterpol.eIndices_.end()){
            //get element corner coords
            const Elem* curE = myGrid_->GetElem(*eIdx);
            if(!eInterpol.matrixComputed_){
              myFunct->GetFeSpace()->GetFe(curE->elemNum)->ComputeMonomialCoefficients(eInterpol.exponents_,eInterpol.coefs_);
              eInterpol.matrixComputed_ = true;
            }
            UInt numVert = Elem::shapes[curE->type].numVertices;
            StdVector<UInt> elemNodes = curE->connect;
            Vector<Double> curCoord;
            std::stringstream coordstream;

            for(UInt i = 0;i<numVert;i++){
              myGrid_->GetNodeCoordinate3D(curCoord,elemNodes[i]);
              for(UInt d = 0;d<curCoord.GetSize(); d++){
                coordstream << curCoord[d] << ",";
              }
              coordstream << std::endl;
            }
            std::string cSTR = coordstream.str();
            //erase last two characters
            cSTR.erase(--cSTR.end());
            cSTR.erase(--cSTR.end());

            (*outFile) << eType << "(" << std::endl;
            (*outFile) << cSTR << std::endl;
            (*outFile) << "){" << std::endl;
            //remember the current position
            WriteDummyResults(curE,eInterpol,myFunct,outFile);
            (*outFile) << "};" << std::endl;
            eIdx++;
          }
          (*outFile) << "TIME{" ;
          std::stringstream timeStr;
          for(UInt step=0;step<eInterpol.saveEnd_;step+=eInterpol.saveInc_){
            timeStr << "t_" << step+1 << ",";
          }
          std::string tStr  = timeStr.str();
          tStr.erase(--tStr.end());
          (*outFile) << tStr << "};" << std::endl;
          std::string interpol;
          GetInterpolationString(interpol,eInterpol);
          (*outFile) << interpol << std::endl;
          elemIt++;
        }
        (*outFile) << "};" << std::endl;
        curIt++;
      }
      resFilePosition_[solIt->first] = 0;
      solIt++;
    }
  }
  void SimOutputParsed::WriteDummyResults(const Elem* elem,
                                               ElemInterpolation& eInterpol,
                                               shared_ptr<BaseFeFunction> feFnc,
                                               std::fstream* out){
     StdVector<Integer> eqns;
     feFnc->GetFeSpace()->GetElemEqns(eqns,elem);
     std::stringstream oStream;

     if(feFnc->GetResultInfo()->entryType == ResultInfo::SCALAR || feFnc->GetResultInfo()->dofNames.GetSize() == 3){
       for(UInt step=eInterpol.saveBegin_;step<eInterpol.saveEnd_;step+=eInterpol.saveInc_){
         for(UInt i =0;i<eqns.GetSize();i++){
           oStream << "x_" << eqns[i] << "_" << step+1 << ",";
         }
       }
     }else if (feFnc->GetResultInfo()->entryType == ResultInfo::VECTOR && feFnc->GetResultInfo()->dofNames.GetSize() == 2){
       UInt d = 0;
       for(UInt step=eInterpol.saveBegin_;step<eInterpol.saveEnd_;step+=eInterpol.saveInc_){
         for(UInt i =0;i<eqns.GetSize();i++){
           if(d == 2){
             oStream << "0" << ",";
             i--;
             d=0;
           }else{
             oStream << "x_" << eqns[i] << "_" << step+1 << ",";
             d++;
           }
         }
         //put final zero
         oStream << "0" << ",";
         d=0;
       }
     }else{
       EXCEPTION("TENSOR RESULTS NOT CLEANLY SUPPORTED in 2D RIGHT NOW");
     }

     std::string oSTR = oStream.str();
     oSTR.erase(--oSTR.end());
     (*out) << oSTR << std::endl;
  }

  void SimOutputParsed::GetInterpolationString(std::string & interp,
                                                     ElemInterpolation eInterpol){
    std::stringstream cStream;
    std::stringstream expStream;

    cStream << std::endl;
    cStream.precision(10);
    for(UInt i = 0; i<eInterpol.coefs_.GetNumRows();i++ ){
      cStream << "{";
      for(UInt j = 0; j<eInterpol.coefs_.GetNumCols();j++ ){
        cStream << std::scientific << eInterpol.coefs_[i][j] << ",";
      }
      cStream.seekp(-1,std::ios::end);
      cStream << "}," << std::endl;
    }
    cStream.seekp(-2,std::ios::end);
    cStream << " " << std::endl;

    for(UInt i = 0; i<eInterpol.exponents_.GetNumRows();i++ ){
      expStream << "{";
      for(UInt j = 0; j<eInterpol.exponents_.GetNumCols();j++ ){
        expStream << eInterpol.exponents_[i][j] << ",";
      }
      expStream.seekp(-1,std::ios::end);
      expStream << "}," << std::endl;
    }
    expStream.seekp(-2,std::ios::end);
    expStream << " " << std::endl;

    interp = "INTERPOLATION_SCHEME{" + cStream.str() + "}{\n" + expStream.str() + "};";
  }

  void SimOutputParsed::PutVarsToResultFile(std::fstream* outfile, std::fstream* infile,
                                            std::string vars,long& destination){
    //this is based on the code found at
    //http://www.codeproject.com/KB/cs/InsertTextInCSharp.aspx

    long DEFAULT_BUFFER = 32 * 1024;

     //determine size of string
    long varsSize = vars.size();

    //determine filelength before resize
    outfile->seekp(0,std::ios::end);
    long origFileSize = outfile->tellp();
    long lengthMove = origFileSize - destination;

    //determine the end of the range to be inserted. i.e. new destination for the next time
    long newDest = destination+varsSize;

    //extend output file
    //outfile->seekp(varsSize,std::ios::end);
    //(*outfile) << std::endl;
    outfile->flush();
    //start copying file contents
    //define current read write positions
    long readPos = destination;// +Length;
    long writePos = newDest;// +Length;

    //define the buffersize, for first test we hardcode to 32KB
    int bfrS = (lengthMove < DEFAULT_BUFFER)? lengthMove : DEFAULT_BUFFER;
    //this second test is to determine if the number of bytes inserted
    //exceeds 2times the buffer length. if so, we set the bufS to the next available
    //power of two
    if(bfrS < (varsSize)){
      //determine next power of two available
      //cut float values
      long exp = long(std::ceil(log2(varsSize)));
      bfrS = int(std::pow(2.0,static_cast<Double>(exp)));
    }

    //number of read operations
    long numRead = lengthMove / bfrS;

    //create buffers
    char* buff = new char[bfrS];
    char* buff2 = new char[bfrS];

    int remainingBytes = (int)lengthMove % bfrS;

    infile->seekg(readPos,std::ios::beg);
    outfile->seekp(writePos,std::ios::beg);
    infile->read(buff2,bfrS);
    //long in=0,out=0;
    //in = infile->tellg();
    //out = outfile->tellp();


    //now loop over and move the contents
    for(long i = 1; i < numRead; i++){
      std::memcpy(buff,buff2,bfrS);
      infile->read(buff2,bfrS);
      outfile->write(buff,bfrS);
    }
    std::memcpy(buff,buff2,bfrS);
    if(remainingBytes>0){
      delete[] buff2;
      buff2 = NULL;
      buff2 = new char[remainingBytes];
      infile->read(buff2,remainingBytes);
      outfile->write(buff,bfrS);
      bfrS = remainingBytes;
      delete[] buff;
      buff = NULL;
      buff = new char[bfrS];
      std::memcpy(buff,buff2,bfrS);
    }
    outfile->write(buff,bfrS);

    // ok everything is shifted we write our string
    outfile->seekg(destination,std::ios::beg);
    (*outfile) << vars;
    //update our destination
    destination = newDest;
    //clear pointers
    delete[] buff;
    delete[] buff2;
    buff = NULL;
    buff2 = NULL;
  }
}
