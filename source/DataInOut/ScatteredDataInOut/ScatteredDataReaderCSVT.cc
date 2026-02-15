#include <fstream>
#include <filesystem>
#include <boost/tokenizer.hpp>

#include "General/Exception.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "PDE/BasePDE.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "ScatteredDataReaderCSVT.hh"
#include "Utils/mathParser/mathParser.hh"


namespace CoupledField 
{
  ScatteredDataReaderCSVT::ScatteredDataReaderCSVT(PtrParamNode& scatteredDataNode,
                                                 bool verbose) :
    ScatteredDataReader(scatteredDataNode, verbose)
  {
	  mp_ = domain->GetMathParser();
	  mHandleStep_ = mp_->GetNewHandle(true);
	  readerMode_ = TF;
	  std::string var;
	  if(domain->GetDriver()->GetAnalysisType() == BasePDE::HARMONIC)
	    var = "f";
	  else
	    var = "t";
	  mp_->SetExpr(mHandleStep_, var);

	  step_  = -1;
  }
  
  ScatteredDataReaderCSVT::~ScatteredDataReaderCSVT()
  {
  }

  void ScatteredDataReaderCSVT::ParseParamNode() // hier wird das XML geparst
  {
    fileName_ = myParamNode_->Get("fileName")->As<std::string>();
    std::string id = myParamNode_->Get("id")->As<std::string>();

    ParamNodeList coordList;
    std::set<std::string> dofs;
    coordList = myParamNode_->Get("coordinates")->GetList("comp");
    for(UInt i=0, n=coordList.GetSize(); i<n; i++) {
      std::string dof = coordList[i]->Get("dof")->As<std::string>();

      if(dofs.find(dof) != dofs.end()) 
      {
        EXCEPTION("Dof '" << dof << "' specified multiple times for "
                  "coordinates of scattered data reader id '" << id << "'." );
      }
      else
      {
        dofs.insert(dof);
      }
        
      if( dof == "x" ) {
        dof2CoordColumn_[0] = coordList[i]->Get("col")->As<UInt>();
      }
      if( dof == "y" ) {
        dof2CoordColumn_[1] = coordList[i]->Get("col")->As<UInt>();
      }
      if( dof == "z" ) {
        dof2CoordColumn_[2] = coordList[i]->Get("col")->As<UInt>();
      }
    }
      
    ParamNodeList quantityList;
    quantityList = myParamNode_->GetList("quantity");
    for(UInt j=0, m=quantityList.GetSize(); j<m; j++) {
      std::string qid = quantityList[j]->Get("id")->As<std::string>();

      dofs.clear();
      ParamNodeList valueList;
      valueList = quantityList[j]->GetList("comp");

      if(valueList.GetSize() == 1){
        qidDof2Column_[qid][0] = valueList[0]->Get("col")->As<UInt>();
      }else{

        for(UInt i=0, n=valueList.GetSize(); i<n; i++) {


          std::string dof = valueList[i]->Get("dof")->As<std::string>();


          if(dofs.find(dof) != dofs.end())
          {
            EXCEPTION("Dof '" << dof << "' specified multiple times for scattered "
                      << "data quantity id '" << qid << "'." );
          }
          else
          {
            dofs.insert(dof);
          }

          if( dof == "x" ) {
            qidDof2Column_[qid][0] = valueList[i]->Get("col")->As<UInt>();
          }
          if( dof == "y" ) {
            qidDof2Column_[qid][1] = valueList[i]->Get("col")->As<UInt>();
          }
          if( dof == "z" ) {
            qidDof2Column_[qid][2] = valueList[i]->Get("col")->As<UInt>();
          }
        }
      }
    }

    time2Col_ = myParamNode_->Get("stepValues")->Get("col")->As<UInt>();
    file2Col_ = myParamNode_->Get("stepFiles")->Get("col")->As<UInt>();

  }

  void ScatteredDataReaderCSVT::ReadData()
  {
    if(registeredQuantities_.empty())
    {
      return;
    }

    UInt curStep = domain->GetBasePDE()->GetSolveStep()->GetActStep();
    UInt totLines = 0;
    //return if reader is up to date
    if(step_==(Integer)curStep){
      return;
    }

    if(step_==-1){
      ParseParamNode();
      
      // Organiztion file
      typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
      std::string row;
      UInt line = 0;

      // Open CSV file
      std::fstream orgfile(fileName_.c_str(), std::ios::in);

      if(!orgfile)
      {
        EXCEPTION("Scattered data organization file '" << fileName_ << "' could not be opened!")
      }

      totLines = 0;
      while(std::getline(orgfile, row))
      {
        ++totLines;
      }

      valuetime_.Resize(totLines-skipLines_-1);
      valueFileName_.Resize(totLines-skipLines_-1);

      orgfile.clear();
      orgfile.seekg(0,std::ios::beg);

      // Iterate over lines in CSV file
      while(std::getline(orgfile, row))
      {
        line++;

        // Skip header lines.
        if(line <= skipLines_) continue;
        if(line <= skipLines_+1){
          fileNameGeom_ = row;
          continue;
        }

        // Tokenize row into string tokens.
        Tokenizer tokens(row, boost::escaped_list_separator<char>('\\', ',', '\"'));

        if(!(std::distance(tokens.begin(), tokens.end())== std::max(time2Col_+1,file2Col_+1)))
        {
          EXCEPTION("Columns in file '" << fileName_ << "' are defect, or are not associated correctly in the XML!")
        }

        Tokenizer::iterator tkIt(tokens.begin());
        std::advance(tkIt,time2Col_);
        valuetime_[(line - skipLines_-2)] = boost::lexical_cast<Double,std::string>(*tkIt);
        std::advance(tkIt,file2Col_);
        valueFileName_[line - skipLines_-2] = (*tkIt);
      }

      orgfile.close();

      // Geometry // can be done scattered data reader CSV // split it in dof and quantities
      line = 0;

      // Open CSV file
      std::fstream geomfile(fileNameGeom_.c_str(),std::ios::in);

      if(!geomfile)
      {
        EXCEPTION("Scattered data geometry file '" << fileNameGeom_ << "' could not be opened!")
      }

      totLines = 0;
      while(std::getline(geomfile, row))
      {
        ++totLines;
      }

      // Variable for doubles values in a single line
      std::vector<double> vec;
      std::vector<double> coord(3);

      geomfile.clear();
      geomfile.seekg(0,std::ios::beg);
      // Iterate over lines in CSV file
      while(std::getline(geomfile, row))
      {
        line++;

        // Skip header lines.
        if(line <= skipLines_) continue;

        // Tokenize row into string tokens and convert them to doubles.
        Tokenizer tokens(row, boost::escaped_list_separator<char>('\\', ',', '\"'));
        if(vec.empty())
        {
          vec.resize(std::distance(tokens.begin(), tokens.end()));
        }

        Tokenizer::iterator tkIt(tokens.begin());

        for (UInt i=0; tkIt!=tokens.end(); ++tkIt, i++)
        {
          vec[i] = boost::lexical_cast<Double,std::string>(*tkIt);
        }

        std::map<UInt, UInt>::iterator dofIt, dofEnd;

        // Get coordinates from row.
        dofIt = dof2CoordColumn_.begin();
        dofEnd = dof2CoordColumn_.end();

        for( ; dofIt != dofEnd; dofIt++ )
        {
          coord[dofIt->first] = vec[dofIt->second];
        }

        coordinates_.push_back(coord);

      }

      geomfile.close();

    } // end if Step==1


    // Quantities // can be done scattered data reader CSV // split it in dof and quantities
    Double curTStep = this->mp_->Eval(mHandleStep_);
    // Iterator over a string vector
    StdVector< std::string >::iterator vfIt(valueFileName_.Begin());

    std::string valueFileName;

    for(UInt i=0;vfIt!=valueFileName_.End();++vfIt,++i){

      if(!(std::abs(valuetime_[i]-curTStep)<1e-10)) continue;
      typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
      std::string row;
      UInt line = 0;
      // Variable for doubles values in a single line
      std::vector<double> vec;
      std::vector<double> vecImag;
      std::vector<double> coord(3);

      // Open CSV file
      std::fstream datafile((*vfIt).c_str(),std::ios::in);

      if(!datafile)
      {
        EXCEPTION("Scattered data value file '" << (*vfIt) << "' could not be opened!")
      }

      while(std::getline(datafile, row))
      {
        line++;

        // Skip header lines.
        if(line <= skipLines_) continue;

        // Tokenize row into string tokens and convert them to doubles.
        Tokenizer tokens(row, boost::escaped_list_separator<char>('\\', ',', '\"'));
        if(vec.empty())
        {
          vec.resize(std::distance(tokens.begin(), tokens.end()));
        }
        if(vecImag.empty())
        {
          vecImag.resize(std::distance(tokens.begin(), tokens.end()));
        }

        Tokenizer::iterator tkIt(tokens.begin());

        for (UInt i=0; tkIt!=tokens.end(); ++tkIt, i++)
        {
          if((*tkIt).find("i")!= std::string::npos){
            Tokenizer tokComplex(*tkIt, boost::escaped_list_separator<char>('\\', 'i', '\"'));
            Tokenizer::iterator tkItComp(tokComplex.begin());

            std::string realPart = (*(tkItComp)).substr(0,(*(tkItComp)).size()-1);


            std::string signImag = (*(tkItComp)).substr((*(tkItComp)).size()-1,(*(tkItComp)).size());
            std::string imagPart = signImag + *(++tkItComp);

            vec[i] = boost::lexical_cast<Double,std::string>(realPart);
            vecImag[i] = boost::lexical_cast<Double,std::string>(imagPart);

          }else{
            vec[i] = boost::lexical_cast<Double,std::string>(*tkIt);
            vecImag[i] = 0.0;
          }
        }

        std::map<UInt, UInt>::iterator dofIt, dofEnd;

        // Get quantities from row
        std::set<std::string>::iterator qIt, qEnd;
        qIt = registeredQuantities_.begin();
        qEnd = registeredQuantities_.end();

        for( ; qIt != qEnd; ++qIt)
        {
          dofIt = qidDof2Column_[*qIt].begin();
          dofEnd = qidDof2Column_[*qIt].end();
          std::vector<double> qdofs(std::distance(dofIt, dofEnd));
          std::vector<double> qdofsImag(std::distance(dofIt, dofEnd));

          if(scatteredDataPerQuantity_.find(*qIt) == scatteredDataPerQuantity_.end()){
            std::vector< std::vector<double> > & ref = scatteredDataPerQuantity_[*qIt];
            std::vector< std::vector<double> > & refI = scatteredDataPerQuantityImag_[*qIt];
            ref.resize(totLines,std::vector<double>(qdofs.size()));
            refI.resize(totLines,std::vector<double>(qdofs.size()));
          }

          for( ; dofIt != dofEnd; dofIt++ )
          {
            qdofs[dofIt->first] = vec[dofIt->second];
            qdofsImag[dofIt->first] = vecImag[dofIt->second];
          }

          scatteredDataPerQuantity_[*qIt][line-1] = qdofs;
          scatteredDataPerQuantityImag_[*qIt][line-1] = qdofsImag;

        }

      }

      datafile.close();
      step_ = curStep;

    }




  }
}
