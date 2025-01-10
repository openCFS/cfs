#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/tokenizer.hpp>

#include "General/Exception.hh"

#include "ScatteredDataReaderCSV.hh"


namespace CoupledField 
{
  ScatteredDataReaderCSV::ScatteredDataReaderCSV(PtrParamNode& scatteredDataNode,
                                                 bool verbose) :
    ScatteredDataReader(scatteredDataNode, verbose)
  {
  }
  
  ScatteredDataReaderCSV::~ScatteredDataReaderCSV()
  {
  }

  void ScatteredDataReaderCSV::ParseParamNode() 
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

  void ScatteredDataReaderCSV::ReadData()
  {
    if(registeredQuantities_.empty())
    {
      return;
    }

    ParseParamNode();
    
    typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
    std::string row;
    UInt line = 0;

    // Open CSV file
    boost::filesystem::ifstream myfile(fileName_);

    if(!myfile)
    {
      EXCEPTION("Scattered data file '" << fileName_ << "' could not be opened!")
    }
    
    // Variable for doubles values in a single line
    std::vector<double> vec;
    std::vector<double> vecImag; // we only initialize it and set it to zero so that harmonic stuff works, although we just ignore it
    std::vector<double> coord(3);

    // Iterate over lines in CSV file
    while(std::getline(myfile, row))
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
        Double value;
        std::stringstream sstr;
        
        sstr << (*tkIt);
        sstr >> value;
        
        vec[i] = value;
        vecImag[i] = 0.0;
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

      // Get quantities from row
      std::set<std::string>::iterator qIt, qEnd;
      qIt = registeredQuantities_.begin();
      qEnd = registeredQuantities_.end();
      
      for( ; qIt != qEnd; qIt++ ) 
      {
        dofIt = qidDof2Column_[*qIt].begin();
        dofEnd = qidDof2Column_[*qIt].end();
        std::vector<double> qdofs(std::distance(dofIt, dofEnd));
        std::vector<double> qdofsImag(std::distance(dofIt, dofEnd));

        for( ; dofIt != dofEnd; dofIt++ ) 
        {
          qdofs[dofIt->first] = vec[dofIt->second];
          qdofsImag[dofIt->first] = vecImag[dofIt->second];
        }
        
        scatteredDataPerQuantity_[*qIt].push_back(qdofs);
        scatteredDataPerQuantityImag_[*qIt].push_back(qdofsImag);
      }
    }

    myfile.close();
    
  }
}
