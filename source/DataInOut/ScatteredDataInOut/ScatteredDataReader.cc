#include <iostream>
#include <fstream>

#include "General/Exception.hh"
#include "ScatteredDataReader.hh"

#include "def_use_ccmio.hh"

#include "ScatteredDataReaderCSV.hh"
#ifdef USE_CCMIO
#include "ScatteredDataReaderCCM.hh"
#endif

using namespace std;

namespace CoupledField 
{

  std::map< std::string,
            boost::shared_ptr<ScatteredDataReader> > ScatteredDataReader::readers_;

  std::map<std::string, std::string > ScatteredDataReader::quantities2Readers_;

  void ScatteredDataReader::CreateReaders(PtrParamNode& scatteredDataNode)
  {
    // Check wether readers have already been created.
    if(!readers_.empty())
    {
      return;
    }

    std::set<std::string> readerIds;
    std::set<std::string> quantityIds;

    // Let's first check if we have unique reader ids.
    ParamNodeList scatteredNodes;
    scatteredNodes = scatteredDataNode->GetChildren();
    for(UInt i=0, n=scatteredNodes.GetSize(); i<n; i++) {
      std::string id = scatteredNodes[i]->Get("id")->As<std::string>();

      if(readerIds.find(id) != readerIds.end()) 
      {
        EXCEPTION("Id '" << id << "' for scattered data readers is not unique!");
      }
      else 
      {
        readerIds.insert(id);
      }
        
      ParamNodeList quantityNodes;
      quantityNodes = scatteredNodes[i]->GetList("quantity");
      // Let's first check if we have unique ids.
      for(UInt j=0, m=quantityNodes.GetSize(); j<m; j++) {
        std::string qid = quantityNodes[j]->Get("id")->As<std::string>();

        if(quantityIds.find(qid) != quantityIds.end()) 
        {
          EXCEPTION("Quantity id '" << qid << "' for scattered data is not unique!");
        }
        else 
        {
          quantityIds.insert(qid);
        }
      }    
    }

    // Now, that we  are sure we have  unique reader and quantity  ids, we can
    // actually create the readers.
    scatteredNodes = scatteredDataNode->GetChildren();
    for(UInt i=0, n=scatteredNodes.GetSize(); i<n; i++) {
      std::string id = scatteredNodes[i]->Get("id")->As<std::string>();
      std::string fn = scatteredNodes[i]->Get("fileName")->As<std::string>();
      std::string type = scatteredNodes[i]->GetName();

      boost::shared_ptr<ScatteredDataReader> reader;
      if(type == "csv") 
      {
        ScatteredDataReaderCSV* SCRCSV = new ScatteredDataReaderCSV(scatteredNodes[i]);
        SCRCSV->SetNumSkipLines(1);
        reader.reset(SCRCSV);
      } else if (type == "ccm") {
#ifdef USE_CCMIO
        ScatteredDataReaderCCM* SCRCCM = new ScatteredDataReaderCCM(scatteredNodes[i]);
        reader.reset(SCRCCM);
#else
        EXCEPTION("STAR-CCM+ files not supported! Compile with USE_CCMIO=ON.");
#endif
      } else {
        EXCEPTION("Unknown type '" << type << "' for scattered data file!");
      }
      readers_[id] = reader;

      ParamNodeList quantityNodes;
      quantityNodes = scatteredNodes[i]->GetList("quantity");
      for(UInt j=0, m=quantityNodes.GetSize(); j<m; j++) {
        std::string qid = quantityNodes[j]->Get("id")->As<std::string>();
        quantities2Readers_[qid] = id;
      }    
    }
  }

  void ScatteredDataReader::RegisterQuantity(const std::string& quantity)
  {
    if(readers_.empty()) 
    {
      EXCEPTION("Readers for scattered data have not been created yet.");
    }

    if(quantities2Readers_.find(quantity) == quantities2Readers_.end()) 
    {
      EXCEPTION("Desired quantity '" << quantity << "' is not available "
                << "through any of the readers for scattered data.");
    }
    
    readers_[quantities2Readers_[quantity]]->registeredQuantities_.insert(quantity);
  }

  void ScatteredDataReader::Read() 
  {
    std::map<std::string, boost::shared_ptr<ScatteredDataReader> >::iterator it, end;

    it = readers_.begin();
    end = readers_.end();

    for( ; it != end; it++ ) 
    {
      it->second->ReadData();
    }
  }

  void ScatteredDataReader::GetQuantity(const std::string& quantity,
                                        std::vector< std::vector<double> >& coordinates,
                                        std::vector< std::vector<double> >& scatteredData)
  {
    if(readers_.empty()) 
    {
      EXCEPTION("Readers for scattered data have not been created yet.");
    }

    if(quantities2Readers_.find(quantity) == quantities2Readers_.end()) 
    {
      EXCEPTION("Desired quantity '" << quantity << "' is not available "
                << "through any of the readers for scattered data.");
    }
    
    coordinates = readers_[quantities2Readers_[quantity]]->coordinates_;
    scatteredData = readers_[quantities2Readers_[quantity]]->scatteredDataPerQuantity_[quantity];
  }
}
