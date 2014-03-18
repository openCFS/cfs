#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/tokenizer.hpp>

#include "General/Exception.hh"

#include "ScatteredDataReaderCSV.hh"


namespace CoupledField 
{

  ScatteredDataReaderCSV::~ScatteredDataReaderCSV()
  {
  }

  void ScatteredDataReaderCSV::Read(std::vector< std::vector<double> >& scatteredData)
  {
    typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
    std::string row;
    UInt line = 0;
    boost::filesystem::ifstream myfile(fileName_);
    
    while(std::getline(myfile, row))
    {
      line++;
      
      if(line <= skipLines_) continue;
      
      Tokenizer tokens(row, boost::escaped_list_separator<char>('\\', ',', '\"'));
      scatteredData.push_back(std::vector<double>());
      scatteredData.rbegin()->resize(std::distance(tokens.begin(), tokens.end()));
      std::vector<double>& vec = (*scatteredData.rbegin());
      
      Tokenizer::iterator tkIt(tokens.begin());
      
      for (UInt i=0; tkIt!=tokens.end(); ++tkIt, i++) 
      {
        Double value;
        std::stringstream sstr;
        
        sstr << (*tkIt);
        sstr >> value;
        
        vec[i] = value;
        
        //          scatteredData.rbegin()->push_back(value);
      }
    }

    myfile.close();
  }

}
