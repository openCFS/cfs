#include <boost/filesystem/fstream.hpp>
#include <boost/tokenizer.hpp>

#ifdef __MINGW64__
#include <intrin.h>
#endif

#include "CoefFunctionScatteredData.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField{

  template<class T>
  class CSVReader
  {
  public:
    CSVReader(std::istream& myfile, UInt skipLines) :
      myfile_(myfile),
      skipLines_(skipLines)
    {
    }
    
    ~CSVReader() 
    {
    }

  private:
    std::istream& myfile_;
    UInt skipLines_;

  public:

    std::istream& ReadCSV(std::vector<std::vector<T> >& data)
    {
      typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
      std::string row;
      UInt line = 0;
      
      while(std::getline(myfile_, row))
      {
        line++;

        if(line <= skipLines_) continue;
        
        Tokenizer tokens(row, boost::escaped_list_separator<char>('\\', ',', '\"'));
        data.push_back(std::vector<T>());
        data.rbegin()->resize(std::distance(tokens.begin(), tokens.end()));
        std::vector<T>& vec = (*data.rbegin());
        
        Tokenizer::iterator tkIt(tokens.begin());

        for (UInt i=0; tkIt!=tokens.end(); ++tkIt, i++) 
        {
          Double value;
          std::stringstream sstr;
          
          sstr << (*tkIt);
          sstr >> value;
          
          vec[i] = value;
          
          //          data.rbegin()->push_back(value);
        }
      }

      return myfile_;
    }
  };

  template<typename T, UInt DOFS>
  CoefFunctionScatteredData<T,DOFS>::CoefFunctionScatteredData(const std::string& fileName)
    : CoefFunction(),
      counter_(0)
  {
    dimType_ = VECTOR;

    boost::filesystem::ifstream* in = NULL;
    in = new boost::filesystem::ifstream(fileName);

    if(in) 
    {
      CSVReader<double> csvReader((*in), 1);
      csvReader.ReadCSV(scatteredData_);
      in->close();
      in = NULL;
    }

    std::list<Point> points;

    for(UInt i=0, n=scatteredData_.size(); i<n; i++)
    {      
      points.push_back(Point(scatteredData_[i][3],
                             scatteredData_[i][4],
                             scatteredData_[i][5],
                             scatteredData_[i][0],
                             scatteredData_[i][1],
                             scatteredData_[i][2]));
    }
    searchTree_.reset(new Tree(points.begin(), points.end()));
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::GetVector( Vector<T>& vec, 
                                                     const LocPointMapped& lpm ) {
    vec.Resize(DOFS);
    //    UInt size = scatteredData_.size();
    
    Vector<double> globPoint(DOFS);
    Vector<T> diff(DOFS);
    lpm.shapeMap->Local2Global(globPoint, lpm.lp);

    Point query(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    if(DOFS == 2) 
    {
      Point query2(globPoint[0], globPoint[1], 0.0, 0.0, 0.0, 0.0);
      query = query2;
    }
    else 
    {
      Point query3(globPoint[0], globPoint[1], globPoint[2], 0.0, 0.0, 0.0);
      query = query3;
    }

    K_neighbor_search search(*searchTree_.get(), query, 1);

#if 0
    // report the N nearest neighbors and their distance
    // This should sort all N points by increasing distance from origin
    for(K_neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
      std::cout << "neighbours: " << it->first.x() << " " << it->first.y() << " " << it->first.z() << " "<< std::sqrt(it->second) << std::endl << it->first.vx() << " " << it->first.vy() << " " << it->first.vz() << std::endl;
    }

    EXCEPTION("Global point: " << globPoint << "\nquery: " << query.x() << " " << query.y() << " " << query.z());
#endif
    

    K_neighbor_search::iterator it = search.begin();

    if(DOFS == 2) 
    {
      vec[0] = it->first.vx();
      vec[1] = it->first.vy();
    }
    else 
    {
      vec[0] = it->first.vx();
      vec[1] = it->first.vy();
      vec[2] = it->first.vz();
    }
  }

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionScatteredData<Double,2>;
  template class CoefFunctionScatteredData<Complex,2>;

  template class CoefFunctionScatteredData<Double,3>;
  template class CoefFunctionScatteredData<Complex,3>;
#endif
}// end of namespace

