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
  CoefFunctionScatteredData<T,DOFS>::CoefFunctionScatteredData(PtrParamNode& scatteredDataNode)
    : CoefFunction()
  {
    std::string fn = scatteredDataNode->Get("fileName")->As<std::string>();
          
    std::map<UInt, UInt> dof2CoordColumn;
    std::map<UInt, UInt> dof2ValueColumn;
    ParamNodeList coordList;
    coordList = scatteredDataNode->Get("coordinates")->GetList("comp");
    for(UInt i=0, n=coordList.GetSize(); i<n; i++) {
      std::string dof = coordList[i]->Get("dof")->As<std::string>();
      
      if( dof == "x" ) {
        dof2CoordColumn[0] = coordList[i]->Get("col")->As<UInt>();
      }
      if( dof == "y" ) {
        dof2CoordColumn[1] = coordList[i]->Get("col")->As<UInt>();
      }
      if( dof == "z" ) {
        dof2CoordColumn[2] = coordList[i]->Get("col")->As<UInt>();
      }
    }

    ParamNodeList valueList;
    valueList = scatteredDataNode->Get("values")->GetList("comp");
    for(UInt i=0, n=valueList.GetSize(); i<n; i++) {
      std::string dof = valueList[i]->Get("dof")->As<std::string>();
      
      if( dof == "x" ) {
        dof2ValueColumn[0] = valueList[i]->Get("col")->As<UInt>();
      }
      if( dof == "y" ) {
        dof2ValueColumn[1] = valueList[i]->Get("col")->As<UInt>();
      }
      if( dof == "z" ) {
        dof2ValueColumn[2] = valueList[i]->Get("col")->As<UInt>();
      }      
    }

    Double factor = scatteredDataNode->Get("factor")->As<Double>();

#ifdef USE_CGAL
    dimType_ = VECTOR;


    boost::filesystem::ifstream* in = NULL;
    in = new boost::filesystem::ifstream(fn);

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
    switch(dof2CoordColumn.size())
    {
      case 2:
        points.push_back(Point(scatteredData_[i][dof2CoordColumn[0]],
                               scatteredData_[i][dof2CoordColumn[1]],
                               0.0,
                               scatteredData_[i][dof2ValueColumn[0]] * factor,
                               scatteredData_[i][dof2ValueColumn[1]] * factor,
                               0.0));
        break;        
      case 3:
        points.push_back(Point(scatteredData_[i][dof2CoordColumn[0]],
                               scatteredData_[i][dof2CoordColumn[1]],
                               scatteredData_[i][dof2CoordColumn[2]],
                               scatteredData_[i][dof2ValueColumn[0]] * factor,
                               scatteredData_[i][dof2ValueColumn[1]] * factor,
                               scatteredData_[i][dof2ValueColumn[2]] * factor));
        break;        
    }
    
    }
    searchTree_.reset(new Tree(points.begin(), points.end()));
#else
    // We do not want to throw an exception from a constructor.
    WARN("CGAL is needed for nearest-neighbor mapping.")
#endif
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::GetVector( Vector<T>& vec, 
                                                     const LocPointMapped& lpm ) {
#ifdef USE_CGAL
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
#else
    EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif
  }

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionScatteredData<Double,2>;
  template class CoefFunctionScatteredData<Complex,2>;

  template class CoefFunctionScatteredData<Double,3>;
  template class CoefFunctionScatteredData<Complex,3>;
#endif
}// end of namespace

