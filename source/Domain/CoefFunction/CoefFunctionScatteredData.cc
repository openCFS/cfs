#include <boost/filesystem/fstream.hpp>
#include "boost/filesystem/operations.hpp"
#include <boost/tokenizer.hpp>

#ifdef __MINGW64__
#include <intrin.h>
#endif

#include "CoefFunctionScatteredData.hh"
#include "FeBasis/FeSpace.hh"

#include "def_use_ccmio.hh"
#include "DataInOut/ScatteredDataInOut/ScatteredDataReaderCSV.hh"
#ifdef USE_CCMIO
#include "DataInOut/ScatteredDataInOut/ScatteredDataReaderCCM.hh"
#endif


namespace CoupledField{

  template<typename T, UInt DOFS>
  CoefFunctionScatteredData<T,DOFS>::CoefFunctionScatteredData(PtrParamNode& scatteredDataNode)
    : CoefFunction(),
      fileName_(""),
      factor_(1.0)
  {
    dimType_ = VECTOR;
    dependType_ = CoefFunction::GENERAL;

    fileName_ = scatteredDataNode->Get("fileName")->As<std::string>();

    format_ = scatteredDataNode->Get("format")->As<std::string>();

    if(format_ == "csv") 
    {
      ParamNodeList coordList;
      coordList = scatteredDataNode->Get("coordinates")->GetList("comp");
      for(UInt i=0, n=coordList.GetSize(); i<n; i++) {
        std::string dof = coordList[i]->Get("dof")->As<std::string>();
        
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
      
      ParamNodeList valueList;
      valueList = scatteredDataNode->Get("values")->GetList("comp");
      for(UInt i=0, n=valueList.GetSize(); i<n; i++) {
        std::string dof = valueList[i]->Get("dof")->As<std::string>();
        
        if( dof == "x" ) {
          dof2ValueColumn_[0] = valueList[i]->Get("col")->As<UInt>();
        }
        if( dof == "y" ) {
          dof2ValueColumn_[1] = valueList[i]->Get("col")->As<UInt>();
        }
        if( dof == "z" ) {
          dof2ValueColumn_[2] = valueList[i]->Get("col")->As<UInt>();
        }      
      }
    } else 
    {
      dof2CoordColumn_[0] = 0;
      dof2CoordColumn_[1] = 1;
      dof2CoordColumn_[2] = 2;
      
      dof2ValueColumn_[0] = 3;
      dof2ValueColumn_[1] = 4;
      dof2ValueColumn_[2] = 5;
    }    
    
    factor_ = scatteredDataNode->Get("factor")->As<Double>();
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::Read() 
  {
    if(!boost::filesystem::exists(fileName_))
    {
      EXCEPTION("Scattered data file '" << fileName_ << "' does not exist!")
    }

#ifdef USE_CGAL
    ScatteredDataReaderPtr sdataReader;

    if(format_ == "csv") 
    {
      ScatteredDataReaderCSV* SCRCSV = new ScatteredDataReaderCSV(fileName_);
      SCRCSV->SetNumSkipLines(1);
      sdataReader.reset(SCRCSV);
    } else if (format_ == "ccm") {
#ifdef USE_CCMIO
      ScatteredDataReaderCCM* SCRCCM = new ScatteredDataReaderCCM(fileName_);
      std::vector<std::string> componentShortNames(3);
      componentShortNames[0] = "SU";
      componentShortNames[1] = "SV";
      componentShortNames[2] = "SW";
      SCRCCM->SetComponentShortNames(componentShortNames);
      sdataReader.reset(SCRCCM);
#else
      EXCEPTION("STAR-CCM+ files not supported! Compile with USE_CCMIO=ON.");
#endif
    } else {
      EXCEPTION("No format for scattered data file specified!");
    }      
    sdataReader->Read(scatteredData_);

    std::list<Point> points;

    for(UInt i=0, n=scatteredData_.size(); i<n; i++)
    {
    switch(dof2CoordColumn_.size())
    {
      case 2:
        points.push_back(Point(scatteredData_[i][dof2CoordColumn_[0]],
                               scatteredData_[i][dof2CoordColumn_[1]],
                               0.0,
                               scatteredData_[i][dof2ValueColumn_[0]] * factor_,
                               scatteredData_[i][dof2ValueColumn_[1]] * factor_,
                               0.0));
        break;        
      case 3:
        points.push_back(Point(scatteredData_[i][dof2CoordColumn_[0]],
                               scatteredData_[i][dof2CoordColumn_[1]],
                               scatteredData_[i][dof2CoordColumn_[2]],
                               scatteredData_[i][dof2ValueColumn_[0]] * factor_,
                               scatteredData_[i][dof2ValueColumn_[1]] * factor_,
                               scatteredData_[i][dof2ValueColumn_[2]] * factor_));
        break;        
    }
    
    }
    searchTree_.reset(new Tree(points.begin(), points.end()));
#else
    EXCEPTION("Switch on USE_CGAL for nearest-neighbor mapping!")
#endif
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::GetVector( Vector<T>& vec, 
                                                     const LocPointMapped& lpm ) {
    if(!scatteredData_.size())
    {
      Read();
    }    

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

    if(it == search.end())
    {
      EXCEPTION("Could not find a nearest neighbor for " << globPoint << "!");
    }
    

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

  template<typename T, UInt DOFS>
  std::string CoefFunctionScatteredData<T,DOFS>::ToString() const {
  std::stringstream out;
  out << "CoefFunctionScatteredData<" << typeid(T).name() << ", " << DOFS << ">";
  return out.str();
}


#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionScatteredData<Double,2>;
  template class CoefFunctionScatteredData<Complex,2>;

  template class CoefFunctionScatteredData<Double,3>;
  template class CoefFunctionScatteredData<Complex,3>;
#endif
}// end of namespace

